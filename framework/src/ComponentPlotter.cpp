/*!
 * @class ComponentPlotter
 *
 * A class for plotting PDF projections onto histograms
 *
 * @author Robert Currie rcurrie@cern.ch
 */

///	ROOT Headers
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TFolder.h"
#include "TTree.h"
#include "THStack.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TPad.h"
///	RapidFit Headers
#include "ComponentPlotter.h"
#include "StatisticsFunctions.h"
#include "EdStyle.h"
#include "StringProcessing.h"
#include "PhaseSpaceBoundary.h"
#include "ClassLookUp.h"
#include "ComponentRef.h"
#include "NormalisedSumPDF.h"
#include "ObservableDiscreteConstraint.h"
///	System Headers
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <cstdlib>
#include <algorithm>
#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

//Constructor with correct arguments
ComponentPlotter::ComponentPlotter( IPDF * NewPDF, IDataSet * NewDataSet, TString PDFStr, TFile* filename, string ObservableName, CompPlotter_config* config, int PDF_Num, const DebugClass* Debug ) :
	observableName( ObservableName ), weightName(), plotPDF( ClassLookUp::CopyPDF(NewPDF) ), plotData( NewDataSet ),
	pdfIntegrator( NULL ), weightsWereUsed(false), weight_norm(1.),
	discreteNames(), continuousNames(), full_boundary( new PhaseSpaceBoundary(*(NewDataSet->GetBoundary())) ), PlotFile( filename ),
	total_points( (config!=NULL)?config->PDF_points:128 ), data_binning( (config!=NULL)?config->data_bins:100 ), pdfStr( PDFStr ),
	logY( (config!=NULL)?config->logY:false ), logX( (config!=NULL)?config->logX:false ), this_config( config ), boundary_min( -99999 ), boundary_max( -99999 ), step_size( -99999 ),
	onlyZero( (config!=NULL)?config->OnlyZero:false ), combination_integral(vector<double>()), ratioOfIntegrals(1,1.), wanted_weights(), format(),
	data_subsets(), allCombinations(), combinationWeights(), combinationDescriptions(), observableValues(), binned_data(), total_components(),
	chi2(), N(), allPullData(), debug(NULL), PDFNum(PDF_Num)
{
	TH1::SetDefaultSumw2(true);

	//pdfIntegrator->SetNumThreads( config!=NULL?config->numThreads:4 );

	if( Debug != NULL ) debug = new DebugClass(*Debug);
	else debug = new DebugClass(false);
	plotPDF->TurnCachingOff();
	plotPDF->SetComponentStatus( true );

	//plotPDF->ChangePhaseSpace( full_boundary );

	//plotPDF->GetPhysicsParameters()->Print();

	pdfIntegrator = new RapidFitIntegrator( plotPDF );

	pdfIntegrator->SetPDF( plotPDF );
	pdfIntegrator->SetDebug(debug);

	if( pdfIntegrator->GetUseGSLIntegrator() ) cout << "Using GSL for projections" << endl;

	RapidFitIntegratorConfig* projectionIntegratorConfig = NULL;

	if( config != NULL ) projectionIntegratorConfig = config->integratorConfig;

	if( projectionIntegratorConfig != NULL )
	{
		RapidFitIntegratorConfig* thisIntConfig = new RapidFitIntegratorConfig( *projectionIntegratorConfig );
		pdfIntegrator->SetUpIntegrator( thisIntConfig );
		plotPDF->SetUpIntegrator( thisIntConfig );
		delete thisIntConfig;
	}

	vector<string> disc_contr = full_boundary->GetDiscreteNames();
	vector<string> pdf_const = plotPDF->GetPrototypeDataPoint();
	for( unsigned int i=0; i< disc_contr.size(); ++i )
	{
		if( StringProcessing::VectorContains( &pdf_const, &(disc_contr[i]) ) == -1 )
		{
			full_boundary->RemoveConstraint( disc_contr[i] );
		}
	}
	plotPDF->ChangePhaseSpace( full_boundary );
	NewDataSet->SetBoundary( full_boundary );

	//////////////////////////////////////////////////////////////////
	//Make the histogram of this observable
	format = new EdStyle();
	format->SetStyle();
	//gStyle->SetMarkerStyle(15);
	//gStyle->SetMarkerSize(Size_t(0.8));

	gStyle->SetPadLeftMargin( (Float_t)0.15 );
	gStyle->SetTitleOffset((Float_t)0.9,"Y");

	boundary_min = full_boundary->GetConstraint( observableName )->GetMinimum();
	boundary_max = full_boundary->GetConstraint( observableName )->GetMaximum();

	step_size = ( boundary_max - boundary_min ) / (double)( total_points - 1 );

	//Work out what to plot
	vector<DataPoint*> allCombinations_input;// = full_boundary->GetDiscreteCombinations();

	if( !config->plotAllCombinations )
	{
		vector<string> DiscreteObs = full_boundary->GetDiscreteNames();
		for( unsigned int i=0; i< DiscreteObs.size(); ++i )
		{
			if( full_boundary->GetConstraint( DiscreteObs[i] )->GetValues().size() > 1 )
			{
				full_boundary->RemoveConstraint( DiscreteObs[i] );
				IConstraint* thisConstraint = (IConstraint*) new ObservableDiscreteConstraint( DiscreteObs[i], vector<double>(1,config->defaultCombinationValue), " ", "" );
				full_boundary->AddConstraint( DiscreteObs[i], thisConstraint );
				delete thisConstraint;
			}
		}
		plotPDF->ChangePhaseSpace( full_boundary );
	}

	allCombinations_input = full_boundary->GetDiscreteCombinations();

	for( unsigned int i=0; i< allCombinations_input.size(); ++i )
	{
		allCombinations.push_back( new DataPoint( *(allCombinations_input[i]) ) );
	}

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "ComponentPlotter" ) )
		{
			cout << "ComponentPlotter:: Setting up Combination Descriptions based on XML" << endl;
		}
	}

	this->SetupCombinationDescriptions();

	if( config != NULL )
	{
		if( !config->combination_names.empty() ) combinationDescriptions = config->combination_names;
	}

	for( unsigned int i=0; i< allCombinations.size(); ++i )
	{
		if( debug != NULL )
		{
			if( debug->DebugThisClass( "ComponentPlotter" ) )
			{
				cout << "ComponentPlotter:: Calculating Test integrals:\t" << i << endl;
			}
		}
		pdfIntegrator->ForceTestStatus( false );
		allCombinations[i]->SetPhaseSpaceBoundary( full_boundary );
		double thisIntegral = 0.;
		try
		{
			allCombinations[i]->SetPhaseSpaceBoundary( full_boundary );
			thisIntegral = pdfIntegrator->NumericallyIntegrateDataPoint( allCombinations[i], full_boundary, plotPDF->GetDoNotIntegrateList() );
		}
		catch(...)
		{
			thisIntegral = 1.;
			cout << endl << "CANNOT PROPERLY NORMALISE WHOLE PDF, THIS WILL LEAD TO NORMALISATION ISSUES OVER THE WHOLE PDF" << endl << endl;
		}
                if( debug != NULL )
                {
                        if( debug->DebugThisClass( "ComponentPlotter" ) )
                        {
				cout << "ComponentPlotter:: Finished Inetgral" << endl;
			}
		}

		combination_integral.push_back( thisIntegral );

		if( plotPDF->GetNumericalNormalisation() == true ) ratioOfIntegrals.push_back( 1. );
		else
		{
			if( config->ScaleNumerical ) ratioOfIntegrals.push_back( pdfIntegrator->GetRatioOfIntegrals() );
			else ratioOfIntegrals.push_back( 1./pdfIntegrator->GetRatioOfIntegrals() );
		}
	}

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "ComponentPlotter" ) )
		{
			cout << "ComponentPlotter:: Calculated the ratio of Integrals (analytic/numeric) to be:" << endl;
			for( unsigned int i=0; i< ratioOfIntegrals.size(); ++i )
			{
				cout << ratioOfIntegrals[i] << endl;
			}
		}
	}
	if( ratioOfIntegrals.size() == 0 || ratioOfIntegrals.empty() ) ratioOfIntegrals =vector<double> ( 1, 1. );

	vector<double> minimum, maximum;
	vector<int> binNumber;

	(void) minimum; (void) maximum; (void) binNumber;

	if( config != NULL )
	{
		if( config->component_names.empty() )
		{
			if( plotPDF->GetName() == plotPDF->GetLabel() )
			{
				config->component_names.push_back( plotPDF->GetName() );
			}
			else
			{
				config->component_names.push_back( "Total" );
			}
			vector<string> pdfComponents = plotPDF->PDFComponents();
			pdfComponents = StringProcessing::MoveElementToStart( pdfComponents, "0" );
			for( unsigned int i=0; i< pdfComponents.size(); ++i )
			{
				ComponentRef* thisRef = new ComponentRef( pdfComponents[i], observableName );
				if( pdfComponents[i] != "0" ) config->component_names.push_back( plotPDF->GetComponentName( thisRef ) );
				delete thisRef;
			}
		}
	}

	TH1::SetDefaultSumw2( true );

	pdfIntegrator->ForceTestStatus( true );
}

//Destructor
ComponentPlotter::~ComponentPlotter()
{
	if( plotPDF != NULL )	delete plotPDF;
	if( pdfIntegrator != NULL ) delete pdfIntegrator;
	while( !allCombinations.empty() )
	{
		if( allCombinations.back() != NULL ) delete allCombinations.back();
		allCombinations.pop_back();
	}
	delete format;
	if( debug != NULL ) delete debug;
}

//Create a root file containing a projection plot over one observable
void ComponentPlotter::ProjectObservable()
{
	vector<string> doNotIntegrate = plotPDF->GetDoNotIntegrateList();

	plotPDF->UnsetCache();
	//Check the observable can be plotted
	bool continuous = !( plotData->GetBoundary()->GetConstraint( observableName )->IsDiscrete() );
	bool doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrate, &(observableName) ) == -1 );
	if( continuous && doIntegrate )
	{
		//Do the projecting of the pdfs
		cout << "Projecting " << observableName << endl;
		this->GenerateProjectionData();
	}
	else
	{
		cerr << "CANNOT PROJECT: " << observableName << endl;
		return;
	}
}

//
//	Returns a vector of arrays of dimention		AllCombinations (+ 0 combination) * total_points
//	This data corresponds to the X axis of the projected functions from the PDF
//
//	This is a simple loop of adding numbers onto a running total and recording the running total between minima and maxima in a vector
//
vector<vector<double>* >* ComponentPlotter::MakeXProjectionData( unsigned int num_combinations )
{
	vector<vector<double>* >* new_dataarray = new vector<vector<double>* >();

	unsigned int total = num_combinations;
	if( total > 1 ) ++total;

	for (unsigned int combinationIndex = 0; combinationIndex < total; ++combinationIndex )
	{
		vector<double>* observableValueArray = new vector<double>();
		observableValueArray->resize((unsigned)total_points);

		for(int pointIndex = 0; pointIndex < total_points; ++pointIndex )
		{
			//	Start at minimum of observable in phase-space
			//
			//	move a step equal to the size of the interval that you should take for this subset of data
			//
			//	due to the fact that the stepsize corresponding to the tag=-1 (combinationIndex=0)
			//	may be different ot the stepsize corresponding to the tag=1 (combinationIndex=1)
			//
			(*observableValueArray)[ unsigned(pointIndex) ] = boundary_min + ( step_size * (pointIndex) );
		}
		new_dataarray->push_back( observableValueArray );
	}

	return new_dataarray;
}

//
//      This returns the Y component of the projections for all combinations and for all components for each projection
//
//      This calls ComponentPlotter::ProjectObservableComponent for each combination
//
//      Combination 0 is constricted as the total of the individual components once all components for all sub component have been calculated
//
vector<vector<double>* >* ComponentPlotter::MakeYProjectionData( string component_str )
{
	vector<vector<double>* >* new_dataarray = new vector<vector<double>* >();

	//	Loop over all discrete combinations for this PDF configuration
	for( unsigned int combinationIndex = 0; combinationIndex < allCombinations.size(); ++combinationIndex )
	{
		cout << "Constructing PDF Integral of: " << observableName << " Combination: " << combinationIndex+1 << " of " << allCombinations.size() << ". For component: " << component_str <<" .\t";

		allCombinations[combinationIndex]->SetPhaseSpaceBoundary( full_boundary );

		//cout << "Combination: " << combinationIndex+1 << " is defined as:" << endl; allCombinations[combinationIndex]->Print();

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Get the raw probability distribution for this component from the PDF

		//Calculate the projection for this combination
		vector<double>* projectionValueVector =

			//	DataPoint with all information on this configuration,
			this->ProjectObservableComponent( allCombinations[combinationIndex],

					//	minimum of range in Observable,	num of points in range,	plot_interval,	component of interest
					observableName, boundary_min, total_points, step_size, component_str );

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		//Normalise the PDF such that the final decision of bin number is the only variable missing.

		//Update the data average values, and make the projection graph arrays

		vector<double>* projectionValueArray = projectionValueVector;

		if( debug != NULL )
		{
			if( debug->DebugThisClass( "ComponentPlotter" ) )
			{
				cout << "ComponentPlotter:: Normalising the PDF component: " << combinationIndex << endl;
				cout << "ComponentPlotter:: Normalisation: " << this->PDF2DataNormalisation( combinationIndex );
			}
		}

		double PDFNormalisation = this->PDF2DataNormalisation( combinationIndex );


		//	Perform Normalisation
		for( unsigned int pointIndex = 0; pointIndex < (unsigned)total_points; ++pointIndex )
		{
			(*projectionValueArray)[ pointIndex ] *= PDFNormalisation;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		cout << "Finished Combination: " << combinationIndex+1 << " of " << allCombinations.size() <<"."<< endl;

		new_dataarray->push_back( projectionValueArray );
	}

	new_dataarray = this->GenerateComponentZero( new_dataarray );

	return new_dataarray;
}

//	When we have more than 1 discrete component we need to create component 0 which contains the total PDF result at this coordinate
vector<vector<double>* >* ComponentPlotter::GenerateComponentZero( vector<vector<double>* >* new_dataarray )
{
	//      If we are asked to look at more than 1 component with this PDF we need to return each component indepedently as well as a total of all of them
	if( new_dataarray->size() > 1 )
	{
		//      Construct 0th combination and fill it as the sum of all sub-combinations of the dataset
		vector<double>* zero_component = new vector<double>();
		zero_component->resize((unsigned)total_points);
		for( unsigned int i=0; i< (unsigned)total_points; ++i )
		{
			//	Zeroth component is defined as the total of each of the sub components at this point
			(*zero_component)[i]=0;
			for( unsigned int j=0; j< new_dataarray->size(); ++j )
			{
				(*zero_component)[i]+= (*(*new_dataarray)[j])[i];
			}
		}

		//      Construct a new object to store 0th and all additional components to be returned
		vector<vector<double>* >* full_array = new vector<vector<double>* >();
		full_array->push_back( zero_component );
		for( unsigned int i=0; i< new_dataarray->size(); ++i )
		{
			full_array->push_back( (*new_dataarray)[i] );
		}
		//      Replace object that is to be returned with this one :D
		if( new_dataarray != NULL ) delete new_dataarray;
		new_dataarray = full_array;
	}

	return new_dataarray;
}

//	Generate Projection Data for a given Observable
//
//	When this function starts it collects statistical information and infromation on the boundary of the observable
//
//	From here it loops over all components and combinations for this PDFWithData object
//
//	For configurations with more than 1 combination I create combination 0 at the end to contain the total sum of all of the unique combination datasets and projections
//
void ComponentPlotter::GenerateProjectionData()
{
	if( debug != NULL )
	{
		if( debug->DebugThisClass( "ComponentPlotter" ) )
		{
			cout << "ComponentPlotter: Counting Number of Events " << endl;
		}
	}

	//	Now we have:
	//			plot phase-space
	//	Lets plot the projection

	//	Check if the 0'th component has been defined and add it if it's missing.
	vector<string> PDF_Components;
	if( onlyZero == true )
	{
		PDF_Components = vector<string>(1,"0");
	}
	else
	{
		PDF_Components = plotPDF->PDFComponents();
		if( StringProcessing::VectorContains( PDF_Components, string("0") ) == -1 )
		{
			vector<string> temp_PDF_Components( 1, "0" );
			for(vector<string>::iterator comb_i = PDF_Components.begin(); comb_i != PDF_Components.end(); ++comb_i )
			{
				temp_PDF_Components.push_back( *comb_i );
			}
			PDF_Components = temp_PDF_Components;
		}
		PDF_Components = StringProcessing::MoveElementToStart( PDF_Components, "0" );
	}


	cout << endl << "Components: " << PDF_Components.size() << endl;



	//	Yes I know these are pointers to templates, but this is faster and was quicker than coding up a struct to contain the information easily
	//	I apologise for any difficult you may have reading it, but this is easier to debug (imo)
	vector<vector<vector<double>* >* >* X_values = new vector<vector<vector<double>* >* >();
	vector<vector<vector<double>* >* >* Y_values = new vector<vector<vector<double>* >* >();



	//	Loop Over ALL components that are provided by this PDF
	for( unsigned int i=0; i< PDF_Components.size(); ++i )
	{
		cout << "\n\t\tCOMPONENT: " << i+1 << " of: " << PDF_Components.size() << "\t\tPDF: "<< plotPDF->GetLabel() << endl <<endl;

		//	Generate the Y values of the projection plot
		//								component_of_interest
		vector<vector<double>* >* Y_values_local = MakeYProjectionData( PDF_Components[i] );

		//	Generate the X values of the projection plot
		vector<vector<double>* >* X_values_local = MakeXProjectionData( (unsigned)allCombinations.size() );

		//	Store the Projection Data for this Component
		X_values->push_back( X_values_local );	Y_values->push_back( Y_values_local );
	}

	//	Write the output to the output file
	TDirectory* here = gDirectory;

	TString PDFDesc; PDFDesc+=PDFNum;
	vector<string> combDescs=combinationDescriptions;
	for( unsigned int i=0; i< combDescs.size(); ++i )
	{
		combDescs[i].append("_PDF_");
		combDescs[i].append(PDFDesc.Data());
	}

	cout << endl << "Writing Combinations:" << endl;
	for( unsigned int i=0; i< combinationDescriptions.size(); ++i ) cout << combinationDescriptions[i] << endl;
	cout << endl;

	this->WriteOutput( X_values, Y_values, combinationDescriptions );
	here->cd();

	for( unsigned int i=0; i< X_values->size(); ++i )
	{
		for( unsigned int j=0; j< (*X_values)[i]->size(); ++j )
		{
			delete (*(*X_values)[i])[j];
		}
		delete (*X_values)[i];
	}
	delete X_values;

	for( unsigned int i=0; i< Y_values->size(); ++i )
	{
		for( unsigned int j=0; j< (*Y_values)[i]->size(); ++j )
		{
			delete (*(*Y_values)[i])[j];
		}
		delete (*Y_values)[i];
	}
	delete Y_values;

	return;
}

//	This routine bins the data for a requested combinationNumber into a TH1, the number of bins is determined through the data_binning int
//
TH1* ComponentPlotter::FormatData( unsigned int combinationNumber )
{
	vector<DataPoint*> wanted_points = data_subsets[combinationNumber];
	TString component_num_str;component_num_str+=combinationNumber;
	TH1* returnable = new TH1D( "Data_For_Component_"+component_num_str, "Data_For_Component_"+component_num_str, data_binning, boundary_min, boundary_max );
	vector<double> wanted_data, this_wanted_weights;
	ObservableRef* new_ref = new ObservableRef( observableName );

	for( vector<DataPoint*>::iterator point_i = wanted_points.begin(); point_i != wanted_points.end(); ++point_i )
	{
		wanted_data.push_back( (*point_i)->GetObservable( *new_ref )->GetValue() );
	}
	delete new_ref;

	if( wanted_weights.empty() && weightsWereUsed ) wanted_weights = vector<double>( wanted_data.size(), 1. );	//	'Should' never occur... but better than a segfault

	this_wanted_weights = wanted_weights;

	//	Fill a histogram with the weighted data, with weights normalised to 1
	//	Think of this as multiplying by a very clever version of the number 1 ;)
	//	This gives the correct normalisation
	vector<double>::iterator point_i = wanted_data.begin();
	vector<double>::iterator weight_i = this_wanted_weights.begin();
	if( weightsWereUsed )
	{
		for( ; point_i != wanted_data.end(); ++point_i, ++weight_i )
		{
			returnable->Fill( *point_i, *weight_i );
		}
	}
	else
	{
		for( ; point_i != wanted_data.end(); ++point_i )
		{
			returnable->Fill( *point_i, 1. );
		}
	}

	return returnable;
}

//	This routine plots all of the data for all combinations and stores the data int TTree objects in a .root file
void ComponentPlotter::WriteOutput( vector<vector<vector<double>* >* >* X_values, vector<vector<vector<double>* >* >* Y_values, vector<string> CombinationDescriptions )
{
	//	Make sure we start in the correct place
	PlotFile->cd();

	//	Seperate directories per PDF in the XML
	if( gDirectory->GetDirectory( pdfStr ) == 0 )	gDirectory->mkdir( pdfStr );
	gDirectory->cd( pdfStr );

	//	Seperate directories per observable in the XML
	if( gDirectory->GetDirectory( observableName.c_str() ) == 0 )	gDirectory->mkdir( observableName.c_str() );
	gDirectory->cd( observableName.c_str() );

	//	Pointer to here
	TDirectory* PlotDirectory = gDirectory;

	//	Loop over all components
	for( unsigned int componentIndex=0; componentIndex < X_values->size(); ++componentIndex )
	{
		TString componentName("Component_");componentName+=componentIndex;
		if( gDirectory->GetDirectory( componentName ) == 0 )	gDirectory->mkdir( componentName );
		gDirectory->cd( componentName );
		TDirectory* componentDir = gDirectory;

		//	Loop over all combinations for this component
		for( unsigned int combinationIndex=0; combinationIndex < (*X_values)[componentIndex]->size(); ++combinationIndex )
		{
			TString combinationName("Combination_");combinationName+=combinationIndex;
			if( gDirectory->GetDirectory( combinationName ) == 0 )	gDirectory->mkdir( combinationName );
			gDirectory->cd( combinationName );

			//	Save all of the X and Y data in the TTree
			TString TTree_name( "TTree_" ); TTree_name+=componentIndex; TTree_name.Append("_"); TTree_name+=combinationIndex;
			TString TTree_title( TTree_name );
			TTree_name.Append("_");TTree_name+=plotPDF->GetRandomFunction()->Rndm();

			string TTreeName( TTree_name.Data() );
			replace( TTreeName.begin(), TTreeName.end(), '.', '_' );

			TTree* this_data = new TTree( TTreeName.c_str(), TTree_title );

			vector<double>* unnormalised = new vector<double>();
			vector<double>* normalisation = new vector<double>();

			double Normalisation = 0.;

			if( (*X_values)[componentIndex]->size() > 1 )
			{
				if( combinationIndex == 0 )
				{
					for( unsigned int i=0; i< (*X_values)[componentIndex]->size()-1; ++i )
					{
						Normalisation += this->PDF2DataNormalisation( i );
					}
				}
				else
				{
					Normalisation = this->PDF2DataNormalisation( combinationIndex-1 );
				}
			}
			else
			{
				Normalisation = this->PDF2DataNormalisation( 0 );
			}

			for( vector<double>::iterator num_i = (*(*Y_values)[componentIndex])[combinationIndex]->begin(); num_i != (*(*Y_values)[componentIndex])[combinationIndex]->end(); ++num_i )
			{
				unnormalised->push_back( (*num_i) / Normalisation );
				normalisation->push_back( Normalisation );
			}

			this->WriteBranch( this_data, "X_data", (*(*X_values)[componentIndex])[combinationIndex] );
			this->WriteBranch( this_data, "Y_data", (*(*Y_values)[componentIndex])[combinationIndex] );
			this->WriteBranch( this_data, "Y_data_unNormalised", unnormalised );
			this->WriteBranch( this_data, "Y_data_Normalisation", normalisation );

			this_data->Write( "", TObject::kOverwrite );

			delete unnormalised;
			delete normalisation;

			TString Data_Name("Raw_Data_"); Data_Name+=componentIndex; Data_Name.Append("_");Data_Name+=combinationIndex;
			TString Data_Title( Data_Name );
			Data_Name.Append("_"); Data_Name+=plotPDF->GetRandomFunction()->Rndm();

			TTree* raw_data = new TTree( Data_Name, Data_Title );

			vector<double>* real_raw_data = new vector<double>();

			ObservableRef tempref( observableName );

			for( vector<DataPoint*>::iterator point_i = data_subsets[combinationIndex].begin(); point_i != data_subsets[combinationIndex].end(); ++point_i )
			{
				real_raw_data->push_back( (*point_i)->GetObservable( tempref )->GetValue() );
			}

			this->WriteBranch( raw_data, "Value", real_raw_data );

			raw_data->Write( "", TObject::kOverwrite );

			delete real_raw_data;

			//	Bin the data for this combination
			TH1* data_plot = this->FormatData( combinationIndex );
			TString ext("_"); ext+=componentIndex; ext.Append("_"); ext+=combinationIndex;
			ext.Append("_"); ext+=plotPDF->GetRandomFunction()->Rndm();
			string cleanExt( ext.Data() );
			replace( cleanExt.begin(), cleanExt.end(), '.', '_' );
			data_plot->SetName( data_plot->GetName() + TString(cleanExt.c_str()) );
			data_plot->SetTitle( "" );

			//	Plot the component on this combination and save the TCanvas
			TString Graph_Name("TGraph_");Graph_Name+=componentIndex;
			Graph_Name.Append("_");Graph_Name+=combinationIndex;
			Graph_Name.Append("_");Graph_Name+=plotPDF->GetRandomFunction()->Rndm();

			string graphCleanName( Graph_Name.Data() );
			replace( graphCleanName.begin(), graphCleanName.end(), '.', '_' );

			TGraph* data_graph = new TGraph( total_points,  &((*(*(*X_values)[componentIndex])[combinationIndex])[0]),  &((*(*(*Y_values)[componentIndex])[combinationIndex])[0]) );
			data_graph->SetTitle( "" ); data_graph->SetName( graphCleanName.c_str() );

			TString Canvas_Name("TCanvas_");Canvas_Name+=componentIndex;
			Canvas_Name.Append("_");Canvas_Name+=combinationIndex;
			Canvas_Name.Append("_");Canvas_Name+=plotPDF->GetRandomFunction()->Rndm();

			string CanvasCleanName( Canvas_Name.Data() );
			replace( CanvasCleanName.begin(), CanvasCleanName.end(), '.', '_' );

			TCanvas* c1 = EdStyle::RapidFitCanvas( CanvasCleanName.c_str(), "" );

			data_plot->Draw();
			data_graph->Draw("PC9 SAME");
			c1->Update();
			double Y_min = -99999.;
			double Y_max = -99999.;
			if( this_config != NULL )
			{
				Y_min = this_config->ymin;
				Y_max = this_config->ymax;
			}
			c1->Update();
			if( Y_min <= -99999. ) Y_min = data_graph->GetYaxis()->GetXmin();
			if( Y_max <= -99999. ) Y_max = data_graph->GetYaxis()->GetXmax();
			data_graph->GetYaxis()->SetRangeUser( Y_min, Y_max );
			c1->Update();
			c1->Write();

			componentDir->cd();
		}

		PlotDirectory->cd();
	}

	cout << "Data Written to File, Making Graphs" << endl;

	//	Starting at the top of the file again
	PlotFile->cd();

	gDirectory->cd( pdfStr );

	if( gDirectory->GetDirectory( "overlay_graphs" ) == 0 )	gDirectory->mkdir( "overlay_graphs" );
	gDirectory->cd( "overlay_graphs" );

	//	Overlay all components for each combination
	//
	//	This is painful as your looping over the upper 2D of a 3D vector inside out...
	//	Sorry but I'm not rewriting the whole class or inverting the vector of vectors as we will just get more mistakes
	//
	//	This definitely works and isn't quite so bad as I feares but knowing which component we are using is a pain
	//
	for( unsigned int combinationIndex=0; combinationIndex < (*X_values)[0]->size(); ++combinationIndex )
	{
		TString combinationIndexstr;combinationIndexstr+=combinationIndex;

		vector<TGraph*> these_components;

		//	Objects must have unqiue names, even though they exist in different memory locations, THANK YOU ROOT GARBAGE COLLECTION!!! (this is the singally worst idea ever to grace c++!)
		for( unsigned int componentIndex=0; componentIndex < X_values->size(); ++componentIndex )
		{
			TGraph* data_graph = new TGraph( total_points,  &((*(*(*X_values)[componentIndex])[combinationIndex])[0]),  &((*(*(*Y_values)[componentIndex])[combinationIndex])[0]) );
			TString data_name("data_"+combinationIndexstr+"_");data_name+=componentIndex;
			data_name.Append("_"); data_name+=plotPDF->GetRandomFunction()->Rndm();

			string dataCleanName( data_name );
			replace( dataCleanName.begin(), dataCleanName.end(), '.', '_' );

			data_graph->SetName( dataCleanName.c_str() ); data_graph->SetTitle( dataCleanName.c_str() );
			data_graph->SetLineColor( (Color_t)(componentIndex+1) );
			data_graph->SetLineWidth( (Width_t)3 );
			data_graph->SetMarkerColor( (Color_t)(componentIndex+1) );
			these_components.push_back( data_graph );
		}
		total_components.push_back( these_components );

		TH1* data_plot = this->FormatData( combinationIndex );
		TString ext("_"); ext+=combinationIndex; ext.Append("_"); ext+=plotPDF->GetRandomFunction()->Rndm();
		string cleanExt( ext.Data() );
		replace( cleanExt.begin(), cleanExt.end(), '.', '_' );
		data_plot->SetName( data_plot->GetName() + TString(cleanExt.c_str()) );
		data_plot->SetTitle( "" );

		//	Object to hold the binned data for the life of this plot
		binned_data.push_back( new TGraphErrors( data_plot ) );

		binned_data.back()->Write("",TObject::kOverwrite);

		string desc;

		if( combinationIndex> 0 ) desc = CombinationDescriptions[combinationIndex-1];
		if( combinationIndex == 0 )
		{
			desc = "All_Data";
			/*
			   for( vector<string>::iterator desc_i = CombinationDescriptions.begin(); desc_i != CombinationDescriptions.end(); ++desc_i )
			   {
			   desc.append( *desc_i );
			   }
			   */
		}

		TString PDFDesc;PDFDesc+=PDFNum;

		desc.append("_PDF_"); desc.append(PDFDesc.Data());

		//	For the moment haven't decided if I should pass the global config to ALL sub plots, I shall get user input on this
		CompPlotter_config* temp = new CompPlotter_config( *this_config );
		temp->logY = logY;
		temp->logX = logX;

		//	Static function so has to be told everything about what you want to plot!
		this->OutputPlot( binned_data.back(), these_components, observableName, desc, plotData->GetBoundary(), plotPDF->GetRandomFunction(), temp, debug );

		if( this_config != NULL )
		{
			if( combinationIndex == 0 && this_config->CalcChi2 == true )
			{
				cout << endl;
				cout << "Calculating Chi^2:" << endl;

				TF1* fitting_function = new TF1( "total_PDF", this, boundary_min, boundary_max, 1, "" );        //      I REFUSE to pass the class name
				//chi2 = binned_data.back()->Chisquare( fitting_function );

				//double chi2_2 = 0.;
				for( unsigned int i=0; i< (unsigned) binned_data[0]->GetN(); ++i )
				{
					double bin_center = binned_data[combinationIndex]->GetX()[i];
					double this_bin = fitting_function->Eval( bin_center );
					double actual_bin = binned_data[combinationIndex]->GetY()[i];
					double this_chi = actual_bin-this_bin;
					this_chi/=binned_data[combinationIndex]->GetEY()[i];//sqrt(fabs(actual_bin));
					//cout << endl << actual_bin << "\t" << this_bin << endl;
					this_chi*=this_chi;
					if( this_chi < 1./0. ) chi2+=this_chi;
				}

				//chi2 = chi2_2;

				//	dof = num Populated Bins - ndof in PDF
				N =  binned_data.back()->GetXaxis()->GetNbins();//plotData->GetDataNumber(NULL);
				for( unsigned int i=0; i< (unsigned) binned_data[0]->GetN(); ++i ) if( fabs(binned_data[combinationIndex]->GetY()[i]) <= 0 ) --N;

				double n = (double) plotPDF->GetPhysicsParameters()->GetAllFloatNames().size();

				//cout << endl << chi2 << "\t" << chi2_2 << "\t" << N << "\t" << n << endl;

				double denominator = (double)(N - n - 1. );
				cout << endl << "Chi^2/ndof :\t" << setprecision(10) << chi2 / denominator << endl;

				/*
				   TString graphName("Chi2Graph");
				   TString graphTitle(graphName);
				   graphName+=plotPDF->GetRandomFunction()->Rndm();
				   TGraph* tf1_graph = new TGraph( fitting_function );
				   tf1_graph->SetTitle( graphTitle );
				   tf1_graph->SetName( graphName );

				   TString Chi2Title("Chi2Title");
				   TString Chi2Name(Chi2Title); Chi2Name+=plotPDF->GetRandomFunction()->Rndm();
				   TCanvas* Chi2test = EdStyle::RapidFitCanvas( Chi2Name, Chi2Title );
				   binned_data.back()->Draw("AP");
				   tf1_graph->Draw("C");
				   Chi2test->Print(Chi2Name+".pdf");
				   Chi2test->Write("",TObject::kOverwrite);
				   */
			}

			if( combinationIndex == 0 && this_config->DrawPull == true )
			{
				cout << endl;
				cout << "Calculating Pull Values" << endl;

				TString TF1_Name( "total_PDF_" ); TF1_Name.Append( desc );
				TF1* PDF_function = new TF1( TF1_Name, this, boundary_min, boundary_max, 1, "" );        //      I REFUSE to pass the class name

				if( allPullData.size() != (unsigned)binned_data[combinationIndex]->GetN() ) allPullData.resize( (unsigned)binned_data[combinationIndex]->GetN(), 0. );


				//cout << desc << endl;
				for( unsigned int i=0; i< (unsigned)binned_data[combinationIndex]->GetN(); ++i )
				{
					double bin_center = binned_data[combinationIndex]->GetX()[i];

					//cout << "bin_center:\t" << bin_center << endl;
					double this_bin = PDF_function->Eval( bin_center );// cout << endl;
					//cout << "this_bin:\t" << this_bin << "\t";
					allPullData[i]+= this_bin;
					//cout << allPullData[i] << "\t" << endl;
				}

				TString desc_pull( desc );
				desc_pull.Append("_pull");

				cout << endl << "Making Plots" << endl;

				TGraphErrors* pullPlot = ComponentPlotter::PullPlot1D( allPullData, binned_data.back(), observableName, desc_pull.Data(), plotPDF->GetRandomFunction(), this_config );

				/*
				   for( unsigned int i=0; i< (unsigned)binned_data[combinationIndex]->GetN(); ++i )
				   {
				   cout << pullPlot->GetY()[i] <<  endl;
				   }
				   */

				this->OutputPlot( binned_data.back(), these_components, observableName, string(desc_pull.Data()), plotData->GetBoundary(), plotPDF->GetRandomFunction(), temp, debug, allPullData );
				(void)pullPlot;
			}
		}
		delete temp;
	}

	//	If there is more than 1 combination it's useful to plot the total's on the same graph with total component 0
	if( this->GetComponents().size() > 1 )
	{
		vector<TGraph*> allZerothComponents;
		vector<vector<TGraph*> > allGraphs = this->GetComponents();
		for( unsigned int i=0; i< allGraphs.size(); ++i )
		{
			allZerothComponents.push_back( allGraphs[i][0] );
		}

		string desc("_All_Combinations");
		TString PDFDesc;PDFDesc+=PDFNum;
		desc.append("_PDF_"); desc.append(PDFDesc.Data());

		vector<string> combDescs;
		if( this_config->combination_names.empty() )
		{
			combDescs.push_back( "Total");
			for( unsigned int i=0; i< combinationDescriptions.size(); ++i )
			{
				combDescs.push_back( combinationDescriptions[i] );
			}
		}

		CompPlotter_config* temp = new CompPlotter_config( *this_config );

		vector<int> colorVec;

		for( unsigned int i=0; i< combinationDescriptions.size(); ++i )
		{
			colorVec.push_back( (int)(i+2) );
		}

		if( temp->color_key.empty() )
		{
			temp->color_key = colorVec;
		}

		if( temp->combination_names.empty() )	temp->component_names = combDescs;
		else temp->component_names = temp->combination_names;

		ComponentPlotter::OutputPlot( binned_data[0], allZerothComponents, observableName, desc, plotData->GetBoundary(), plotPDF->GetRandomFunction(), temp, debug );

		if( !allPullData.empty() )
		{
			desc.append("_wPulls");
			ComponentPlotter::OutputPlot( binned_data[0], allZerothComponents, observableName, desc, plotData->GetBoundary(), plotPDF->GetRandomFunction(), temp, debug, allPullData );
		}

		delete temp;
	}

}

//	Again I won't be bothered coding this to check the input, this is the USERS problem!
vector<TGraph*> ComponentPlotter::MergeComponents( vector<vector<TGraph*> > input, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;

	if( input.size() == 1 ) return input[0];

	vector<TGraph*> output_components;
	for( unsigned int component_i=0; component_i< input[0].size(); ++component_i )
	{
		//	Collect all component 0 into a vector or component 1 into a seperate vector etc...
		vector<TGraph*> this_component;
		for( unsigned int combination_i = 0; combination_i < input.size(); ++combination_i )
		{
			this_component.push_back( input[combination_i][component_i] );
		}

		//	Copy all of the data into a vector of vector<double>
		vector<vector<double> > X_val, Y_val;
		for( vector<TGraph*>::iterator graph_i = this_component.begin(); graph_i != this_component.end(); ++graph_i )
		{
			int data_num = (*graph_i)->GetN();

			double* x_pointer = (*graph_i)->GetX();
			double* y_pointer = (*graph_i)->GetY();

			vector<double> this_x, this_y;

			for( int i=0; i < data_num; ++i )
			{
				this_x.push_back( x_pointer[i] );
				this_y.push_back( y_pointer[i] );
			}

			X_val.push_back( this_x );
			Y_val.push_back( this_y );
		}


		vector<double> final_X_val(X_val[0].size()), final_Y_val(Y_val[0].size());

		//      Sum the X and Y values
		for( unsigned int i=0; i< X_val[0].size(); ++i )
		{
			final_X_val[i] = X_val[0][i];
			final_Y_val[i] = 0.;
			for( unsigned int j=0; j< X_val.size(); ++j )
			{
				final_Y_val[i] += Y_val[j][i];
			}
		}

		//      Objects must have unqiue names, even though they exist in different memory locations, THANK YOU ROOT GARBAGE COLLECTION!!! (this is the singally worst idea ever to grace c++!)
		TString TGraphName("TGraph_");TGraphName+=rand->Rndm();

		string TGRaphCleanName( TGraphName.Data() );
		replace( TGRaphCleanName.begin(), TGRaphCleanName.end(), '.', '_' );

		TGraph* output_graph = new TGraph( (Int_t)X_val[0].size(), &(final_X_val[0]), &(final_Y_val[0]) );
		output_graph->SetName( TGRaphCleanName.c_str() );
		output_graph->SetTitle("");
		output_graph->SetLineColor( this_component[0]->GetLineColor() );
		output_graph->SetLineWidth( this_component[0]->GetLineWidth() );
		output_graph->SetLineStyle( this_component[0]->GetLineStyle() );
		output_graph->SetMarkerColor( this_component[0]->GetMarkerColor() );
		output_graph->SetMarkerStyle( this_component[0]->GetMarkerStyle() );
		output_components.push_back( output_graph );
	}

	return output_components;
}

//	I won't be bothered coding this in order to check that I have been provided with compatible input. The USER should check this!
TGraphErrors* ComponentPlotter::MergeBinnedData( vector<TGraphErrors*> input, TRandom* rand )
{
	if( rand == NULL ) rand = gRandom;

	if( input.size() == 1 ) return input[0];

	vector<vector<double> > X_val, X_err, Y_val, Y_err;
	for( vector<TGraphErrors*>::iterator graph_i = input.begin(); graph_i != input.end(); ++graph_i )
	{
		int data_num = (*graph_i)->GetN();

		double* x_pointer = (*graph_i)->GetX();
		double* y_pointer = (*graph_i)->GetY();
		double* x_err_pointer = (*graph_i)->GetEX();
		double* y_err_pointer = (*graph_i)->GetEY();

		vector<double> this_x, this_y, this_ex, this_ey;

		for( int i=0; i < data_num; ++i )
		{
			this_x.push_back( x_pointer[i] );
			this_y.push_back( y_pointer[i] );
			this_ex.push_back( x_err_pointer[i] );
			this_ey.push_back( y_err_pointer[i] );
		}

		X_val.push_back( this_x );
		Y_val.push_back( this_y );
		X_err.push_back( this_ex );
		Y_err.push_back( this_ey );
	}

	vector<double> final_X_val(X_val[0].size()), final_X_err(X_val[0].size()), final_Y_val(X_val[0].size()), final_Y_err(X_val[0].size());

	//	Sum the X and Y values and add the errors in quadrature
	for( unsigned int i=0; i< X_err[0].size(); ++i )
	{
		double y_err_sq = 0.;
		final_X_val[i] = X_val[0][i];
		final_X_err[i] = X_err[0][i];
		final_Y_val[i] = 0.;
		for( unsigned int j=0; j< Y_err.size(); ++j )
		{
			final_Y_val[i] += Y_val[j][i];
			y_err_sq += Y_err[j][i]*Y_err[j][i];
		}
		final_Y_err[i] = sqrt( y_err_sq );
	}

	//      Objects must have unqiue names, even though they exist in different memory locations, THANK YOU ROOT GARBAGE COLLECTION!!! (this is the singally worst idea ever to grace c++!)
	TString TGraphErrorsName("TGraphErrors_");TGraphErrorsName+=rand->Rndm();

	string TGraphErrorsCleanName( TGraphErrorsName.Data() );
	replace( TGraphErrorsCleanName.begin(), TGraphErrorsCleanName.end(), '.', '_' );

	TGraphErrors* output_graph = new TGraphErrors( (Int_t)X_val[0].size(), &(final_X_val[0]), &(final_Y_val[0]), &(final_X_err[0]), &(final_Y_err[0]) );
	output_graph->SetName( TGraphErrorsCleanName.c_str() );
	output_graph->SetTitle("");

	output_graph->SetLineColor( input[0]->GetLineColor() );
	output_graph->SetLineWidth( input[0]->GetLineWidth() );
	output_graph->SetLineStyle( input[0]->GetLineStyle() );
	output_graph->SetMarkerColor( input[0]->GetMarkerColor() );
	output_graph->SetMarkerStyle( input[0]->GetMarkerStyle() );

	return output_graph;
}


//	Plot all components on this combinations and print and save the canvas
void ComponentPlotter::OutputPlot( TGraphErrors* input_data, vector<TGraph*> input_components, string observableName,
		string CombinationDescription, PhaseSpaceBoundary* total_boundary, TRandom* rand, CompPlotter_config* conf, DebugClass* debug, vector<double> input_bin_theory_data )
{
	if( rand == NULL ) rand = gRandom;
	TString TCanvas_Name("Overlay_"+observableName+"_"+CombinationDescription+"_");TCanvas_Name+=rand->Rndm();

	string TCanvasCleanName( TCanvas_Name.Data() );
	replace( TCanvasCleanName.begin(), TCanvasCleanName.end(), '.', '_' );

	bool logy=false;
	bool logx=false;

	TCanvas* c1 = EdStyle::RapidFitCanvas( TCanvasCleanName.c_str(), "" );

	if( conf != NULL )
	{
		if( conf->logY )
		{
			logy=true;
			c1->SetLogy( true );
			//c1->Update();
		}
		if( conf->logX )
		{
			logx=true;
			c1->SetLogx( true );
		}
	}

	TPad* pad1=NULL;	//		Holds Projection Plot
	TPad* pad2=NULL;	//		Holds Pull Plot

	if( !input_bin_theory_data.empty() )
	{
		pad1 = new TPad("pad1","pad1", 0., 0.3, 1., 1.);
		pad2 = new TPad("pad2","pad2", 0., 0.0, 1., 0.3);
		pad1->Draw();
		pad2->Draw();

		pad1->cd();
	}

	TString plotTitle;

	vector<int> Style_Key, Color_Key, Width_Key;

	vector<string> component_names;

	double X_min=-99999, X_max=-99999;
	double Y_min=-99999, Y_max=-99999;

	TString X_Title, Y_Title;

	double final_chi2=-99999;
	double legend_size=0.1;

	bool addLHCb=false;
	bool addRightLHCb=false;
	bool limitPulls=false;
	bool drawSpline=true;

	double XaxisTitleScale=1.;
	double XaxisLabelScale=1.;
	double YaxisTitleScale=1.;
	double YaxisLabelScale=1.;
	if( conf != NULL )
	{
		if( conf->logY )
		{
			logy=true;
			if( !input_bin_theory_data.empty() )
			{
				pad1->SetLogy( true );
			}
			else
			{
				c1->SetLogy( true );
			}
			//c1->Update();
		}
		if( conf->logX )
		{
			logx=true;
			if( !input_bin_theory_data.empty() )
			{
				pad1->SetLogx( true );
				pad2->SetLogx( true );
			}
			else
			{
				c1->SetLogx( true );
			}
		}
		plotTitle = conf->PlotTitle;
		Style_Key = conf->style_key;
		Color_Key = conf->color_key;
		Width_Key = conf->width_key;
		component_names = conf->component_names;
		legend_size = conf->LegendTextSize;
		addLHCb = conf->addLHCb;
		addRightLHCb = conf->addRightLHCb;
		X_min = conf->xmin;
		X_max = conf->xmax;
		Y_min = conf->ymin;
		Y_max = conf->ymax;
		X_Title = conf->xtitle;
		Y_Title = conf->ytitle;
		final_chi2 = conf->Chi2Value;
		limitPulls = conf->LimitPulls;
		drawSpline = conf->useSpline;
		XaxisTitleScale = conf->XaxisTitleScale;
		XaxisLabelScale = conf->XaxisLabelScale;
		YaxisTitleScale = conf->YaxisTitleScale;
		YaxisLabelScale = conf->YaxisLabelScale;
	}

	input_data->SetTitle( plotTitle );
	input_data->Draw("AP9");
	if( !input_bin_theory_data.empty() )
	{
		pad1->Modified();
		pad1->Update();
	}
	if( !input_bin_theory_data.empty() )
	{
		if( logy ) pad1->SetLogy( true );
		if( logx )
		{
			pad1->SetLogx( true );
			pad2->SetLogx( true );
		}
	}
	else
	{
		if( logy ) c1->SetLogy( true );
		if( logx ) c1->SetLogx( true );
	}
	c1->Update();

	if( X_min <= -99999 ) X_min = total_boundary->GetConstraint( observableName )->GetMinimum();
	if( X_max <= -99999 ) X_max = total_boundary->GetConstraint( observableName )->GetMaximum();
	if( Y_min <= -99999 ) Y_min = input_data->GetYaxis()->GetXmin();//logy==true?0.5:0.;
	if( Y_max <= -99999 ) Y_max = input_data->GetYaxis()->GetXmax();

	TString Y_ext(" / ( ");
	stringstream thisStream;
	if( conf != NULL )
	{
		thisStream << setw(4) << setprecision(3) /*<< scientific*/ << fabs(X_max-X_min)/((double)(conf->data_bins));
	}
	else
	{
		thisStream << setw(4) << setprecision(3) /*<< scientific*/ << fabs(X_max-X_min)/100.;
	}
	Y_ext.Append( thisStream.str() );
	Y_ext.Append(" ");
	TString Unit; Unit.Append( EdStyle::GetParamRootUnit( observableName ) );
	if( !StringProcessing::is_empty(Unit) ) Unit.Append(" ");
	Y_ext.Append(Unit); Y_ext.Append(")");

	if( StringProcessing::is_empty( X_Title ) )
	{
		X_Title = EdStyle::GetParamRootName( observableName );
		TString unit = EdStyle::GetParamRootUnit( observableName );
		if( !StringProcessing::is_empty( unit ) )
		{
			X_Title.Append( " " );
			X_Title.Append( unit );
		}
	}
	if( StringProcessing::is_empty( Y_Title ) )
	{
		Y_Title = "Candidates";
		Y_Title.Append( Y_ext );
	}

	input_data->GetYaxis()->SetRangeUser( Y_min, Y_max );
	input_data->GetYaxis()->SetTitle( Y_Title );
	input_data->GetYaxis()->SetTitleSize( (Float_t)YaxisTitleScale*input_data->GetYaxis()->GetTitleSize() );
	input_data->GetYaxis()->SetLabelSize( (Float_t)YaxisLabelScale*input_data->GetYaxis()->GetLabelSize() );
	input_data->GetXaxis()->SetRangeUser( X_min, X_max );
	input_data->GetXaxis()->SetTitle( X_Title );
	input_data->GetXaxis()->SetTitleSize( (Float_t)XaxisTitleScale*input_data->GetXaxis()->GetTitleSize() );
	input_data->GetXaxis()->SetLabelSize( (Float_t)XaxisLabelScale*input_data->GetXaxis()->GetLabelSize() );

	if( !input_bin_theory_data.empty() )
	{
		pad1->Modified();
		pad1->Update();
	}

	c1->Modified();
	c1->Update();

	//input_data->GetXaxis()->CenterTitle( true );

	if( !input_bin_theory_data.empty() )
	{
		pad1->Modified();
		pad1->Update();

		if( logy ) pad1->SetLogy( true );
		if( logx )
		{
			pad1->SetLogx( true );
			pad2->SetLogx( true );
		}
	}
	else
	{
		if( logy ) c1->SetLogy( true );
		if( logx ) c1->SetLogx( true );
	}

	c1->Update();


	TPaveText* myLatex=NULL;
	if( addLHCb )
	{
		myLatex = EdStyle::LHCbLabel();
	}
	if( addRightLHCb )
	{
		myLatex = EdStyle::RightLHCbLabel();
	}

	if( myLatex != NULL ) myLatex->Draw();

	TLegend* leg = NULL;
	if( conf->useLegend )
	{
		if( conf->TopRightLegend )
		{
			leg = EdStyle::LHCbLegend();
		} else if( conf->TopLeftLegend )
		{
			leg = EdStyle::LHCbLeftLegend();
		} else if( conf->BottomRightLegend )
		{
			leg = EdStyle::LHCbBottomLegend();
		} else if( conf->BottomLeftLegend )
		{
			leg = EdStyle::LHCbBottomLeftLegend();
		}
	}

	if( leg != NULL ) leg->SetTextSize( (Float_t) legend_size );

	if( leg != NULL ) leg->AddEntry( input_data, "Data", "pl" );

	unsigned int num=0;
	for( vector<TGraph*>::iterator comp_i = input_components.begin(); comp_i != input_components.end(); ++comp_i, ++num )
	{
		if( !Style_Key.empty() )
		{
			if( num < Style_Key.size() )
			{
				(*comp_i)->SetLineStyle( (Style_t)Style_Key[num] );
			}
		}
		if( !Color_Key.empty() )
		{
			if( num < Color_Key.size() )
			{
				(*comp_i)->SetLineColor( (Color_t)Color_Key[num] );
			}
		}
		if( !Width_Key.empty() )
		{
			if( num < Width_Key.size() )
			{
				(*comp_i)->SetLineWidth( (Width_t)Width_Key[num] );
			}
		}

		if( (*comp_i)->GetLineWidth() != 0 )
		{
			if( drawSpline )
			{
				(*comp_i)->Draw("C9");
			}
			else
			{
				(*comp_i)->Draw("L9");
			}

			if( !component_names.empty() )
			{
				if( num < component_names.size() )
				{
					if( leg != NULL ) leg->AddEntry( (*comp_i), TString(component_names[num]), "l" );
				}
				else
				{
					if( leg != NULL ) leg->AddEntry( (*comp_i), "Unnamed", "l" );
				}
			}
		}
	}

	if( final_chi2 > 0 )
	{
		TString Chi2Text( "#chi^{2}/ndof : " );
		stringstream chi_stream; chi_stream << setprecision(4) << final_chi2;
		Chi2Text.Append( chi_stream.str() );
		if( leg !=NULL ) leg->AddEntry( (TObject*)NULL, Chi2Text, "" );
	}

	if( !input_bin_theory_data.empty() )
	{
		pad1->Modified();
		pad1->Update();
	}

	c1->Update();

	if( !component_names.empty() ) if( leg != NULL ) leg->Draw();

	if( !input_bin_theory_data.empty() )
	{
		pad1->Modified();
		pad1->Update();
	}

	if( !input_bin_theory_data.empty() )
	{
		if( logy ) pad1->SetLogy( true );
		if( logx )
		{
			pad1->SetLogx( true );
			pad2->SetLogx( true );
		}
	}
	else
	{
		if( logy ) c1->SetLogy( true );
		if( logx ) c1->SetLogx( true );
	}

	c1->Update();


	if( !input_bin_theory_data.empty() )
	{
		pad2->cd();

		vector<double> pull_value, pull_error_value,  x_values,  x_errs;

		for( unsigned int i=0; i< input_bin_theory_data.size(); ++i )		//	Explicit Assumption that theory.size() == input_data->GetN()
		{
			double theory_y = input_bin_theory_data[i];

			double data_y = input_data->GetY()[i];
			double data_err = input_data->GetErrorY( (int)i );

			double pull = ( data_y -theory_y ) / data_err;

			/*!
			 * From what I understand in comparing the PDF to data in this way the residuals are normalised with an error of 1.
			 *
			 * I don't implicitly trust this, however I don't have a better solution to calculate a sensible error at the time of writing
			 */
			double pull_err = 1.;

			if( pull >= DBL_MAX || data_err < 1E-10 )	pull = 0.;
			if( fabs(data_y) < 1E-10 ) pull_err = 0.;

			pull_value.push_back( pull );
			pull_error_value.push_back( pull_err );

			x_values.push_back( input_data->GetX()[i] );
			x_errs.push_back( input_data->GetErrorX( (int)i ) );

			//cout << pull << "\t" << data_y << "-" << theory_y << "/" << data_err << "\t\t\t" << pull_err << "\t\t\t" << x_values.back() << "\t" << x_errs.back() << endl;
		}
		//cout << endl;

		TGraphErrors* pullGraph = new TGraphErrors( (int)pull_value.size(), &(x_values[0]), &(pull_value[0]), &(x_errs[0]), &(pull_error_value[0]) );
		pullGraph->Draw("AP9");
		pad2->Modified();
		pad2->Update();
		c1->Update();

		pullGraph->GetYaxis()->SetTitle( "Pull" );
		pullGraph->GetYaxis()->SetTitleSize( input_data->GetYaxis()->GetTitleSize() );
		pullGraph->GetYaxis()->SetLabelSize( input_data->GetYaxis()->GetLabelSize() *(Float_t)(1./0.5) );
		pullGraph->GetXaxis()->SetTitleSize( input_data->GetXaxis()->GetTitleSize() );
		pullGraph->GetXaxis()->SetLabelSize( input_data->GetXaxis()->GetLabelSize() *(Float_t)(1./0.5) );
		pullGraph->GetXaxis()->SetRangeUser( X_min, X_max );
		//pullGraph->GetXaxis()->CenterTitle( true );
		//pullGraph->GetXaxis()->SetTitle( X_Title );
		pad2->Modified();
		pad2->Update();
		pullGraph->GetYaxis()->SetTitleOffset((Float_t)(pullGraph->GetYaxis()->GetTitleOffset()/3.));
		pullGraph->GetYaxis()->SetTitleSize((Float_t)(pullGraph->GetYaxis()->GetTitleSize()*2.));
		pullGraph->GetXaxis()->SetTitleOffset((Float_t)(pullGraph->GetXaxis()->GetTitleOffset()/2.));
		pullGraph->GetXaxis()->SetTitleSize((Float_t)(pullGraph->GetXaxis()->GetTitleSize()*2.));
		if( limitPulls )
		{
			pullGraph->GetYaxis()->SetNdivisions( 3 );
			pullGraph->GetYaxis()->SetRangeUser( -5., 5. );
			pullGraph->SetMaximum( 5. );
			pullGraph->SetMinimum( -5. );
		}
		pad2->Modified();
		pad2->Update();
		c1->Update();
	}

	c1->Write("",TObject::kOverwrite);

	TString Clean_Description = StringProcessing::Clean( CombinationDescription.c_str() );

	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".C") );
	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".pdf") );
	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".png") );

}

//      This is a slimmed down version of the AddBranch function I have written elsewhere
void ComponentPlotter::WriteBranch( TTree* input_tree, TString Branch_Name, vector<double>* branch_data )
{
	if( input_tree->GetEntries() != 0 )
	{
		if( branch_data->size() != (unsigned)input_tree->GetEntries() )
		{
			return;
		}
	}

	Double_t X_object=-1.;

	TBranch* this_X_branch = input_tree->Branch( Branch_Name, &X_object, TString(Branch_Name+"/D") );

	input_tree->SetEntries( (int)branch_data->size() );

	vector<double>::iterator dat_x_i = branch_data->begin();

	for( ; dat_x_i != branch_data->end(); ++dat_x_i )
	{
		X_object = *dat_x_i;
		this_X_branch->Fill();
	}
}

//	Initialize the Desciptions of the Relevant Combinations in the DataSet
void ComponentPlotter::SetupCombinationDescriptions()
{
	cout << "Optimally:" << endl << endl;

	//      If we have 'no discrete combinations' then get the whole dataset for this pdf
	vector<DataPoint*> whole_dataset;

	cout << plotData->GetDataNumber(NULL) << endl;

	for( int i=0; i< plotData->GetDataNumber(NULL); ++i )
	{
		whole_dataset.push_back( plotData->GetDataPoint( i ) );
	}

	data_subsets.push_back( whole_dataset );

	ObservableRef ObservableName( observableName );

	vector<double> discrete_observableValues_i;
	//      Get all values of the wanted parameter from this subset
	for( unsigned int j=0; j< whole_dataset.size(); ++j )
	{
		discrete_observableValues_i.push_back( whole_dataset[j]->GetObservable( ObservableName )->GetValue() );
	}

	cout << "Combination: 1 " << "Min: " << StatisticsFunctions::Minimum( discrete_observableValues_i );
	cout << "\tMax: " << StatisticsFunctions::Maximum( discrete_observableValues_i ) << "\tBinNum: " << StatisticsFunctions::OptimumBinNumber( discrete_observableValues_i ) << endl;

	combinationDescriptions.push_back( "All_Data" );

	discreteNames = full_boundary->GetDiscreteNames();

	//      If we have discrete combinations then we need to loop over each set seperatley to split up the data and information about the data
	if( allCombinations.size() > 1 )
	{
		combinationDescriptions.clear();
		//      For all discrete combinations that have been calculated elsewhere
		for( unsigned int combinationIndex=0; combinationIndex< allCombinations.size(); ++combinationIndex )
		{
			cout << "Combination: " << combinationIndex+2 << " ";
			//      For this combination work out what values correspond to the discrete parameters
			vector<double> values_for_this_combination;
			for( unsigned int paramIndex=0; paramIndex< discreteNames.size(); ++paramIndex )
			{
				values_for_this_combination.push_back( allCombinations[combinationIndex]->GetObservable( discreteNames[paramIndex] )->GetValue() );
			}

			//              Get all data which falls within this paramreter set
			data_subsets.push_back( plotData->GetDiscreteSubSet( discreteNames, values_for_this_combination ) );

			vector<double> this_discrete_observableValues_i;
			//      Get all values of the wanted parameter from this subset
			for( unsigned int j=0; j< data_subsets.back().size(); ++j )
			{
				this_discrete_observableValues_i.push_back( (data_subsets.back())[j]->GetObservable( ObservableName )->GetValue() );
			}

			cout << "Min: " << StatisticsFunctions::Minimum( this_discrete_observableValues_i );
			cout << "\tMax: " << StatisticsFunctions::Maximum( this_discrete_observableValues_i );
			cout << "\tBinNum: " << StatisticsFunctions::OptimumBinNumber( this_discrete_observableValues_i ) << endl;


			TString ThisName;
			for( unsigned int paramIndex=0; paramIndex< discreteNames.size(); ++paramIndex )
			{
				ThisName.Append( discreteNames[paramIndex] );
				ThisName.Append( "=" );
				ThisName+= values_for_this_combination[paramIndex];
				ThisName.Append( "_" );
			}
			cout << "Adding This: " << ThisName << endl;
			combinationDescriptions.push_back( ThisName.Data() );
		}
	}

	return;
}

//Return the values tracing the PDF projection in the Observable of choice
vector<double>* ComponentPlotter::ProjectObservableComponent( DataPoint* InputPoint, string ObservableName, double Minimum, int PlotNumber, double PlotInterval, string component )
{
	//	Move initializer(s) outside of for loop
	double integralvalue=-1., observableValue=-1.;
	Observable* newObs=NULL;

	//Find the value of the observable projection at each data point
	vector<double>* pointValues = new vector<double>();

	cout << ObservableName << ": " << Minimum << " <-> " << Minimum + ( PlotInterval * (PlotNumber-1) ) << endl;

	//	This class object has been created to speed up the communication between this class and the basePDF as it may pass through several PDF wrappers
	ComponentRef* comp_obj = new ComponentRef( component, ObservableName );

	string unit = InputPoint->GetObservable( string(ObservableName) )->GetUnit();

	ObservableRef thisObservableRef( ObservableName );

	//	Step over the whole observable range
	for (int pointIndex = 0; pointIndex < PlotNumber; ++pointIndex )
	{
		//	Inform the user of how far we have got :D
		cout << left << setw(5) << setprecision(3)  << (pointIndex+1.)/PlotNumber*100 << "\% Complete" << setw(20) << " " << "\r" << flush;

		//	Value of Observable we want to evaluate for this step
		observableValue = Minimum + ( PlotInterval * pointIndex );

		//	Set the value of Observable to the new step value
		//	All other CONTINUOUS observables set to average of range
		newObs = new Observable( string(ObservableName), observableValue, unit );
		DataPoint* thisPoint = new DataPoint( *InputPoint );
		thisPoint->SetObservable( thisObservableRef, newObs );
		delete newObs; newObs = NULL;

		if( debug != NULL )
		{
			if( debug->DebugThisClass( "ComponentPlotter" ) )
			{
				cout << endl << "ComponentPlotter: Calling RapidFitIntegrator::ProjectObservable " << pointIndex << " of " << PlotNumber << "\t" << ObservableName << "==" << observableValue << "\t";
				//InputPoint->Print();
				pdfIntegrator->SetDebug( debug );
			}
		}

		//	perform actual evaluation of the PDF with the configuration in the InputPoint, in the whole boundary with the given name and for selected component
		integralvalue = pdfIntegrator->ProjectObservable( thisPoint, full_boundary, ObservableName, comp_obj );

		if( debug != NULL )
		{
			if( debug->DebugThisClass( "ComponentPlotter" ) )
			{
				cout << "ComponentPlotter: \tI got: " << integralvalue << endl;
				InputPoint->Print();
				cout << "ComponentPlotter: Using: " << plotPDF->GetLabel() << "\tProjecting Component: " << component << "\tObservable: " << observableName << endl;
				full_boundary->Print();
			}
		}

		pointValues->push_back( integralvalue );
		delete thisPoint;
	}

	cout << "Finished Projecting" << endl;

	if( debug != NULL )
	{
		if( debug->DebugThisClass("ComponentPlotter") )
		{
			cout << endl << "ComponentPlotter: Performing Sanity Check" << endl;
		}
	}

	this->Sanity_Check( pointValues, component );

	//cout << "Returning pointValues" << endl;

	delete comp_obj;

	return pointValues;
}

//	Check that the passed vector of points aren't all NULL if they are print some debug information
void ComponentPlotter::Sanity_Check( vector<double>* pointValues, TString component )
{
	bool result = true;
	for( vector<double>::iterator inf = pointValues->begin(); inf != pointValues->end(); ++inf )
	{
		if( fabs( *inf ) > 1E-6 )
		{
			result = false;
		}
	}
	if( result )
	{
		cerr << "Failed to get anything other than 0 for PDF!!" << endl;
		full_boundary->Print();
		cerr << "component: " << component << endl;
		cerr << "\nCheck your Plots!" << endl;
		//exit(-99);
	}
}

void ComponentPlotter::SetWeightsWereUsed( string input )
{
	if( debug != NULL )
	{
		if( debug->DebugThisClass( "ComponentPlotter" ) )
		{
			cout << "ComponentPlotter: Setting Weights Were Used:" << endl;
			cout << "ComponentPlotter: WeightName: " << input << endl;

		}
	}

	double sum=0.;
	for( int i=0 ; i < plotData->GetDataNumber(NULL); ++i )
	{
		sum+=plotData->GetDataPoint( i )->GetObservable( input )->GetValue();
	}

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "ComponentPlotter" ) )
		{
			cout << "ComponentPlotter: Sum of Weights: " << sum << endl;
			cout << "ComponentPlotter: DataSet Size: " << plotData->GetDataNumber(NULL) << endl;
		}
	}
	weightsWereUsed = true;
	weightName = input;

	double alpha=1.;

	if( fabs(plotData->GetAlpha()-1.)<1E-6 )
	{
		alpha = plotData->GetSumWeights()/plotData->GetSumWeightsSq();
	}
	else
	{
		alpha = plotData->GetAlpha();
	}

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "ComponentPlotter" ) )
		{
			cout << "ComponentPlotter: Alpha: " << alpha << endl;
		}
	}

	ObservableRef* weight_ref = new ObservableRef( weightName );
	for( int i=0 ; i < plotData->GetDataNumber(NULL); ++i )
	{
		wanted_weights.push_back( plotData->GetDataPoint( i )->GetObservable( *weight_ref )->GetValue() );
	}
	delete weight_ref;

	weight_norm=0.;
	double weight_sum=plotData->GetSumWeights();

	weight_norm = weight_sum / plotData->GetDataNumber(NULL);

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "ComponentPlotter" ) )
		{
			cout << "Weight Norm: " << weight_norm << endl;
		}
	}
}

//	Functions to access the internal results in this class
vector<TGraphErrors*> ComponentPlotter::GetBinnedData()
{
	return binned_data;
}

vector<vector<TGraph*> > ComponentPlotter::GetComponents()
{
	return total_components;
}

double ComponentPlotter::operator() (double *x, double *p)
{
	(void) p;			//	NOT USED OR UNDERSTOOD IN THIS CONTEXT
	double fixed_obs_val = x[0];

	Observable* oldObs=NULL;
	Observable* newObs=NULL;

	double integral_value=0.;

	ComponentRef* comp_obj = new ComponentRef( "0", observableName );

	unsigned int comb_num=0;

	for( vector<DataPoint*>::iterator comb_i = allCombinations.begin();  comb_i != allCombinations.end(); ++comb_i, ++comb_num )
	{
		DataPoint* InputPoint = new DataPoint( *( (*comb_i) ) );

		oldObs = InputPoint->GetObservable( observableName );
		newObs = new Observable( observableName, fixed_obs_val, oldObs->GetUnit() );
		InputPoint->SetObservable( observableName, newObs );
		//delete oldObs;

		double PDFProjection = pdfIntegrator->ProjectObservable( InputPoint, full_boundary, observableName, comp_obj );

		double PDFNormalisation = this->PDF2DataNormalisation( comb_num );

		integral_value += PDFProjection*PDFNormalisation;

		delete InputPoint;
	}

	delete comp_obj;

	cout << "Value At: " << left << setw(5) << setprecision(3) << x[0] << "\tis:\t" << setprecision(4) << integral_value << setw(20) << " " <<  "\r" << flush;

	return integral_value;
}

double ComponentPlotter::PDF2DataNormalisation( const unsigned int combinationIndex ) const
{
	double normalisation=1.;

	normalisation *= fabs(ratioOfIntegrals[ combinationIndex ]);			//	Attempt to correct for Numerical != analytical due to any constant factor due to numerical inaccuracy
	//	(some constant close to 1. exactly 1. for numerical PDFs)

	double dataNum = plotData->GetDataNumber( allCombinations[combinationIndex] );
	normalisation *= dataNum / (double) data_binning;				//	Normalise to this density of events	(Num of events per bin in flatPDF)

	double range = fabs( boundary_max-boundary_min );
	normalisation *= range;								//	Correct for the range of the dataset	(absolute range of Observable being projected)

	normalisation /= combination_integral[ combinationIndex ];			//	Total Integral of the PDF	(We're plotting prob of PDF being at this point for a non-unitary PDF)

	normalisation *= weight_norm;							//	Correct for the effect of non-unitary weights used in the fit

	return normalisation;
}

pair<double,double> ComponentPlotter::GetChi2Numbers()
{
	return make_pair( chi2, N );
}

void ComponentPlotter::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
}


vector<double> ComponentPlotter::GetFunctionEval()
{
	return allPullData;
}

TGraphErrors* ComponentPlotter::PullPlot1D( vector<double> input_bin_theory_data, TGraphErrors* input_data, string observableName, string CombinationDescription, TRandom* rand, CompPlotter_config* conf, DebugClass* debug )
{
	if( rand == NULL ) rand = gRandom;

	vector<double> pull_value;
	vector<double> pull_error_value;

	vector<double> x_values;
	vector<double> x_errs;

	//cout << endl;
	for( unsigned int i=0; i< input_bin_theory_data.size(); ++i )		//	Explicit Assumption that theory.size() == input_data->GetN()
	{
		double theory_y = input_bin_theory_data[i];

		double data_y = input_data->GetY()[i];
		double data_err = input_data->GetErrorY( (int)i );

		double pull = ( data_y - theory_y ) / data_err;

		/*!
		 * From what I understand in comparing the PDF to data in this way the residuals are normalised with an error of 1.
		 *
		 * I don't implicitly trust this, however I don't have a better solution to calculate a sensible error at the time of writing
		 */
		double pull_err = 1.;

		if( pull >= DBL_MAX || data_err < 1E-10 )	pull = 0.;
		if( fabs(data_y) < 1E-10 ) pull_err = 0.;

		pull_value.push_back( pull );
		pull_error_value.push_back( pull_err );

		x_values.push_back( input_data->GetX()[i] );
		x_errs.push_back( input_data->GetErrorX( (int)i ) );

		//cout << pull << "\t" << data_y << "-" << theory_y << "/" << data_err << "\t\t\t" << pull_err << "\t\t\t" << x_values.back() << "\t" << x_errs.back() << endl;
	}
	//cout << endl;

	TGraphErrors* pullGraph = new TGraphErrors( (int)pull_value.size(), &(x_values[0]), &(pull_value[0]), &(x_errs[0]), &(pull_error_value[0]) );

	TString pullGraphName="PullGraph_"; pullGraphName+=rand->Rndm();

	string pullGraphCleanName( pullGraphName.Data() );
	replace( pullGraphCleanName.begin(), pullGraphCleanName.end(), '.', '_' );

	pullGraph->SetName( pullGraphCleanName.c_str() ); pullGraph->SetTitle( pullGraphCleanName.c_str() );

	TString canvas_name = "PullPlot_"; canvas_name.Append(observableName); canvas_name.Append("_"); canvas_name+=rand->Rndm();

	string canvas_clean_name( canvas_name.Data() );
	replace( canvas_clean_name.begin(), canvas_clean_name.end(), '.', '_' );

	TCanvas* c1 = EdStyle::RapidFitCanvas( canvas_clean_name.c_str(), canvas_name );
	pullGraph->Draw("AP9");
	c1->Update();

	if( conf != NULL )
	{
		if( conf->LimitPulls )
		{
			pullGraph->GetYaxis()->SetRangeUser(-5.,5.);
			c1->Update();
		}
	}

	//pullGraph->GetYaxis()->SetRangeUser(-5.,5.);

	TString X_Title;
	TString Y_Title;

	if( StringProcessing::is_empty( X_Title ) )
	{
		X_Title = EdStyle::GetParamRootName( observableName );
		TString unit = EdStyle::GetParamRootUnit( observableName );
		if( !StringProcessing::is_empty( unit ) )
		{
			X_Title.Append( " " );
			X_Title.Append( unit );
		}
	}
	if( StringProcessing::is_empty( Y_Title ) ) Y_Title = "Pull";

	input_data->GetYaxis()->SetTitle( Y_Title );
	//input_data->GetXaxis()->SetRangeUser( X_min, X_max );
	//input_data->GetXaxis()->CenterTitle( true );
	input_data->GetXaxis()->SetTitle( X_Title );

	c1->Update();

	TString Clean_Description = StringProcessing::Clean( CombinationDescription.c_str() );

	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".C") );
	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".pdf") );
	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".png") );

	TCanvas* c2 = EdStyle::RapidFitCanvas( "pull", "pull");
	TString pull_name = "pull_"; pull_name.Append(observableName);
	TH1D * pull_histograms = new TH1D(pull_name, pull_name, StatisticsFunctions::OptimumBinNumber( pull_value ), -5, 5);
	for ( int i = 0; i < (int)input_bin_theory_data.size() ; i++ ) pull_histograms->Fill(pull_value[(unsigned)i]);
	pull_histograms->Draw();
	gStyle->SetOptFit(1);
	pull_histograms->Fit("gaus");
	c2->Update();
	c2->Print( TString("pull"+observableName+"_"+Clean_Description+".pdf") );
    delete pull_histograms;

	return pullGraph;
}


void ComponentPlotter::WriteData( TGraphErrors* Total_BinnedData, vector<TGraph*> Total_Components, TString destination )
{
	TTree* outputTree = new TTree( destination, destination );

	unsigned int points = (unsigned)Total_BinnedData->GetN();

	vector<double> X_val;
	vector<double> X_err;
	vector<double> Y_val;
	vector<double> Y_err;
	for( unsigned int i=0; i< points; ++i )
	{
		X_val.push_back( Total_BinnedData->GetX()[i] );
		X_err.push_back( Total_BinnedData->GetEX()[i] );
		Y_val.push_back( Total_BinnedData->GetY()[i] );
		Y_err.push_back( Total_BinnedData->GetEY()[i] );
	}

	ComponentPlotter::WriteBranch( outputTree, "BinnedData_Value_X", &X_val  );
	ComponentPlotter::WriteBranch( outputTree, "BinnedData_Error_X", &X_err  );
	ComponentPlotter::WriteBranch( outputTree, "BinnedData_Value_Y", &Y_val  );
	ComponentPlotter::WriteBranch( outputTree, "BinnedData_Error_Y", &Y_err  );

	for( unsigned int i=0; i< Total_Components.size(); ++i )
	{
		vector<double> this_X_val;
		vector<double> this_Y_val;
		unsigned int this_points = (unsigned)Total_Components[i]->GetN();
		for( unsigned int j=0; j< this_points; ++j )
		{
			this_X_val.push_back( Total_Components[i]->GetX()[j] );
			this_Y_val.push_back( Total_Components[i]->GetY()[j] );
			TString ComponentName="Component_";
			ComponentName+=i;
			TString ComponentXName=ComponentName;ComponentXName.Append("_X");
			TString ComponentYName=ComponentName;ComponentYName.Append("_Y");
			ComponentPlotter::WriteBranch( outputTree, ComponentXName, &this_X_val );
			ComponentPlotter::WriteBranch( outputTree, ComponentYName, &this_Y_val );
		}
	}

	outputTree->Write("",TObject::kOverwrite);
}

string ComponentPlotter::XML( int projectionType )
{
	stringstream xml;
	string projectionOuterTag="ComponentProjection";
	if( projectionType == 2 ) projectionOuterTag = "Projection";

	xml << "<" << projectionOuterTag << ">" << endl;

	xml << "\t" << "<Name>" << "someObservable" << "</Name>" << endl;
	xml << "\t" << "<CompNames>" << "component1Name:component2Name:component3Name:..." << "</CompNames> # or <ComponentNames> refers to the PDF Components eg CP-Even:CP-Odd:Swave" << endl;
	xml << "\t" << "<CombinationNames>" << "combination1Name:combination2Name:combination3Name:..." << "</CombinationNames> # refers to the different Discrete Combinations eg tag=+1:tag=-1:untagged" << endl;
	xml << "\t" << "<Xmax>" << "XaxisMax" << "</Xmax>" << endl;
	xml << "\t" << "<Xmin>" << "XaxisMin" << "</Xmin>" << endl;
	xml << "\t" << "<Ymax>" << "YaxisMax" << "</Ymax>" << endl;
	xml << "\t" << "<Ymin>" << "YaxisMin" << "</Ymin>" << endl;
	xml << "\t" << "<XTitle>" << "XaxisTitle" << "</XTitle>" << endl;
	xml << "\t" << "<YTitle>" << "YaxisTitle" << "</YTitle>" << endl;
	xml << "\t" << "<TrustNumerical>" << "True/False" << "</TrustNumerical> # Correct for Numerical vs analytical difference" << endl;
	xml << "\t" << "<CalcChi2>" << "True/False" << "</CalcChi2> # Calculate Chi2 for the Projection" << endl;
	xml << "\t" << "<DrawPull>" << "True/False" << "</DrawPull> # Draw Residual Plot" << endl;
	xml << "\t" << "<LimitPulls>" << "True/False" << "</LimitPulls> # limit pull plot to +/- 5" << endl;
	xml << "\t" << "<AddLHCb>" << "True/False" << "</AddLHCb> # Add LHCb logo on top left" << endl;
	xml << "\t" << "<AddRightLHCb>" << "True/False" << "</AddRightLHCb> # Add LHCb logo on top right" << endl;
	xml << "\t" << "<LegendTextSize>" << "Size_t" << "</LegendTextSize> # Legend font size default: 0.05" << endl;
	xml << "\t" << "<TopRightLegend>" << "True/False" << "</TopRightLegend> # Add Legend on top right of plot" << endl;
	xml << "\t" << "<TopLeftLegend>" << "True/False" << "</TopLeftLegend> # Add Legend on top left of plot" << endl;
	xml << "\t" << "<BottomRightLegend>" << "True/False" << "</BottomRightLegend> # Add Legend on bottom right of plot" << endl;
	xml << "\t" << "<BottomLeftLegend>" << "True/False" << "</BottomLeftLegend> # Add Legend on bottom left of plot" << endl;
	xml << "\t" << "<NoLegend>" << "True/False" << "</NoLegend> # Do Not draw the Legend" << endl;
	xml << "\t" << "<UseSpline>" << "True/False" << "</UseSpline> # Use a Spline to interpolate between PDF points" << endl;
	xml << "\t" << "<Threads>" << "NumberOfThreads" << "</Threads> # Number of Threads that GSL will use when multi-threading" << endl;
	xml << "\t" << "<FixedIntegrationPoints>" << "numberOfPointsPerGSLIntegral" << "</FixedIntegrationPoints> # Set the Number of points for GSL to use to be non default" << endl;
	xml << "\t" << "<UseGSLNumericalIntegration>" << "True/False" << "</UseGSLNumericalIntegration> # Use the GSL Integrator for projections, this is multi-threaded so better" << endl;
	xml << "\t" << "<StyleKey>" << "LineStyle1:LineStyle2:LineStyle3:..." << "</StyleKey> # Styles to use for different lines" << endl;
	xml << "\t" << "<ColorKey>" << "LineColor1:LineColor2:LineColor3:..." << "</ColorKey> # Colors to use for different lines" << endl;
	xml << "\t" << "<WidthKey>" << "LineWidth1:LineWidth2:LineWidth3:..." << "</WidthKey> # Widths of lines on plot, width 0 means a line is not drawn or added to Legend" << endl;
	xml << "\t" << "<LogY>" << "True/False" << "</LogY> # Set Log on Y axis" << endl;
	xml << "\t" << "<LogX>" << "True/False" << "</LogX> # Set Log on X axis" << endl;
	xml << "\t" << "<PDFpoints>" << "numberOfPDFIntegralPoints" << "</PDFpoints> # Number of bins that the PDF will be evaluated(integrated) at in Xaxis" << endl;
	xml << "\t" << "<DataBins>" << "numberOfBinsInDataHisto" << "</DataBins> # Number of bins that will be used in constructing the histogram of data" << endl;
	xml << "\t" << "<PlotAllCombinatons>" << "True/False" << "</PlotAllCombinatons> # Do we need to project tag+/-1 separately?" << endl;
	xml << "\t" << "<DefaultDiscreteValue>" << "someValue" << "</DefaultDiscreteValue> # Advanced usage to replace a series of discrete values with a non-zero default" << endl;
	xml << "\t" << "<XaxisTitleScale>" << "someScale" << "</XaxisTitleScale> # Scale the Font Size of the X axis Title" << endl;
	xml << "\t" << "<XaxisLabelScale>" << "someScale" << "</XaxisLabelScale> # Scale the Font Size of the X axis Labels" << endl;
	xml << "\t" << "<YaxisTitleScale>" << "someScale" << "</YaxisTitleScale> # Scale the Font Size of the Y axis Title" << endl;
	xml << "\t" << "<YaxisLabelScale>" << "someScale" << "</YaxisLabelScale> # Scale the Font Size of the Y axis Labels" << endl;

	xml << "</" << projectionOuterTag << ">" << endl;

	return xml.str();
}

