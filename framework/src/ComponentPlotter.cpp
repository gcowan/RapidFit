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
///	RapidFit Headers
#include "ComponentPlotter.h"
#include "StatisticsFunctions.h"
#include "EdStyle.h"
#include "StringProcessing.h"
#include "PhaseSpaceBoundary.h"
#include "ClassLookUp.h"
#include "ComponentRef.h"
///	System Headers
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <cstdlib>
#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

//Constructor with correct arguments
ComponentPlotter::ComponentPlotter( IPDF * NewPDF, IDataSet * NewDataSet, TString PDFStr, TFile* filename, string ObservableName, CompPlotter_config* config ):
	observableName( ObservableName ), weightName(), plotPDF( ClassLookUp::CopyPDF(NewPDF) ), plotData( NewDataSet ), pdfIntegrator( new RapidFitIntegrator(NewPDF) ), weightsWereUsed(false), weight_norm(1.),
	discreteNames(), continuousNames(), full_boundary( new PhaseSpaceBoundary(*(NewDataSet->GetBoundary())) ), PlotFile( filename ),
	total_points( (config!=NULL)?config->PDF_points:128 ), data_binning( (config!=NULL)?config->data_bins:100 ), pdfStr( PDFStr ), logY( (config!=NULL)?config->logY:false ),
	this_config( config ), boundary_min( -999 ), boundary_max( -999 ), step_size( -999 ), onlyZero( (config!=NULL)?config->OnlyZero:false )
{
	plotPDF->TurnCachingOff();

	pdfIntegrator->SetPDF( plotPDF );
	pdfIntegrator->ProjectionSettings();
	//////////////////////////////////////////////////////////////////
	//Make the histogram of this observable
	format = new EdStyle();
	format->SetStyle();
	//gStyle->SetMarkerStyle(15);
	//gStyle->SetMarkerSize(Size_t(0.8));

	boundary_min = full_boundary->GetConstraint( observableName )->GetMinimum();
	boundary_max = full_boundary->GetConstraint( observableName )->GetMaximum();

	step_size = ( boundary_max - boundary_min ) / (double)( total_points - 1 );

	//Work out what to plot
	allCombinations = this->GetDiscreteCombinations( combinationWeights, combinationDescriptions );

	for( unsigned int i=0; i< allCombinations.size(); ++i )
	{
		pdfIntegrator->ForceTestStatus( false );
		combination_integral.push_back( pdfIntegrator->Integral( allCombinations[i], plotData->GetBoundary() ) );
		if( plotPDF->GetNumericalNormalisation() == true ) ratioOfIntegrals.push_back( 1. );
		else ratioOfIntegrals.push_back( pdfIntegrator->GetRatioOfIntegrals() );
	}

	vector<double> minimum, maximum;
	vector<int> binNumber;

	//Get the data needed to make the histogram
	vector<vector<double> > observableValues = GetCombinationStatistics( observableName, minimum, maximum, binNumber );

	(void) minimum; (void) maximum; (void) binNumber;
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
}

//Create a root file containing a projection plot over one observable
void ComponentPlotter::ProjectObservable()
{
	vector<string> doNotIntegrate = plotPDF->GetDoNotIntegrateList();

	for(vector<string>::iterator comb_i = combinationDescriptions.begin(); comb_i != combinationDescriptions.end(); ++comb_i ) cout << *comb_i << endl;

	plotPDF->UnsetCache();
	//Check the observable can be plotted
	bool continuous = !( plotData->GetBoundary()->GetConstraint( observableName )->IsDiscrete() );
	bool doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrate, &(observableName) ) == -1 );
	if ( continuous && doIntegrate )
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
		vector<double>* projectionValueArray = new vector<double>();
		projectionValueArray->resize((unsigned)total_points);
		double range = ( boundary_max - boundary_min );

		//	Perform Normalization							THIS IS COMPLEX AND IS COPIED LARGELY FROM PLOTTER CLASS
		for (int pointIndex = 0; pointIndex < total_points; ++pointIndex )
		{
			//Projection graph
			//
			(*projectionValueArray)[(unsigned)pointIndex] =
				//						ratio Numerical/Analytical		n-points		Observable range	full Integral
				(*projectionValueVector)[(unsigned)pointIndex] * ratioOfIntegrals[combinationIndex] * (double)plotData->GetDataNumber( allCombinations[combinationIndex] ) * range / combination_integral[combinationIndex];

			//	'Fraction of this Combination in Combination 0'
			(*projectionValueArray)[(unsigned)pointIndex] *= combinationWeights[combinationIndex];

			//	This is due to ROOT binning the data into 100 bins by default
			(*projectionValueArray)[(unsigned)pointIndex] /= double(data_binning);
		}

		delete projectionValueVector;

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


//	Generate all unique discrete phase-space boundaries and initialize plotInterval object
vector<PhaseSpaceBoundary*> ComponentPlotter::GeneratePhaseSpaceBoundary( vector<DataPoint*> AllCombinations, vector<double> minimum, vector<double> maximum, vector<double>& plotInterval )
{
	vector<PhaseSpaceBoundary*> discrete_boundaries;
	//      for each combination of components in the PDF phasespace
	for( unsigned int combinationIndex=0; combinationIndex< AllCombinations.size(); ++combinationIndex )
	{
		//      Copy the phasespace so we can safely alter and destroy it
		PhaseSpaceBoundary* new_boundary = new PhaseSpaceBoundary( *full_boundary );

		//      Store the plotInterval (stepsize) for this observable combination
		plotInterval.push_back( ( maximum[combinationIndex] - minimum[combinationIndex] ) / (double)( total_points - 1 ) );

		//      for each discrete observable
		for( unsigned int discreteIndex=0; discreteIndex< discreteNames.size(); ++discreteIndex )
		{
			//      Read this unique combination
			double val = AllCombinations[combinationIndex]->GetObservable( discreteNames[discreteIndex] )->GetValue();
			vector<double> constraint_val; constraint_val.push_back( val );
			//      Store the combination so that we have set the phase-space up
			new_boundary->SetConstraint( discreteNames[discreteIndex], constraint_val,
					AllCombinations[combinationIndex]->GetObservable( discreteNames[discreteIndex] )->GetUnit() );
		}
		//      Store the boundary
		discrete_boundaries.push_back( new_boundary );
	}

	return discrete_boundaries;
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
	//	Not interested in the max/min per subset as ROOT is unable to correctly add histograms over different ranges/binnings (almost an understandable issue)
	//if( plotData->GetBoundary() != NULL )
	//{
	//	if( minimum.size() != observableValues.size() ) minimum.resize( observableValues.size() );
	//	if( maximum.size() != observableValues.size() ) maximum.resize( observableValues.size() );
	//	for( unsigned int i=0; i<= minimum.size(); ++i )
	//	{
	//		minimum[i] = boundary_min;
	//		maximum[i] = boundary_max;
	//	}
	//}
	//else
	//{
	//	cerr << "Cannot normalise without a valid Boundary in the PhaseSpace!!!" << endl;
	//	exit(12375);
	//}

	vector<int> num_observables;
	for( unsigned int i=0; i< observableValues.size(); ++i )
	{
		num_observables.push_back( (int)observableValues[i].size() );
		cout << "Combination: " << i+1 << " has " << num_observables.back() << " events." << endl;
	}


	//vector<double> plotInterval, projectionIntegral;

	//vector<PhaseSpaceBoundary*> discrete_boundaries = this->GeneratePhaseSpaceBoundary( allCombinations, minimum, maximum, plotInterval );


	//if( allCombinations.size() > 1 ) plotInterval.push_back( ( maximum.back() - minimum.back() ) / (double)( total_points - 1 ) );

	//	Now we have:
	//			plot phase-space
	//	Lets plot the projection

	//	Check if the 0'th component has been defined and add it if it's missing.
	vector<string> PDF_Components = plotPDF->PDFComponents();
	if( StringProcessing::VectorContains( PDF_Components, string("0") ) == -1 )
	{
		vector<string> temp_PDF_Components( 1, "0" );
		for(vector<string>::iterator comb_i = PDF_Components.begin(); comb_i != PDF_Components.end(); ++comb_i )
		{
			temp_PDF_Components.push_back( *comb_i );
		}
		PDF_Components = temp_PDF_Components;
	}

	if( onlyZero == true ) PDF_Components = vector<string>(1,"0");

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
	this->WriteOutput( X_values, Y_values, combinationDescriptions );
	here->cd();

	return;
}

//	This routine bins the data for a requested combinationNumber into a TH1, the number of bins is determined through the data_binning int
//
TH1* ComponentPlotter::FormatData( unsigned int combinationNumber )
{
	vector<DataPoint*> wanted_points = data_subsets[combinationNumber];
	TString component_num_str;component_num_str+=combinationNumber;
	TH1* returnable = new TH1D( "Data_For_Component_"+component_num_str, "Data_For_Component_"+component_num_str, data_binning, boundary_min, boundary_max );
	vector<double> wanted_data, wanted_weights;
	ObservableRef* new_ref = new ObservableRef( observableName );
	ObservableRef* weight_ref = new ObservableRef( weightName );
	for( vector<DataPoint*>::iterator point_i = wanted_points.begin(); point_i != wanted_points.end(); ++point_i )
	{
		wanted_data.push_back( (*point_i)->GetObservable( *new_ref )->GetValue() );
		if( weightsWereUsed ) wanted_weights.push_back( (*point_i)->GetObservable( *weight_ref )->GetValue() );
		else wanted_weights.push_back( 1. );
	}
	delete new_ref;

	//	Fill a histogram with the weighted data, with weights normalised to 1
	//	Think of this as multiplying by a very clever version of the number 1 ;)
	//	This gives the correct normalisation
	vector<double>::iterator point_i = wanted_data.begin();
	vector<double>::iterator weight_i = wanted_weights.begin();
	for( ; point_i != wanted_data.end(); ++point_i, ++weight_i )
	{
		returnable->Fill( *point_i, *weight_i/weight_norm );
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

			TTree* this_data = new TTree( TTree_name, TTree_title );

			vector<double>* unnormalised = new vector<double>();

			for( vector<double>::iterator num_i = (*(*Y_values)[componentIndex])[combinationIndex]->begin(); num_i != (*(*Y_values)[componentIndex])[combinationIndex]->end(); ++num_i )
			{
				unnormalised->push_back( (*num_i) * data_binning );
			}

			this->WriteBranch( this_data, "X_data", (*(*X_values)[componentIndex])[combinationIndex] );
			this->WriteBranch( this_data, "Y_data", (*(*Y_values)[componentIndex])[combinationIndex] );
			this->WriteBranch( this_data, "Y_data_unNormalised", unnormalised );

			this_data->Write( "", TObject::kOverwrite );

			delete unnormalised;

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
			data_plot->SetName( data_plot->GetName() + ext );
			data_plot->SetTitle( "" );

			//	Plot the component on this combination and save the TCanvas
			TString Graph_Name("TGraph_");Graph_Name+=componentIndex;
			Graph_Name.Append("_");Graph_Name+=combinationIndex;
			Graph_Name.Append("_");Graph_Name+=plotPDF->GetRandomFunction()->Rndm();
			TGraph* data_graph = new TGraph( total_points,  &((*(*(*X_values)[componentIndex])[combinationIndex])[0]),  &((*(*(*Y_values)[componentIndex])[combinationIndex])[0]) );
			data_graph->SetTitle( "" ); data_graph->SetName( Graph_Name );

			TString Canvas_Name("TCanvas_");Canvas_Name+=componentIndex;
			Canvas_Name.Append("_");Canvas_Name+=combinationIndex;
			Canvas_Name.Append("_");Canvas_Name+=plotPDF->GetRandomFunction()->Rndm();
			TCanvas* c1 = new TCanvas( Canvas_Name, "", 1680, 1050 );

			data_plot->Draw();
			data_graph->Draw("PL9 SAME");
			c1->Update();
			data_graph->GetYaxis()->SetRangeUser( 0., data_graph->GetYaxis()->GetXmax() );
			c1->Update();
			c1->Write();

			componentDir->cd();
		}

		PlotDirectory->cd();
	}

	//	Starting at the top of the file again
	PlotFile->cd();

	gDirectory->cd( pdfStr );

	if( gDirectory->GetDirectory( "overlay_graphs" ) == 0 )	gDirectory->mkdir( "overlay_graphs" );
	gDirectory->cd( "overlay_graphs" );

	//	Overlay all components for each combination
	//
	//	This is painful as your looping over a the upper 2D of a 3D vector inside out...
	//	Sorry but I'm not rewriting the whole class or inverting the vector of vectors as we will just get more mistakes
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
			data_graph->SetName( data_name ); data_graph->SetTitle( data_name );
			data_graph->SetLineColor( (Color_t)(componentIndex+1) );
			data_graph->SetMarkerColor( (Color_t)(componentIndex+1) );
			these_components.push_back( data_graph );
		}
		total_components.push_back( these_components );

		TH1* data_plot = this->FormatData( combinationIndex );
		TString ext("_"); ext+=combinationIndex; ext.Append("_"); ext+=plotPDF->GetRandomFunction()->Rndm();
		data_plot->SetName( data_plot->GetName() + ext );
		data_plot->SetTitle( "" );

		//	Object to hold the binned data for the life of this plot
		binned_data.push_back( new TGraphErrors( data_plot ) );

		binned_data.back()->Write("",TObject::kOverwrite);

		string desc;

		if( combinationIndex> 0 ) desc = CombinationDescriptions[combinationIndex-1];
		if( combinationIndex == 0 )
		{
			for( vector<string>::iterator desc_i = CombinationDescriptions.begin(); desc_i != CombinationDescriptions.end(); ++desc_i )
			{
				desc.append( *desc_i );
			}
		}

		//	For the moment haven't decided if I should pass the global config to ALL sub plots, I shall get user input on this
		CompPlotter_config* temp = new CompPlotter_config();
		temp->logY = logY;

		//	Static function so has to be told everything about what you want to plot!
		this->OutputPlot( binned_data.back(), these_components, observableName, desc, plotData->GetBoundary(), plotPDF->GetRandomFunction(), temp );

		if( this_config != NULL )
		{
			if( combinationIndex == 0 && this_config->CalcChi2 == true )
			{
				cout << endl;
				cout << "Calculating Chi^2:" << endl;

				TF1* fitting_function = new TF1( "total_PDF", this, boundary_min, boundary_max, 1, "" );        //      I REFUSE to pass the class name
				chi2 = binned_data.back()->Chisquare( fitting_function );

				N = plotData->GetDataNumber();
				double n = (double) plotPDF->GetPhysicsParameters()->GetAllFloatNames().size();
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
				   TCanvas* Chi2test = new TCanvas( Chi2Name, Chi2Title, 1680, 1050 );
				   binned_data.back()->Draw("AP");
				   tf1_graph->Draw("C");
				   Chi2test->Print(Chi2Name+".pdf");
				   Chi2test->Write("",TObject::kOverwrite);
				   */
			}
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
		TGraph* output_graph = new TGraph( (Int_t)X_val[0].size(), &(final_X_val[0]), &(final_Y_val[0]) );
		output_graph->SetName( TGraphName );
		output_graph->SetTitle("");
		output_graph->SetLineColor( this_component[0]->GetLineColor() );
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
	TGraphErrors* output_graph = new TGraphErrors( (Int_t)X_val[0].size(), &(final_X_val[0]), &(final_Y_val[0]), &(final_X_err[0]), &(final_Y_err[0]) );
	output_graph->SetName( TGraphErrorsName );
	output_graph->SetTitle("");

	output_graph->SetLineColor( input[0]->GetLineColor() );
	output_graph->SetLineStyle( input[0]->GetLineStyle() );
	output_graph->SetMarkerColor( input[0]->GetMarkerColor() );
	output_graph->SetMarkerStyle( input[0]->GetMarkerStyle() );

	return output_graph;
}

//	Plot all components on this combinations and print and save the canvas
void ComponentPlotter::OutputPlot( TGraphErrors* input_data, vector<TGraph*> input_components, string observableName, string CombinationDescription, PhaseSpaceBoundary* total_boundary, TRandom* rand, CompPlotter_config* conf )
{
	if( rand == NULL ) rand = gRandom;
	TString TCanvas_Name("Overlay_"+observableName+"_"+CombinationDescription+"_");TCanvas_Name+=rand->Rndm();

	TCanvas* c1 = new TCanvas( TCanvas_Name, "", 1680, 1050 );

	TString plotTitle;

	vector<int> Style_Key, Color_Key, Width_Key;

	vector<string> component_names;

	double X_min=-999, X_max=-999;
	double Y_min=-999, Y_max=-999;

	TString X_Title, Y_Title;

	bool logy=false;

	double final_chi2=-999;

	if( conf != NULL )
	{
		if( conf->logY )
		{
			logy=true;
			c1->SetLogy( true );
			//c1->Update();
		}
		plotTitle = conf->PlotTitle;
		Style_Key = conf->style_key;
		Color_Key = conf->color_key;
		Width_Key = conf->width_key;
		component_names = conf->component_names;
		X_min = conf->xmin;
		X_max = conf->xmax;
		Y_min = conf->ymin;
		Y_max = conf->ymax;
		X_Title = conf->xtitle;
		Y_Title = conf->ytitle;
		final_chi2 = conf->Chi2Value;
	}

	input_data->SetTitle( plotTitle );
	input_data->Draw("AP9");
	c1->Update();

	if( X_min <= -999 ) X_min = total_boundary->GetConstraint( observableName )->GetMinimum();
	if( X_max <= -999 ) X_max = total_boundary->GetConstraint( observableName )->GetMaximum();
	if( Y_min <= -999 ) Y_min = logy==true?0.5:0.;
	if( Y_max <= -999 ) Y_max = input_data->GetYaxis()->GetXmax();

	if( StringProcessing::is_empty( X_Title ) ) X_Title = EdStyle::GetParamRootName( observableName ) + " " + EdStyle::GetParamRootUnit( observableName );
	if( StringProcessing::is_empty( Y_Title ) ) Y_Title = "Events";

	input_data->GetYaxis()->SetRangeUser( Y_min, Y_max );
	input_data->GetYaxis()->SetTitle( Y_Title );
	input_data->GetXaxis()->SetRangeUser( X_min, X_max );
	input_data->GetXaxis()->SetTitle( X_Title );

	c1->Update();

	TLegend* leg = EdStyle::LHCbLegend();

	leg->AddEntry( input_data, "Data", "pl" );

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
			(*comp_i)->Draw("L9");
		}

		if( !component_names.empty() )
		{
			if( num < component_names.size() )
			{
				leg->AddEntry( (*comp_i), TString(component_names[num]), "l" );
			}
			else
			{
				leg->AddEntry( (*comp_i), "Unnamed", "l" );
			}
		}
	}

	if( final_chi2 > 0 )
	{
		TString Chi2Text( "#chi^{2}/ndof : " );
		stringstream chi_stream; chi_stream << setprecision(4) << final_chi2;
		Chi2Text.Append( chi_stream.str() );
		leg->AddEntry( (TObject*)NULL, Chi2Text, "" );
	}

	c1->Update();

	if( !component_names.empty() ) leg->Draw();

	c1->Update();

	c1->Write("",TObject::kOverwrite);

	TString Clean_Description = StringProcessing::Clean( CombinationDescription.c_str() );

	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL, *nullbuf=NULL;
	ofstream filestr;
	filestr.open ("/dev/null");
	//      If the user wanted silence we point the Std Output Streams to the oblivion of NULL
	cout_bak = cout.rdbuf();
	cerr_bak = cerr.rdbuf();
	clog_bak = clog.rdbuf();
	nullbuf = filestr.rdbuf();
	freopen("/dev/null", "w", stderr);
	cout.rdbuf(nullbuf);
	cerr.rdbuf(nullbuf);
	clog.rdbuf(nullbuf);

	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".C") );
	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".pdf") );
	c1->Print( TString("Overlay_"+observableName+"_"+Clean_Description+".png") );

	cout.rdbuf( cout_bak );
	cerr.rdbuf( cerr_bak );
	clog.rdbuf( clog_bak );
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

//Run the statistics functions
vector<double> ComponentPlotter::GetStatistics( string ObservableName, double & Minimum, double & Maximum, int & BinNumber )
{
	//Make a vector of the values found for the observable
	vector<double> observableValues;
	for ( int dataIndex = 0; dataIndex < plotData->GetDataNumber(); ++dataIndex )
	{
		observableValues.push_back( plotData->GetDataPoint(dataIndex)->GetObservable(ObservableName)->GetValue() );
	}

	//Find out number and range of data points needed
	Maximum = StatisticsFunctions::Maximum(observableValues);
	Minimum = StatisticsFunctions::Minimum(observableValues);
	BinNumber = StatisticsFunctions::OptimumBinNumber(observableValues);

	cout << endl << "Optimally:\t" << "Max: " << Maximum << "\tMin: " << Minimum << "\tBinning: " << BinNumber << endl;

	return observableValues;
}

//Return the values tracing the PDF projection in the Observable of choice
vector<double>* ComponentPlotter::ProjectObservableComponent( DataPoint* InputPoint, string ObservableName, double Minimum, int PlotNumber, double PlotInterval, string component )
{
	//	Move initializer(s) outside of for loop
	double integralvalue=-1., observableValue=-1.;
	Observable* oldObs=NULL; Observable* newObs=NULL;

	//Find the value of the observable projection at each data point
	vector<double>* pointValues = new vector<double>();

	cout << ObservableName << ": " << Minimum << " <-> " << Minimum + ( PlotInterval * (PlotNumber-1) ) << endl;

	//	This class object has been created to speed up the communication between this class and the basePDF as it may pass through several PDF wrappers
	ComponentRef* comp_obj = new ComponentRef( component );

	//	Step over the whole observable range
	for (int pointIndex = 0; pointIndex < PlotNumber; ++pointIndex )
	{
		//	Inform the user of how far we have got :D
		cout << setprecision(3) << setw(4) << (pointIndex+1.)/PlotNumber*100 << "\%\tComplete\r\r" << flush;

		//	Value of Observable we want to evaluate for this step
		observableValue = Minimum + ( PlotInterval * pointIndex );

		//	Set the value of Observable to the new step value
		//	All other CONTINUOUS observables set to average of range
		oldObs = InputPoint->GetObservable(ObservableName);
		newObs = new Observable( ObservableName, observableValue, 0, oldObs->GetUnit() );
		InputPoint->SetObservable( ObservableName, newObs );
		delete newObs; newObs = NULL;

		//	perform actual evaluation of the PDF with the configuration in the InputPoint, in the whole boundary with the given name and for selected component
		integralvalue = pdfIntegrator->ProjectObservable( InputPoint, full_boundary, ObservableName, comp_obj );

		pointValues->push_back( integralvalue );
	}
	delete comp_obj;

	this->Sanity_Check( pointValues, component );

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


//	Calculate the statistics and return each data subset corresponding to each combination that has been defined
vector<vector<double> > ComponentPlotter::GetCombinationStatistics( string ObservableName, vector<double>& Min, vector<double>& Max, vector<int>& BinNum )
{
	cout << "Optimally:" << endl << endl;
	//	Object to contain the output data
	vector<vector<double> > temp_discrete_observableValues;

	//	If we have 'no discrete combinations' then get the whole dataset for this pdf
	vector<DataPoint*> whole_dataset;

	cout << plotData->GetDataNumber() << endl;

	for( int i=0; i< plotData->GetDataNumber(); ++i )
	{
		whole_dataset.push_back( plotData->GetDataPoint( i ) );
	}

	data_subsets.push_back( whole_dataset );

	vector<double> discrete_observableValues_i;
	//	Get all values of the wanted parameter from this subset
	for( unsigned int j=0; j< whole_dataset.size(); ++j )
	{
		discrete_observableValues_i.push_back( whole_dataset[j]->GetObservable( ObservableName )->GetValue() );
	}

	//	Save the min/max and optimal bin number for this subset
	Max.push_back( StatisticsFunctions::Maximum( discrete_observableValues_i ) );
	Min.push_back( StatisticsFunctions::Minimum( discrete_observableValues_i ) );
	BinNum.push_back( StatisticsFunctions::OptimumBinNumber( discrete_observableValues_i ) );

	cout << "Combination: 1 " << "Min: " << Min.back() << "\tMax: " << Max.back() << "\tBinNum: " << BinNum.back() << endl;

	//	store all of the values of this observable
	temp_discrete_observableValues.push_back( discrete_observableValues_i );

	//	If we have discrete combinations then we need to loop over each set seperatley to split up the data and information about the data
	if( allCombinations.size() > 1 )
	{
		//	For all discrete combinations that have been calculated elsewhere
		for( unsigned int combinationIndex=0; combinationIndex< allCombinations.size(); ++combinationIndex )
		{
			cout << "Combination: " << combinationIndex+2 << " ";
			//	For this combination work out what values correspond to the discrete parameters
			vector<double> values_for_this_combination;
			for( unsigned int paramIndex=0; paramIndex< discreteNames.size(); ++paramIndex )
			{
				values_for_this_combination.push_back( allCombinations[combinationIndex]->GetObservable( discreteNames[paramIndex] )->GetValue() );
			}

			//		Get all data which falls within this paramreter set
			//vector<DataPoint*> discrete_set_i = plotData->GetDiscreteSubSet( discreteNames, values_for_this_combination );
			data_subsets.push_back( plotData->GetDiscreteSubSet( discreteNames, values_for_this_combination ) );

			//print( discreteNames );
			//print( values_for_this_combination );

			vector<double> discrete_observableValues_i;
			//	Get all values of the wanted parameter from this subset
			for( unsigned int j=0; j< data_subsets.back().size(); ++j )
			{
				discrete_observableValues_i.push_back( (data_subsets.back())[j]->GetObservable( ObservableName )->GetValue() );
			}

			//print( discrete_observableValues_i );

			//	Save the min/max and optimal bin number for this subset
			Max.push_back( StatisticsFunctions::Maximum( discrete_observableValues_i ) );
			Min.push_back( StatisticsFunctions::Minimum( discrete_observableValues_i ) );
			BinNum.push_back( StatisticsFunctions::OptimumBinNumber( discrete_observableValues_i ) );

			cout << "Min: " << Min.back() << "\tMax: " << Max.back() << "\tBinNum: " << BinNum.back() << endl;

			//	store all of the values of this observable
			temp_discrete_observableValues.push_back( discrete_observableValues_i );
		}
	}

	return temp_discrete_observableValues;
}

//Return a list of data points
//Each should take the data average value of each continuous observable
//Each should represent one combination of possible discrete values
vector<DataPoint*> ComponentPlotter::GetDiscreteCombinations( vector<double>& DataPointWeights, vector<string>& DataPointDescriptions )
{
	//Calculate all possible combinations of discrete observables
	vector<string> allNames = full_boundary->GetAllNames();
	vector<vector<double> > discreteValues;
	vector< vector<double> > discreteCombinations = StatisticsFunctions::DiscreteCombinations( &allNames, plotData->GetBoundary(), discreteNames, continuousNames, discreteValues );

	//Initialise the data averaging
	vector<double> continuousSums, continuousStdevs;
	vector<long> combinationCounts;
	for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
	{
		continuousSums.push_back(0.0);
		continuousStdevs.push_back(0.0);
	}
	for (unsigned int combinationIndex = 0; combinationIndex < discreteCombinations.size(); ++combinationIndex )
	{
		combinationCounts.push_back(0);
	}

	//Examine the data set. Find the average value for each continuous observable, and the weight for each discrete combination
	for (int dataIndex = 0; dataIndex < plotData->GetDataNumber(); ++dataIndex )
	{
		DataPoint* readDataPoint = plotData->GetDataPoint(dataIndex);

		//Sum the continuous values, in preparation for taking the average
		for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
		{
			continuousSums[continuousIndex] += readDataPoint->GetObservable( continuousNames[continuousIndex] )->GetValue();
		}

		//Calculate the index for the discrete combination, and increment the corresponding count
		int combinationIndex = 0;
		int incrementValue = 1;
		for (int discreteIndex = int(discreteNames.size() - 1); discreteIndex >= 0; discreteIndex-- )
		{
			double currentValue = readDataPoint->GetObservable( discreteNames[unsigned(discreteIndex)] )->GetValue();

			for (unsigned int valueIndex = 0; valueIndex < discreteValues[unsigned(discreteIndex)].size(); ++valueIndex )
			{
				if ( fabs( discreteValues[unsigned(discreteIndex)][valueIndex] - currentValue ) < DOUBLE_TOLERANCE )
				{
					combinationIndex += ( incrementValue * int(valueIndex) );
					incrementValue *= int(discreteValues[unsigned(discreteIndex)].size());
					break;
				}
			}
		}
		++combinationCounts[unsigned(combinationIndex)];
	}

	//Calculate averages and weights
	vector<double> combinationWeights;
	double dataNumber = (double)plotData->GetDataNumber();
	for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
	{
		continuousSums[continuousIndex] /= dataNumber;
		for( int dataIndex = 0; dataIndex < dataNumber; ++dataIndex )
		{
			DataPoint* readDataPoint = plotData->GetDataPoint(dataIndex);
			double val = readDataPoint->GetObservable( continuousNames[continuousIndex] )->GetValue();
			continuousStdevs[continuousIndex]+=( val - continuousSums[continuousIndex] )*( val -continuousSums[continuousIndex] );
		}
		continuousStdevs[continuousIndex] = sqrt( continuousStdevs[continuousIndex] );
		continuousStdevs[continuousIndex] /= ((double)dataNumber-1.);
	}
	for (unsigned int combinationIndex = 0; combinationIndex < discreteCombinations.size(); ++combinationIndex )
	{
		combinationWeights.push_back( (double)combinationCounts[combinationIndex] / (dataNumber) );
	}

	//Create the data points to return
	vector<DataPoint*> newDataPoints;
	vector<string> allDescriptions;

	for (unsigned int combinationIndex = 0; combinationIndex < discreteCombinations.size(); ++combinationIndex )
	{
		DataPoint* templateDataPoint = new DataPoint( full_boundary->GetAllNames() );

		for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
		{
			Observable* oldValue = templateDataPoint->GetObservable( continuousNames[continuousIndex] );
			Observable* newValue = new Observable( oldValue->GetName(), continuousSums[continuousIndex], continuousStdevs[continuousIndex], oldValue->GetUnit() );
			templateDataPoint->SetObservable( continuousNames[continuousIndex], newValue );
			delete newValue;
		}

		string description = "(";

		//Output the discrete values for this combination
		for (unsigned int discreteIndex = 0; discreteIndex < discreteNames.size(); ++discreteIndex )
		{
			//Set the data point
			Observable* oldValue = templateDataPoint->GetObservable( discreteNames[discreteIndex] );

			Observable* newValue = new Observable( oldValue->GetName(), discreteCombinations[combinationIndex][discreteIndex],
					0, oldValue->GetUnit() );
			templateDataPoint->SetObservable( discreteNames[discreteIndex], newValue );
			delete newValue;

			stringstream valuestream;
			valuestream << discreteCombinations[combinationIndex][discreteIndex];
			string value = valuestream.str();

			description.append(discreteNames[discreteIndex]);
			description.append("=");
			description.append(value);

			if( discreteIndex != discreteNames.size() - 1 ) description.append(")(");
		}

		description.append(")");
		allDescriptions.push_back(description);
		newDataPoints.push_back(templateDataPoint);
	}

	cout << "Unique Combinations: " << newDataPoints.size() << endl << endl;

	//Output the results
	DataPointDescriptions = allDescriptions;
	DataPointWeights = combinationWeights;
	return newDataPoints;
}

void ComponentPlotter::SetWeightsWereUsed( string input )
{
	weightsWereUsed = true;
	weightName = input;

	ObservableRef* weight_ref = new ObservableRef( weightName );
	for( int i=0 ; i < plotData->GetDataNumber(); ++i )
	{
		if( weightsWereUsed ) wanted_weights.push_back( plotData->GetDataPoint( i )->GetObservable( *weight_ref )->GetValue() );
		else wanted_weights.push_back( 1. );
	}
	delete weight_ref;

	weight_norm=0.;
	if( weightsWereUsed )
	{
		for( vector<double>::iterator w_i = wanted_weights.begin(); w_i != wanted_weights.end(); ++w_i )
		{
			weight_norm+=*w_i;
		}
		weight_norm /= plotData->GetDataNumber();
	}
	else
	{
		weight_norm=1.;
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
	(void) p;
	double fixed_obs_val = x[0];

	Observable* oldObs=NULL;
	Observable* newObs=NULL;

	double integral_value=0.;

	ComponentRef* comp_obj = new ComponentRef( "0" );


	unsigned int comb_num=0;

	for( vector<DataPoint*>::iterator comb_i = allCombinations.begin();  comb_i != allCombinations.end(); ++comb_i, ++comb_num )
	{
		DataPoint* InputPoint = new DataPoint( *( (*comb_i) ) );

		oldObs = InputPoint->GetObservable( observableName );
		newObs = new Observable( observableName, fixed_obs_val, 0, oldObs->GetUnit() );
		InputPoint->SetObservable( observableName, newObs );
		//delete oldObs;

		double temp_integral_value = pdfIntegrator->ProjectObservable( InputPoint, full_boundary, observableName, comp_obj );

		temp_integral_value = temp_integral_value * ratioOfIntegrals[comb_num] * (double)plotData->GetDataNumber() * (boundary_max-boundary_min) / combination_integral[comb_num];

		temp_integral_value = temp_integral_value * ( combinationWeights[comb_num] / double(data_binning) );

		integral_value += temp_integral_value;

		//delete InputPoint;
	}

	cout << "Value At: " << setprecision(3) << x[0] << "\tis:\t" << setprecision(4) << integral_value << "\r\r";

	return integral_value;
}

pair<double,double> ComponentPlotter::GetChi2Numbers()
{
	return make_pair( chi2, N );
}

