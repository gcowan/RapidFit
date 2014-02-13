
#include "TH1D.h"
#include "TCanvas.h"
#include "THn.h"
#include "MultiDimChi2.h"
#include "PhaseSpaceBoundary.h"
#include "IDataSet.h"
#include <vector>
#include "StringProcessing.h"
#include "IPDF.h"
#include <string>
#include <cmath>
#include "ClassLookUp.h"
#include "PDFWithData.h"
#include "RapidFitIntegrator.h"
#include "RapidFitIntegratorConfig.h"

using namespace::std;

MultiDimChi2::MultiDimChi2( vector<PDFWithData*> allObjects, PhaseSpaceBoundary* thisBound, vector<string> wantedObservables ) :
	internalHisto(NULL), allBinCenters(NULL) 
{
	//	Populate PDF, DataSets and PhaseSpaceBoundaries to be used in this analysis
	cout << "MultiDimChi2 Hello!" << endl << endl;
	cout << "initializing objects" << endl;
	this->populateAllObjects( allObjects );
	this->ConstructBoundaries( thisBound, wantedObservables );

	cout << "Deciding binning and constructing THnD" << endl;
	//	Construct the Histogram object which contains all of the Data of Interest
	this->ConstructInternalHisto( wantedObservables, thisBound );

	cout << "Calculating Bin Centers" << endl;
	//      Construct the Central coordinates of every bin in each axis, this will be the basis for generating the corrdinates to run the Chi2 test over
	this->ConstructBinCenters();

	cout << "Calculating the Numerical vs Analytical Ratios" << endl;
	//	Construct the Ratios of the Numerical vs analytical Integrals (condensing the differenced between the PDF numerator and denominator to be a single number but what else can we do...)
	this->ConstructIntegralsRatios( wantedObservables );

	cout << "Constructing all of the Coordinates to perform a Chi2 test over" << endl;
	//	Construct the full set of coordinates that the Chi2 test is to be run over
	this->ConstructAllCoordinates();

	cout << "Finished Initializing!" << endl << endl;
}

void MultiDimChi2::populateAllObjects( vector<PDFWithData*> allObjects )
{
	allPDFData = allObjects;	//	Copy Reference, don't take ownership!
	for( unsigned int i=0; i< allObjects.size(); ++i )
	{
		allDataSets.push_back( allObjects[i]->GetDataSet() );	//      Copy Reference, don't take ownership!
		allPDFs.push_back( ClassLookUp::CopyPDF( allObjects[i]->GetPDF() ) );	// 	More than likely going to cause something to change so lets make it on our copy, just incase
		RapidFitIntegratorConfig* thisConf = allPDFs.back()->GetPDFIntegrator()->GetIntegratorConfig();
		thisConf->useGSLIntegrator = true;
		thisConf->numThreads = 16;
		allPDFs.back()->GetPDFIntegrator()->SetUpIntegrator( thisConf );
	}
}

void MultiDimChi2::ConstructIntegralsRatios( vector<string> wanted_params )
{
	(void) wanted_params;
	ratioOfIntegrals.clear();
	combinationIntegrals.clear();
	weightNorms.clear();

	for( unsigned int PDFNum=0; PDFNum< allPDFs.size(); ++PDFNum )
	{
		cout << "\tCalculating Ratios for PDF: " << PDFNum+1 << endl;
		IPDF* thisPDF = allPDFs[PDFNum];
		RapidFitIntegrator* pdfIntegrator = thisPDF->GetPDFIntegrator();
		IDataSet* thisData = allDataSets[PDFNum];
		vector<double> thisratioOfIntegrals;
		vector<double> thisCombination_integral;
		PhaseSpaceBoundary* thisBound = allBoundaries[PDFNum];

		thisPDF->ChangePhaseSpace( thisBound );

		vector<DataPoint*> allCombinations = thisBound->GetDiscreteCombinations();

		//thisBound->Print();


		for( unsigned int i=0; i< allCombinations.size(); ++i )
		{
			thisratioOfIntegrals.push_back( 1. );
			thisCombination_integral.push_back( 1. );

			/*			cout << "\t\tCalculating Ratios for Combination: " << i+1 << " of " << allCombinations.size() << endl;
			//allCombinations[i]->Print();
			//
			//   if( debug != NULL )
			//   {
			//   if( debug->DebugThisClass( "MultiDimChi2" ) )
			//   {
			//   cout << "MultiDimChi2:: Calculating Test integrals:\t" << i << endl;
			//   }
			//   }
			pdfIntegrator->ForceTestStatus( false );

			vector<string> doNotList = thisPDF->GetDoNotIntegrateList();
			//doNotList = StringProcessing::CombineUniques( doNotList, wanted_params );
			double thisIntegral = 0.;
			try
			{
			//allCombinations[i]->SetPhaseSpaceBoundary( thisBound );
			//cout << thisBound->DiscreteDescription( allCombinations[i] ) << endl;
			thisIntegral = pdfIntegrator->NumericallyIntegrateDataPoint( allCombinations[i], thisBound, doNotList );
			}
			catch(...)
			{
			thisIntegral = 1.;
			cout << endl << "CANNOT PROPERLY NORMALISE WHOLE PDF, THIS WILL LEAD TO NORMALISATION ISSUES OVER THE WHOLE PDF" << endl << endl;
			}
			//
			//   if( debug != NULL )
			//   {
			//   if( debug->DebugThisClass( "MultiDimChi2" ) )
			//   {
			//   cout << "MultiDimChi2:: Finished Inetgral" << endl;
			//   }
			//   }

			thisCombination_integral.push_back( thisIntegral );

			if( thisPDF->GetNumericalNormalisation() == true ) thisratioOfIntegrals.push_back( 1. );
			//else
			//{
			//	thisratioOfIntegrals.push_back( pdfIntegrator->GetRatioOfIntegrals() );
			//}
			thisratioOfIntegrals.push_back( 1. );

*/
		}

		/*
		   if( debug != NULL )
		   {
		   if( debug->DebugThisClass( "MultiDimChi2" ) )
		   {
		   cout << "MultiDimChi2:: Calculated the ratio of Integrals (analytic/numeric) to be:" << endl;
		   for( unsigned int i=0; i< thisratioOfIntegrals.size(); ++i )
		   {
		   cout << thisratioOfIntegrals[i] << endl;
		   }
		   }
		   }*/
		if( thisratioOfIntegrals.size() == 0 || thisratioOfIntegrals.empty() ) thisratioOfIntegrals =vector<double> ( 1, 1. );

		ratioOfIntegrals.push_back( thisratioOfIntegrals );
		//combinationIntegrals.push_back( thisCombination_integral );

		double thisWeight_norm = 0.;
		double weight_sum = thisData->GetSumWeights();

		thisWeight_norm = weight_sum / thisData->GetDataNumber(NULL);

		/*
		   if( debug != NULL )
		   {
		   if( debug->DebugThisClass( "MultiDimChi2" ) )
		   {
		   cout << "Weight Norm: " << weight_norm << endl;
		   }
		   } */
		weightNorms.push_back( thisWeight_norm );

		cout << "Alpha: " << thisData->GetAlpha() << endl;
		//weightNorms.back() /= thisData->GetAlpha();

		pdfIntegrator->ForceTestStatus( true );
	}
}

void MultiDimChi2::ConstructInternalHisto( vector<string> wantedObservables, PhaseSpaceBoundary* thisBound )
{
	nDim=0;
	x_min.clear(); x_max.clear(); goodObservables.clear(); x_bins.clear();

	cout << "\tDecicing Binning" << endl;
	data_binning = 1;
	unsigned int x_binning=0;
	for( unsigned int i=0; i< wantedObservables.size(); ++i )
	{
		IConstraint* thisConstraint = thisBound->GetConstraint( wantedObservables[i] );
		if( !thisConstraint->IsDiscrete() )
		{
			++nDim;
			x_min.push_back( thisConstraint->GetMinimum() );
			x_max.push_back( thisConstraint->GetMaximum() );
			goodObservables.push_back( wantedObservables[i] );
			if( wantedObservables[i] == "time" ) x_binning = 4.;
			else x_binning = 4;
			//x_binning = 2;
			x_bins.push_back( x_binning );
			data_binning *= x_binning;
			ThisObsBinning* thisDimension = new ThisObsBinning();
			thisDimension->ObservableName = wantedObservables[i];
			thisDimension->thisMin = x_min.back();
			thisDimension->thisMax = x_max.back();
			thisDimension->theseBins = x_binning;
			thisDimension->thisStepSize = (thisDimension->thisMax-thisDimension->thisMin)/((double)thisDimension->theseBins);
			theseDimensions.push_back( thisDimension );
		}
	}

	cout << "\tConstructing THnD" << endl;
	if( internalHisto != NULL ) delete internalHisto;

	cout << "\t\tinternal_Chi2nDim\t" << nDim << "\t" << x_bins[0] << "\t" << x_min[0] << "\t" << x_max[0] << endl;
	internalHisto = new THnD( "internal_Chi2nDim", "internal_Chi2nDim", (int)nDim, &(x_bins[0]), &(x_min[0]), &(x_max[0]) );

	vector<double> thisValues( goodObservables.size(), 0. );
	double thisWeight = 0.;
	DataPoint* thisPoint=NULL;

	cout << "\tPopulating THnD" << endl;
	for( unsigned int i=0; i< allDataSets.size(); ++i )
	{
		for( unsigned int j=0; j< (unsigned)allDataSets[i]->GetDataNumber(); ++j )
		{
			thisPoint = allDataSets[i]->GetDataPoint( j );
			for( unsigned int k=0; k< goodObservables.size(); ++k )
			{
				thisValues[k] = thisPoint->GetObservable( goodObservables[k] )->GetValue();
			}
			thisWeight = thisPoint->GetEventWeight();
			internalHisto->Fill( &(thisValues[0]), thisWeight );
		}
	}

	for( unsigned int i=0; i< theseDimensions.size(); ++i )
	{
		TH1D* h1 = internalHisto->Projection( i, "E" );
		TString name("c");
		name+=i;
		TCanvas* c1 = new TCanvas( name, name );
		h1->Draw("");
		c1->Print("Output.pdf");
	}
}

void MultiDimChi2::ConstructBinCenters()
{
	double thisMin=0.;
	double thisMax=0.;
	int theseBins=0;
	double thisCenter=0.;
	double thisStepSize=0.;
	double numSteps=0.;
	for( unsigned int i=0; i< nDim; ++i )
	{
		vector<double> thisObservableBins;
		thisMin = x_min[i];
		thisMax = x_max[i];
		theseBins = x_bins[i];
		thisStepSize = (thisMax-thisMin)/((double)theseBins);

		for( unsigned int j=1; j<= (unsigned)theseBins; ++j )
		{
			numSteps=(double)j; numSteps-=0.5;
			thisCenter = thisMin+numSteps*thisStepSize;

			thisObservableBins.push_back( thisCenter );
		}
		theseDimensions[i]->binCenters = thisObservableBins;
	}
}

void MultiDimChi2::PerformMuiltDimTest()
{
	cout << "MultiDimChi2: About to Perform Chi2 Calculation" << endl;
	PhaseSpaceBoundary* thisBoundary=NULL;
	(void) thisBoundary;

	vector<double> expected_events( allBinCenters->size(), 0. );

	vector<double> observed_events( allBinCenters->size(), 0. );

	unsigned int histo_binNum=0;
	vector<double> thisBinCenter( (*allBinCenters)[0].size(), 0. );

	cout << "Looping over: " << allBinCenters->size()+1 << " coordianates!" << endl;
	for( unsigned int binNum=0; binNum< allBinCenters->size(); ++binNum )
	{
		thisBinCenter = allBinCenters->at( binNum );

		histo_binNum = (unsigned)internalHisto->GetBin( &(thisBinCenter[0]) );

		cout << "Determining Number of Observed Events in Bin: " << histo_binNum << endl;

		observed_events[binNum] = internalHisto->GetBinContent( histo_binNum );

		cout << "Calculating Number of Expected Events in Bin: " << histo_binNum << endl;

		expected_events[binNum] = this->CalculateTotalExpected( thisBinCenter );

		cout << "O: " << observed_events[binNum] << "  E: " << expected_events[binNum] << endl;
		//exit(0);
	}

	cout << "Finished Looping over all coordinates!" << endl;

	double TotalChi2 = this->CalcChi2( expected_events, observed_events );

	cout << "The Total Chi2 = " << TotalChi2 << endl << endl;

	return;
}

double MultiDimChi2::CalcChi2( vector<double> expected_events, vector<double> observed_events )
{
	double Chi2Value=0.;
	double thisChi2=0.;
	for( unsigned int i=0; i< expected_events.size(); ++i )
	{
		cout << "O: " << observed_events[i] << "  E: " << expected_events[i] << endl;
		thisChi2 = expected_events[i] - observed_events[i] + observed_events[i] * log( observed_events[i] / expected_events[i] );
		/*
		//thisChi2 = expected_events[i] - observed_events[i] + observed_events[i]*log( observed_events[i] / expected_events[i] );
		vector<double> thisBinCenter = allBinCenters->at( i );
		int histo_binNum = (unsigned)internalHisto->GetBin( &(thisBinCenter[0]) );
		thisChi2 = ( observed_events[i] - expected_events[i] ) / internalHisto->GetBinError( histo_binNum );
		cout << thisChi2*thisChi2 << endl;
		*/
		if( !std::isnan(thisChi2) && !std::isinf(thisChi2) )
		{
			Chi2Value += thisChi2;
		}
	}
	return 2.*Chi2Value;
}

double MultiDimChi2::CalculateTotalExpected( vector<double> thisBinCenter )
{

	double result_AllPDFs = 0.;

	cout << "Calculating Expected Number of Events for this Bin" << endl;

	IPDF* thisPDF=NULL;
	IDataSet* thisDataSet=NULL;
	DataPoint* thisDataPoint=NULL;

	vector<ObservableRef> theseDim;
	vector<string> theseDimName;
	for( unsigned int i=0; i< theseDimensions.size(); ++i )
	{
		theseDim.push_back( ObservableRef( theseDimensions[i]->ObservableName ) );
		theseDimName.push_back( theseDimensions[i]->ObservableName );
	}

	for( unsigned int PDFNum=0; PDFNum< allPDFs.size(); ++PDFNum )
	{
		thisPDF = allPDFs[PDFNum];
		thisDataSet = allDataSets[PDFNum];
		RapidFitIntegrator* thisPDFIntegrator = thisPDF->GetPDFIntegrator();

		thisPDFIntegrator->SetPDF( thisPDF );
		thisPDFIntegrator->ForceTestStatus( false );

		vector<string> doNotIntegrate = thisPDF->GetDoNotIntegrateList();

		doNotIntegrate = StringProcessing::CombineUniques( doNotIntegrate, theseDimName );

		PhaseSpaceBoundary* thisPhaseSpace = new PhaseSpaceBoundary( *allBoundaries[PDFNum] );
		PhaseSpaceBoundary* thisPhaseSpace2 = new PhaseSpaceBoundary( *allBoundaries[PDFNum] );

		vector<DataPoint*> theseDataPoints;
		vector<DataPoint*> tempPoints = thisPhaseSpace->GetDiscreteCombinations();

		for( unsigned int i=0; i< tempPoints.size(); ++i ) theseDataPoints.push_back( new DataPoint( *tempPoints[i] ) );

		double thisResult=0.;

		cout << "\tI have " << theseDataPoints.size() << " integrals to perform for PhaseSpace: " << PDFNum+1 << endl;

		for( unsigned int combinationNum = 0; combinationNum< theseDataPoints.size(); ++combinationNum )
		{
			thisDataPoint = new DataPoint( *theseDataPoints[combinationNum] );

			for( unsigned int i=0; i< theseDim.size(); ++i )
			{
				Observable* newObs = new Observable( string(theseDim[i]), thisBinCenter[i], "noUnits_Chi2" );
				thisDataPoint->SetObservable( theseDim[i], newObs );
				delete newObs;
				double newMin = thisBinCenter[i]-0.5* theseDimensions[i]->thisStepSize;
				double newMax = thisBinCenter[i]+0.5* theseDimensions[i]->thisStepSize;
				thisPhaseSpace->SetConstraint( string(theseDim[i]), newMin, newMax, "noUnits_Chi2" );
				//cout << string(theseDim[i]) << "  " << newMin << "  " << newMax << endl;
			}
			vector<string> discNames = thisPhaseSpace->GetDiscreteNames();
			for( unsigned int i=0; i< discNames.size(); ++i )
			{
				double value = thisDataPoint->GetObservable( discNames[i] )->GetValue();
				thisPhaseSpace->SetConstraint( discNames[i], vector<double>(1,value), "noUnits_Chi2" );
				thisPhaseSpace2->SetConstraint( discNames[i], vector<double>(1,value), "noUnits_Chi2" );
			}

			cout << "\tCalculating Integral: " << combinationNum +1 << endl;

			
			//thisPDFIntegrator->ForceTestStatus( false );

			//double Ratio = thisPDFIntegrator->TestIntegral( thisDataPoint, thisPhaseSpace );

			thisPDFIntegrator->ForceTestStatus( true );

			double Total = thisPDFIntegrator->NumericallyIntegrateDataPoint( thisDataPoint, thisPhaseSpace2, thisPDF->GetDoNotIntegrateList() );

			thisPDFIntegrator->clearGSLIntegrationPoints();
			double Integral = thisPDFIntegrator->NumericallyIntegrateDataPoint( thisDataPoint, thisPhaseSpace, thisPDF->GetDoNotIntegrateList() );
			

			/*
			double Total = thisPDF->Integral( thisDataPoint, thisPhaseSpace2 );

			double Integral = thisPDF->Integral( thisDataPoint, thisPhaseSpace );
			*/
			/*
			thisPDFIntegrator->clearGSLIntegrationPoints();
			//Integral = this->CorrectIntegral( Integral, thisDataPoint, thisPhaseSpace, thisPDFIntegrator );
			double newMin = thisBinCenter[0]-0.5* theseDimensions[0]->thisStepSize;
			thisDataPoint->SetObservable( "time", newMin, "noUnits_Chi2" );
			double Integral2 = thisPDFIntegrator->NumericallyIntegrateDataPoint( thisDataPoint, thisPhaseSpace, doNotIntegrate );

			thisPDFIntegrator->clearGSLIntegrationPoints();
			double newMax = thisBinCenter[0]+0.5* theseDimensions[0]->thisStepSize;
			thisDataPoint->SetObservable( "time", newMax, "noUnits_Chi2" );
			double Integral3 = thisPDFIntegrator->NumericallyIntegrateDataPoint( thisDataPoint, thisPhaseSpace, doNotIntegrate );
			*/
			cout << "\tIntegral = " << Integral << endl;//" : " << Integral2 << " : " << Integral3 << endl;
			cout << "\tTotal = " << Total << endl;

			cout << "\tScaling Integral to DataSet." << endl;

			double PDF2DataNorm = this->CorrectYield( thisDataSet, thisDataPoint );

			cout << "\tScale = " << PDF2DataNorm << endl;

			double thisYield = (Integral/Total) * PDF2DataNorm;

			cout << "\tExpected Yield for thisCombination = " << thisYield << endl;

			thisResult += thisYield;
		}

		result_AllPDFs += thisResult;
	}

	cout << "\tReturning Expected Number of Events for this Coordinate" << endl;

	return result_AllPDFs;
}

double MultiDimChi2::CorrectIntegral( double input_Integral, DataPoint* thisPoint, PhaseSpaceBoundary* thisPhaseSpace, RapidFitIntegrator* thisPDFIntegrator )
{
	double output_integral = input_Integral;

	//	Correct for Non Unitary PDF
	output_integral /= 1.;//hisPDFIntegrator->TestIntegral( thisPoint, thisPhaseSpace );

	return output_integral;
}

double MultiDimChi2::CorrectYield( IDataSet* thisSet, DataPoint* thisPoint )
{
	double total_yield = 0.;
	vector<DataPoint*> thesePoints = thisSet->GetDiscreteSubSet( thisPoint );
	if( thisSet->GetWeightsWereUsed() )
	{
		ObservableRef thisRef( thisSet->GetWeightName() );
		for( unsigned int i=0; i< thesePoints.size(); ++i )
		{
			bool isInRange = true;
			/* for( unsigned int j=0; j< theseDimensions.size(); ++j )
			{
				double center = thisPoint->GetObservable( theseDimensions[j]->ObservableName )->GetValue();
				double StepSize = theseDimensions[j]->thisStepSize;
				double new_min = center-0.5*StepSize;
				double new_max = center+0.5*StepSize;

				double thisPointVal = thesePoints[i]->GetObservable( theseDimensions[j]->ObservableName )->GetValue();

				if( thisPointVal > new_max || thisPointVal < new_min )
				{
					isInRange = false;
					break;
				}
			} */

			if( isInRange )
			{
				total_yield += thesePoints[i]->GetEventWeight();//GetObservable( thisRef )->GetValue();
			}
		}
	}
	else
	{
		for( unsigned int i=0; i< thesePoints.size(); ++i )
		{
			bool isInRange = true;
			for( unsigned int j=0; j< theseDimensions.size(); ++j )
			{
				double center = thisPoint->GetObservable( theseDimensions[j]->ObservableName )->GetValue();
				double StepSize = theseDimensions[j]->thisStepSize;
				double new_min = center-0.5*StepSize;
				double new_max = center+0.5*StepSize;

				double thisPointVal = thesePoints[i]->GetObservable( theseDimensions[j]->ObservableName )->GetValue();

				if( thisPointVal > new_max || thisPointVal < new_min )
				{
					isInRange = false;
					break;
				}
			}

			if( isInRange )
			{
				++total_yield;
			}
		}
	}
	return total_yield;
}

void MultiDimChi2::ConstructAllCoordinates()
{
	vector<double> thisBinCenter( theseDimensions.size(), 0. );

	vector<vector<double> > theseBinCenters;

	unsigned int total_points = 1;
	for( unsigned int i=0; i< theseDimensions.size(); ++i )
	{
		total_points*=theseDimensions[i]->theseBins;
	}

	if( allBinCenters != NULL ) delete allBinCenters;
	allBinCenters = new vector<vector<double> > ( total_points, thisBinCenter );

	for( unsigned int i=0; i< nDim; ++i )
	{
		this->AddCoordinates( i );
	}

	/*
	for( unsigned int i=0; i< allBinCenters->size(); ++i )
	{
		for( unsigned int j=0; j< allBinCenters->at(i).size(); ++j )
		{
			cout << (*allBinCenters)[i][j] << "  ";
		}
		cout << endl;
	}
	exit(0);
	*/
}

void MultiDimChi2::AddCoordinates( unsigned int thisDim )
{
	cout << "Adding Coordinated for dimension: " << thisDim << endl;

	unsigned int number_of_set_repeats = 1;
	for( unsigned int i=0; i< thisDim; ++i )
	{
		number_of_set_repeats *= theseDimensions[i]->theseBins;
	}

	cout << "There are: " << (thisDim-0) << " dimensions 'outside' this one." << endl;

	unsigned int number_of_individual_repeats = 1;

	for( unsigned int i=thisDim+1; i< nDim; ++i )
	{
		cout << i << "  " << theseDimensions[i]->theseBins << endl;
		number_of_individual_repeats *= theseDimensions[i]->theseBins;
	}

	cout << "There are: " << ((nDim-1)-thisDim) << " dimensions 'inside' this one." << endl;
	cout << number_of_set_repeats << "  x  " << number_of_individual_repeats << endl;

	unsigned int global_count=0;
	for( unsigned int i=0; i< number_of_set_repeats; ++i )
	{
		for( unsigned int j=0; j< theseDimensions[thisDim]->theseBins; ++j )
		{
			for( unsigned int k=0; k< number_of_individual_repeats; ++k )
			{
				(*allBinCenters)[global_count][thisDim] = theseDimensions[thisDim]->binCenters[j];
				++global_count;
			}
		}
	}
}

void MultiDimChi2::ConstructBoundaries( PhaseSpaceBoundary* totalPhaseSpace, vector<string> wanted_observables )
{
	for( unsigned int i=0; i< allPDFs.size(); ++i )
	{
		vector<string> allDescribedObservables = allPDFs[i]->GetPrototypeDataPoint();
		allBoundaries.push_back( new PhaseSpaceBoundary( allDescribedObservables ) );
		PhaseSpaceBoundary* thisPhaseSpace = allBoundaries.back();
		for( unsigned int j=0; j< allDescribedObservables.size(); ++j )
		{
			//	if( StringProcessing::VectorContains( wanted_observables, allDescribedObservables[j] ) == -1 )
			{
				IConstraint* thisConst = ClassLookUp::CopyConstraint( totalPhaseSpace->GetConstraint( allDescribedObservables[j] ) );
				thisPhaseSpace->SetConstraint( allDescribedObservables[j], thisConst );
				delete thisConst;
			}
		}
		allPDFs[i]->ChangePhaseSpace( thisPhaseSpace );
		allDataSets[i]->SetBoundary( thisPhaseSpace );
	}
}

double MultiDimChi2::CalculateRange( PhaseSpaceBoundary* thisBound )
{
	double range = 1.;
	for( unsigned int i=0; i< theseDimensions.size(); ++i )
	{
		//IConstraint* thisConst = thisBound->GetConstraint( theseDimensions[i]->ObservableName );
		double this_min = x_min[i];//thisConst->GetMinimum();
		double this_max = x_max[i];//thisConst->GetMaximum();
		range *= fabs( this_max - this_min );
	}
	return range;
}

double MultiDimChi2::PDF2DataNormalisation( unsigned int PDFDataNum, const unsigned int combinationIndex, DataPoint* thisDataPoint )
{
	IPDF* thisPDF = allPDFs[ PDFDataNum ];
	IDataSet* thisDataSet = allDataSets[ PDFDataNum ];
	PhaseSpaceBoundary* thisBound = allBoundaries[ PDFDataNum ];
	double normalisation=1.;


	PhaseSpaceBoundary* thisBound3 = new PhaseSpaceBoundary( *thisBound );
	vector<string> wantedParams;
	for( unsigned int i=0; i< theseDimensions.size(); ++i )
	{
		wantedParams.push_back( theseDimensions[i]->ObservableName );
		double value = thisDataPoint->GetObservable( wantedParams.back() )->GetValue();
		thisBound3->SetConstraint( theseDimensions[i]->ObservableName, vector<double>(1, value), "unitless" );
	}
	vector<string> fixed_param = thisBound3->GetDiscreteNames();
	for( unsigned int i=0; i< fixed_param.size(); ++i )
	{
		double value = thisDataPoint->GetObservable( fixed_param[i] )->GetValue();
		thisBound3->SetConstraint( fixed_param[i], vector<double>(1, value), "unitless" );
	}

	vector<string> doNotList = thisPDF->GetDoNotIntegrateList();
	doNotList = StringProcessing::CombineUniques( doNotList, wantedParams );

	RapidFitIntegrator* pdfIntegrator = thisPDF->GetPDFIntegrator();

	pdfIntegrator->ForceTestStatus( false );
	//doNotList = StringProcessing::CombineUniques( doNotList, wanted_params );
	double thisIntegral = 0.;
	//allCombinations[i]->SetPhaseSpaceBoundary( thisBound );
	//cout << thisBound->DiscreteDescription( allCombinations[i] ) << endl;
	double thisRatio =  pdfIntegrator->TestIntegral( thisDataPoint, thisBound3, doNotList );

	cout << "\t\tScaling Based on Numerical/Analytical Ratio: " << 1./thisRatio << endl;
	normalisation /= thisRatio;                    //      Attempt to correct for Numerical != analytical due to any constant factor due to numerical inaccuracy
	//      (some constant close to 1. exactly 1. for numerical PDFs)

	vector<DataPoint*> allCombinations = thisBound->GetDiscreteCombinations();
	DataPoint* thisCombination = allCombinations[combinationIndex];

	double total_yield=0.;

	if( thisDataSet->GetWeightsWereUsed() )
	{
		vector<DataPoint*> thesePoints = thisDataSet->GetDiscreteSubSet( allCombinations[combinationIndex] );

		cout << thesePoints.size() << endl;

		ObservableRef thisWeight( thisDataSet->GetWeightName() );

		for( unsigned int i=0; i< thesePoints.size(); ++i )
		{
			total_yield += thesePoints[i]->GetObservable( thisWeight )->GetValue();
		}
	}
	else
	{
		total_yield = thisDataSet->GetDataNumber( allCombinations[combinationIndex] );
	}	

	double dataNum = total_yield;
	//cout << "\t\tScaling Based on DataNumber in this Combination: " << dataNum / (double) data_binning << endl;
	//normalisation *= dataNum / (double) data_binning;                               //      Normalise to this density of events     (Num of events per bin in flatPDF)


	//double range = this->CalculateRange( thisBound );                                          //      Correct for the range of the dataset    (absolute range of Observable being projected)
	//cout << "\t\tScaling Based on Range of binning: " << range << endl;
	//normalisation *= range;

	cout << total_yield << endl;

	normalisation *= total_yield;

	PhaseSpaceBoundary* thisBound2 = new PhaseSpaceBoundary( *thisBound );

	for( unsigned int AxisNum=0; AxisNum< theseDimensions.size(); ++AxisNum )
	{
		double ThisStep = theseDimensions[AxisNum]->thisStepSize;

		double centralValue = thisDataPoint->GetObservable( theseDimensions[AxisNum]->ObservableName )->GetValue();

		double newMax = centralValue+0.5*ThisStep;
		double newMin = centralValue-0.5*ThisStep;

		cout << "\t\tConstructing New Constraint for: " << theseDimensions[AxisNum]->ObservableName << "\t" << newMin << " : " << newMax << endl;
		thisBound2->SetConstraint( theseDimensions[AxisNum]->ObservableName, newMin, newMax, "noUnit_MultiDimChi2" );
	}

	//thisBound2->Print();

	cout << "\t\tCalculating Bin Integral" << endl;
	//thisDataPoint->Print();
	//thisPDF->ChangePhaseSpace( thisBound2 );
	//IPDF* newPDF = ClassLookUp::CopyPDF( thisPDF );
	//cout << "Copied PDF" << endl;
	//double someNum = newPDF->Normalisation( thisDataPoint, thisBound3 );
	//cout << "Here" << endl;
	//combinationIntegrals[ PDFDataNum ][ combinationIndex ] = someNum;

	delete thisBound2;
	//delete newPDF;

	double someNum = pdfIntegrator->NumericallyIntegrateDataPoint( thisDataPoint, thisBound3, doNotList );

	//cout << "\t\tScaling Based on the Total Integral of the PDF: " << 1./combinationIntegrals[ PDFDataNum ][ combinationIndex ] << endl;
	//normalisation /= combinationIntegrals[ PDFDataNum ][ combinationIndex ];                      //      Total Integral of the PDF       (We're plotting prob of PDF being at this point for a non-unitary PDF)

	//cout << "\t\tScaling Based on the Total Integral of the PDF: " <<  1./ someNum << endl;
	//normalisation /= someNum;

	//cout << "\t\tScaling Based on the Weights in the DataSet: " << 1. / thisDataSet->GetAlpha() << endl;
	//normalisation /= thisDataSet->GetAlpha();                                                   //      Correct for the effect of non-unitary weights used in the fit

	/*
	   if( debug != NULL )
	   {
	   if( debug->DebugThisClass( "MultiDimChi2" ) )
	   {
	   cout << endl;
	   cout << "fabs(ratioOfIntegrals[ combinationIndex ]) : " << fabs(ratioOfIntegrals[ combinationIndex ]) << endl;
	   cout << "thisDataSet->GetDataNumber( allCombinations[combinationIndex] ) : " << thisDataSet->GetDataNumber( allCombinations[combinationIndex] ) << endl;
	   cout << "combinationIndex of allCombinations.size() : " << combinationIndex << " of " << allCombinations.size() << endl;
	   cout << "data_binning : " << data_binning << endl;
	   cout << "fabs( boundary_max-boundary_min ) : " << fabs( boundary_max-boundary_min ) << endl;
	   cout << "combination_integral[ combinationIndex ] : " << combination_integral[ combinationIndex ] << endl;
	   cout << "weight_norm : " << weight_norm << endl;
	   cout << "normalisation : " << normalisation << endl;
	   cout << endl;
	   }
	   }
	   */

	return normalisation;
}

