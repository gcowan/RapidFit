

#include "THn.h"
#include "MultiDimChi2.h"
#include "PhaseSpaceBoundary.h"
#include "IDataSet.h"
#include <vector>
#include "StringProcessing.h"
#include "IPDF.h"
#include <string>
#include "ClassLookUp.h"
#include "PDFWithData.h"

using namespace::std;

MultiDimChi2::MultiDimChi2( vector<PDFWithData*> allObjects, PhaseSpaceBoundary* thisBound, vector<string> wantedObservables )
{
	//	Populate PDF, DataSets and PhaseSpaceBoundaries to be used in this analysis
	cout << "MultiDimChi2 Hello!" << endl << endl;
	cout << "initializing objects" << endl;
	this->populateAllObjects( allObjects );
	this->ConstructBoundaries( thisBound );

	cout << "Deciding binning and constructing THnD" << endl;
	//	Construct the Histogram object which contains all of the Data of Interest
	this->ConstructInternalHisto( wantedObservables, thisBound );

	cout << "Calculating Bin Centers" << endl;
	//      Construct the Central coordinates of every bin in each axis, this will be the basis for generating the corrdinates to run the Chi2 test over
	this->ConstructBinCenters();

	cout << "Calculating the Numerical vs Analytical Ratios" << endl;
	//	Construct the Ratios of the Numerical vs analytical Integrals (condensing the differenced between the PDF numerator and denominator to be a single number but what else can we do...)
	this->ConstructIntegralsRatios();

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
	}
}

void MultiDimChi2::ConstructIntegralsRatios()
{
	ratioOfIntegrals.clear();
	combinationIntegrals.clear();
	weightNorms.clear();

	for( unsigned int PDFNum=0; PDFNum< allPDFs.size(); ++PDFNum )
	{
		IPDF* thisPDF = allPDFs[PDFNum];
		RapidFitIntegrator* pdfIntegrator = thisPDF->GetPDFIntegrator();
		IDataSet* thisData = allDataSets[PDFNum];
		vector<double> thisratioOfIntegrals;
		vector<double> thisCombination_integral;
		PhaseSpaceBoundary* thisBound = allBoundaries[PDFNum];

		vector<DataPoint*> allCombinations = allBoundaries[PDFNum]->GetDiscreteCombinations();

		for( unsigned int i=0; i< allCombinations.size(); ++i )
		{
			/*
			   if( debug != NULL )
			   {
			   if( debug->DebugThisClass( "MultiDimChi2" ) )
			   {
			   cout << "MultiDimChi2:: Calculating Test integrals:\t" << i << endl;
			   }
			   }*/
			pdfIntegrator->ForceTestStatus( false );
			double thisIntegral = 0.;
			try
			{
				thisIntegral = pdfIntegrator->NumericallyIntegrateDataPoint( allCombinations[i], thisBound, thisPDF->GetDoNotIntegrateList() );
			}
			catch(...)
			{
				thisIntegral = 1.;
				cout << endl << "CANNOT PROPERLY NORMALISE WHOLE PDF, THIS WILL LEAD TO NORMALISATION ISSUES OVER THE WHOLE PDF" << endl << endl;
			}
			/*
			   if( debug != NULL )
			   {
			   if( debug->DebugThisClass( "MultiDimChi2" ) )
			   {
			   cout << "MultiDimChi2:: Finished Inetgral" << endl;
			   }
			   }*/

			thisCombination_integral.push_back( thisIntegral );

			if( thisPDF->GetNumericalNormalisation() == true ) thisratioOfIntegrals.push_back( 1. );
			else
			{
				thisratioOfIntegrals.push_back( pdfIntegrator->GetRatioOfIntegrals() );
			}
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
		combinationIntegrals.push_back( thisCombination_integral );

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
		   }*/
		weightNorms.push_back( thisWeight_norm );
	}
}

void MultiDimChi2::ConstructInternalHisto( vector<string> wantedObservables, PhaseSpaceBoundary* thisBound )
{
	nDim=0;
	x_min.clear(); x_max.clear(); goodObservables.clear(); x_bins.clear();

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
			x_binning = 10;
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

	if( internalHisto != NULL ) delete internalHisto;

	internalHisto = new THnD( "internal_Chi2nDim", "internal_Chi2nDim", (int)nDim, &(x_bins[0]), &(x_min[0]), &(x_max[0]) );

	vector<double> thisValues( goodObservables.size(), 0. );
	double thisWeight = 0.;
	DataPoint* thisPoint=NULL;

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
		thisChi2 = expected_events[i] - observed_events[i] + observed_events[i]*log( observed_events[i] / expected_events[i] );
		Chi2Value += thisChi2;
	}
	Chi2Value*=2.;
	return Chi2Value;
}

double MultiDimChi2::CalculateTotalExpected( vector<double> thisBinCenter )
{

	double result_AllPDFs = 0.;

	cout << "Calculating Expected Number of Events for this Bin" << endl;

	IPDF* thisPDF=NULL;
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
		RapidFitIntegrator* thisPDFIntegrator = thisPDF->GetPDFIntegrator();

		vector<string> doNotIntegrate = thisPDF->GetDoNotIntegrateList();

		doNotIntegrate = StringProcessing::CombineUniques( doNotIntegrate, theseDimName );

		PhaseSpaceBoundary* thisParamSet = allBoundaries[PDFNum];

		vector<DataPoint*> theseDataPoints = thisParamSet->GetDiscreteCombinations();

		double thisResult=0.;

		cout << "I have " << theseDataPoints.size()+1 << " integrals to perform for PhaseSpace: " << PDFNum+1 << endl;

		for( unsigned int combinationNum = 0; combinationNum< theseDataPoints.size(); ++combinationNum )
		{
			thisDataPoint = theseDataPoints[combinationNum];

			for( unsigned int i=0; i< theseDim.size(); ++i )
			{
				Observable* newObs = new Observable( string(theseDim[i]), thisBinCenter[i], "noUnits_Chi2" );
				thisDataPoint->SetObservable( theseDim[i], newObs );
				delete newObs;
			}

			cout << "Calculating Integral: " << combinationNum +1 << endl;

			double Integral = thisPDFIntegrator->NumericallyIntegrateDataPoint( thisDataPoint, thisParamSet, doNotIntegrate ); 

			cout << "Integral = " << Integral << endl;

			cout << "Scaling Integral to DataSet." << endl;

			double PDF2DataNorm = this->PDF2DataNormalisation( PDFNum, combinationNum );

			cout << "Scale = " << PDF2DataNorm << endl;

			double thisYield=Integral*PDF2DataNorm;

			cout << "Expected Yield for thisCombination = " << thisYield << endl;

			thisResult += thisYield;
		}

		result_AllPDFs += thisResult;
	}

	cout << "Returning Expected Number of Events for this Coordinate" << endl;

	return result_AllPDFs;
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
}

void MultiDimChi2::AddCoordinates( unsigned int thisDim )
{
	unsigned int number_of_set_repeats = 1;
	for( unsigned int i=0; i< thisDim; ++i )
	{
		number_of_set_repeats *= theseDimensions[i]->theseBins;
	}

	unsigned int number_of_individual_repeats = 1;

	for( unsigned int i=thisDim; i< nDim; ++i )
	{
		number_of_individual_repeats *= theseDimensions[i]->theseBins;
	}

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

void MultiDimChi2::ConstructBoundaries( PhaseSpaceBoundary* totalPhaseSpace )
{
	for( unsigned int i=0; i< allPDFs.size(); ++i )
	{
		vector<string> allDescribedObservables = allPDFs[i]->GetPrototypeDataPoint();
		allBoundaries.push_back( new PhaseSpaceBoundary( allDescribedObservables ) );
		PhaseSpaceBoundary* thisPhaseSpace = allBoundaries.back();
		for( unsigned int j=0; j< allDescribedObservables.size(); ++j )
		{
			IConstraint* thisConst = totalPhaseSpace->GetConstraint( allDescribedObservables[j] );
			thisPhaseSpace->SetConstraint( allDescribedObservables[j], thisConst );
			delete thisConst;
		}
		allPDFs[i]->ChangePhaseSpace( thisPhaseSpace );
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

double MultiDimChi2::PDF2DataNormalisation( unsigned int PDFDataNum, const unsigned int combinationIndex )
{
	IDataSet* thisDataSet = allDataSets[ PDFDataNum ];
	PhaseSpaceBoundary* thisBound = allBoundaries[ PDFDataNum ];
	double normalisation=1.;

	normalisation *= fabs(ratioOfIntegrals[ PDFDataNum ][ combinationIndex ]);                    //      Attempt to correct for Numerical != analytical due to any constant factor due to numerical inaccuracy
	//      (some constant close to 1. exactly 1. for numerical PDFs)

	vector<DataPoint*> allCombinations = thisBound->GetDiscreteCombinations();

	double dataNum = thisDataSet->GetDataNumber( allCombinations[combinationIndex] );
	normalisation *= dataNum / (double) data_binning;                               //      Normalise to this density of events     (Num of events per bin in flatPDF)

	double range = this->CalculateRange( thisBound );                                          //      Correct for the range of the dataset    (absolute range of Observable being projected)
	normalisation *= range;

	normalisation /= combinationIntegrals[ PDFDataNum ][ combinationIndex ];                      //      Total Integral of the PDF       (We're plotting prob of PDF being at this point for a non-unitary PDF)

	normalisation *= weightNorms[ PDFDataNum ];                                                   //      Correct for the effect of non-unitary weights used in the fit

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

