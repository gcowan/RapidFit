/**
  @class EfficiencyWeightPreCalculator

  A Precalculator for producing the weights for efficiency weighting a dataset

  @author Greig Cowan greig.cowan@cern.ch
  @date 2013-10-18
 */

//	RapidFit Headers
#include "EfficiencyWeightPreCalculator.h"
#include "NormalisedSumPDF.h"
#include "FitAssembler.h"
#include "MemoryDataSet.h"
#include "RapidFitIntegrator.h"
#include "ClassLookUp.h"
#include "StatisticsFunctions.h"
#include "ObservableContinuousConstraint.h"
#include "ObservableRef.h"
//	System Headers
#include <math.h>
#include <stdlib.h>
#include <time.h>

EfficiencyWeightPreCalculator::EfficiencyWeightPreCalculator( FitResult* InputResult, string WeightName, unsigned int Inputconfig ) :
	inputResult(InputResult), signalPDF(NULL), weightName(WeightName), config( Inputconfig)
{
}

void EfficiencyWeightPreCalculator::SetApplyAlphaCorrection( bool b ){}

void EfficiencyWeightPreCalculator::ConfigurePDFs( IPDF* InputPDF )
{
	ParameterSet* inputSet = inputResult->GetResultParameterSet()->GetDummyParameterSet();
	signalPDF->UpdatePhysicsParameters( inputSet );
}

//Destructor
EfficiencyWeightPreCalculator::~EfficiencyWeightPreCalculator()
{
	if( signalPDF != NULL ) delete signalPDF;
	// Can't copy FitResult, so no need to attempt to destroy it here!!!
}

//Calculate the sWeights
IDataSet * EfficiencyWeightPreCalculator::ProcessDataSet( IDataSet * InputData, IPDF* InputPDF )
{
	//this->ConfigurePDFs( InputPDF );

	cout << endl << "Efficiency weighting dataset" << endl;

	//Now loop through input data, calculate sWeights and make a new DataSet
	vector<string> allObservables = InputData->GetBoundary()->GetAllNames();
	vector<string>::iterator observableIterator;
	allObservables.push_back(weightName);

	vector<DataPoint*> allPoints;
	vector<double> allValues;
	
	RapidFitIntegrator * testIntegrator = new RapidFitIntegrator( InputPDF, true, true );
        vector<string> pdfComponents  = InputPDF->PDFComponents();
	vector<string> doNotIntegrate = InputPDF->GetDoNotIntegrateList();
        ComponentRef * thisRef = new ComponentRef( "0", "dummyObservable" );
	double pdf_norm = testIntegrator->NumericallyIntegratePhaseSpace( InputData->GetBoundary(), doNotIntegrate, thisRef );
        	
	cout << "Number of events in non-weighted dataset: " << InputData->GetDataNumber() << endl;
    	double pdf_value(0.);
	double sumOfWeights(0.);
	for ( int eventIndex = 0; eventIndex < InputData->GetDataNumber(); ++eventIndex )
	{
		//Make the new data point
		DataPoint * currentEvent = InputData->GetDataPoint(eventIndex);
		DataPoint * newEvent = new DataPoint( *currentEvent );

        	pdf_value = InputPDF->Evaluate(currentEvent);
		sumOfWeights += pdf_value/pdf_norm;
		//Add the sWeight as an Observable
		newEvent->AddObservable( weightName, pdf_value/pdf_norm, "Unitless" );

		allValues.push_back( pdf_value );
		allPoints.push_back( newEvent );
	}

    PhaseSpaceBoundary* dataSetBoundary = new PhaseSpaceBoundary( *(InputData->GetBoundary()) );
	double min=0., max=0.;
	min = StatisticsFunctions::Minimum( allValues );
	max = StatisticsFunctions::Maximum( allValues );
	min = 0.; max = 1.; //GAC fix
	ObservableContinuousConstraint* weightConstraint = new ObservableContinuousConstraint( weightName, min, max, "Unitless" );
	dataSetBoundary->AddConstraint( weightName, weightConstraint );
    MemoryDataSet * newDataSet = new MemoryDataSet( dataSetBoundary );

	cout << "min, max, allPoints.size(), sum of weights: " << min << " " << max << " " << allPoints.size() << sumOfWeights << endl;
	for( unsigned int i=0; i< allPoints.size(); ++i )
	{
		newDataSet->AddDataPoint( allPoints[i] );
	}

	//Output time information
	time_t timeNow;
	time(&timeNow);
	cout << "Number of events in weighted dataset: " << newDataSet->GetDataNumber() << endl;
	cout << "Finished weighting: " << ctime( &timeNow ) << endl;
	//exit(0);

	return (IDataSet*) newDataSet;
}

