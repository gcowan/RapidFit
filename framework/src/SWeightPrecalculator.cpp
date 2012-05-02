/**
  @class SWeightPrecalculator

  A Precalculator for producing the sWeights for a data set for backgorund subtraction

  @author Benjamin Wynne bwynne@cern.ch
  @date 2009-12-14
 */

//	RapidFit Headers
#include "SWeightPrecalculator.h"
#include "NormalisedSumPDF.h"
#include "FitAssembler.h"
#include "MemoryDataSet.h"
#include "RapidFitIntegrator.h"
#include "ClassLookUp.h"
#include "StatisticsFunctions.h"
#include "ObservableContinuousConstraint.h"
//	System Headers
#include <math.h>
#include <stdlib.h>
#include <time.h>

SWeightPrecalculator::SWeightPrecalculator( FitResult* InputResult, string WeightName, unsigned int config ) : inputResult(InputResult), signalPDF(NULL), backgroundPDF(NULL), weightName(WeightName), fractionName()
{
	if( inputResult->GetPhysicsBottle()->GetResultPDF( 0 )->GetName() == "NormalisedSum" )
	{
		NormalisedSumPDF* inputNorm = (NormalisedSumPDF*)inputResult->GetPhysicsBottle()->GetResultPDF( 0 );
		if( config == 1 )
		{
			signalPDF = ClassLookUp::CopyPDF( inputNorm->GetFirstPDF() );
			backgroundPDF = ClassLookUp::CopyPDF( inputNorm->GetSecondPDF() );
		}
		else if ( config == 2 )
		{
			signalPDF = ClassLookUp::CopyPDF( inputNorm->GetSecondPDF() );
			backgroundPDF = ClassLookUp::CopyPDF( inputNorm->GetFirstPDF() );
		}
		else
		{
			cerr << "Unknown Option, Please Select configuration 1 or 2" << endl;
			exit(-9823);
		}
		fractionName = inputNorm->GetFractionName();
	}
	else
	{
		cerr << "Provided PDF, NOT a NormalisedSumPDF, cannot Proceed with S-Weighting" << endl;
		exit(-3);
	}

}

//Destructor
SWeightPrecalculator::~SWeightPrecalculator()
{
	if( signalPDF != NULL ) delete signalPDF;
	if( backgroundPDF != NULL ) delete backgroundPDF;
	// Can't copy FitResult, so no need to attempt to destroy it here!!!
}

//Calculate the sWeights
IDataSet * SWeightPrecalculator::ProcessDataSet( IDataSet * InputData )
{
	cout << endl << "Calculating sWeights" << endl;

	//Retrieve the correct fraction
	double signalFraction = inputResult->GetResultParameterSet()->GetResultParameter( fractionName )->GetValue();
	long numberSignalEvents = (long)floor( InputData->GetDataNumber() * signalFraction );
	long numberBackgroundEvents = (long)ceil( InputData->GetDataNumber() * ( 1 - signalFraction ) );

	//Debug
	cout << "Signal events: " << numberSignalEvents << endl;
	cout << "Background events: " << numberBackgroundEvents << endl;

	//Calculate the matrix elements to use in the sWeight
	vector<double> signalValues, backgroundValues;
	pair< double, double > matrixElements = CalculateMatrixElements( numberSignalEvents, numberBackgroundEvents, InputData, signalValues, backgroundValues );

	//Now loop through input data, calculate sWeights and make a new DataSet
	vector<string> allObservables = InputData->GetBoundary()->GetAllNames();
	vector<string>::iterator observableIterator;
	allObservables.push_back(weightName);

	vector<DataPoint*> allPoints;
	vector<double> allValues;

	for ( int eventIndex = 0; eventIndex < InputData->GetDataNumber(); ++eventIndex )
	{
		//Calculate the sWeight
		double numerator = double( ( matrixElements.first * signalValues[unsigned(eventIndex)] ) + ( matrixElements.second * backgroundValues[unsigned(eventIndex)] ) );
		double denominator = double( ( double(numberSignalEvents) * signalValues[unsigned(eventIndex)] ) + ( double(numberBackgroundEvents) * backgroundValues[unsigned(eventIndex)] ) );

		//Make the new data point
		DataPoint * currentEvent = InputData->GetDataPoint(eventIndex);
		DataPoint * newEvent = new DataPoint( *currentEvent );

		//Add the sWeight as an Observable
		newEvent->AddObservable( weightName, numerator / denominator, 0.0, "Unitless" );

		allValues.push_back( numerator / denominator );
		allPoints.push_back( newEvent );
	}

        PhaseSpaceBoundary* dataSetBoundary = new PhaseSpaceBoundary( *(InputData->GetBoundary()) );
	double min=0, max=0;
	min = StatisticsFunctions::Minimum( allValues );
	max = StatisticsFunctions::Maximum( allValues );
	ObservableContinuousConstraint* weightConstraint = new ObservableContinuousConstraint( weightName, min, max, "Unitless" );
	dataSetBoundary->AddConstraint( weightName, weightConstraint );
        MemoryDataSet * newDataSet = new MemoryDataSet( dataSetBoundary );

	for( unsigned int i=0; i< allPoints.size(); ++i )
	{
		newDataSet->AddDataPoint( allPoints[i] );
	}

	//Output time information
	time_t timeNow;
	time(&timeNow);
	cout << "Finished calculating sWeights: " << ctime( &timeNow ) << endl;
	//exit(0);

	return (IDataSet*) newDataSet;
}

//A separate method for calculating the martix elements for the weight
pair< double, double > SWeightPrecalculator::CalculateMatrixElements( long NumberSignal, long NumberBackground, IDataSet * InputData, vector<double> & SignalValues, vector<double> & BackgroundValues )
{
	//The three unique components of the matrix (one is used twice)
	double signalSignal = 0.0;
	double signalBackground = 0.0;
	double backgroundBackground = 0.0;

	//Store the PDF evaluations to save time
	vector<double> saveSignalValues, saveBackgroundValues;

	//Make the PDF integrators
	RapidFitIntegrator * signalIntegrator = (RapidFitIntegrator*) ( signalPDF->RequestIntegrator() );
	RapidFitIntegrator * backgroundIntegrator = (RapidFitIntegrator*) ( backgroundPDF->RequestIntegrator() );

	//The matrix is a sum over all events
	for ( int eventIndex = 0; eventIndex < InputData->GetDataNumber(); ++eventIndex )
	{
		//Evaluate signal and background PDFs for the point
		DataPoint * currentEvent = InputData->GetDataPoint(eventIndex);
		double signalValue = signalPDF->Evaluate(currentEvent);
		double backgroundValue = backgroundPDF->Evaluate(currentEvent);

		//Normalise function values
		signalValue /= signalIntegrator->Integral( currentEvent, InputData->GetBoundary() );
		backgroundValue /= backgroundIntegrator->Integral( currentEvent, InputData->GetBoundary() );

		//Store fucntion values
		saveSignalValues.push_back(signalValue);
		saveBackgroundValues.push_back(backgroundValue);

		//Do the matrix calculations
		double sqrtDenominator = ( double(NumberSignal) * signalValue ) + ( double(NumberBackground) * backgroundValue );
		signalSignal += ( signalValue * signalValue ) / ( sqrtDenominator * sqrtDenominator );
		signalBackground += ( signalValue * backgroundValue ) / ( sqrtDenominator * sqrtDenominator );
		backgroundBackground += ( backgroundValue * backgroundValue ) / ( sqrtDenominator * sqrtDenominator );
	}

	//Calculate the required components of the inverse
	double discriminant = ( signalSignal * backgroundBackground ) - ( signalBackground * signalBackground );
	double returnSignalSignal = backgroundBackground / discriminant;
	double returnSignalBackground = -signalBackground / discriminant;

	//Return results
	SignalValues = saveSignalValues;
	BackgroundValues = saveBackgroundValues;
	return make_pair( returnSignalSignal, returnSignalBackground );
}

