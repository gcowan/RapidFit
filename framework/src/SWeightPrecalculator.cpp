/**
  @class SWeightPrecalculator

  A Precalculator for producing the sWeights for a data set for backgorund subtraction

  @author Benjamin Wynne bwynne@cern.ch
  @date 2009-12-14
 */

#include "SWeightPrecalculator.h"
#include "NormalisedSumPDF.h"
#include "FitAssembler.h"
#include "MemoryDataSet.h"
#include "RapidFitIntegrator.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>

//Default constructor
SWeightPrecalculator::SWeightPrecalculator()
{
}

//Constructor with correct arguments
SWeightPrecalculator::SWeightPrecalculator( IPDF * InputSignal, IPDF * InputBackground, ParameterSet * InputParameters, string WeightName ) : signalPDF(InputSignal),
	backgroundPDF(InputBackground), fitParameters(InputParameters), weightName(WeightName)
{
}

//Destructor
SWeightPrecalculator::~SWeightPrecalculator()
{
}

//Calculate the sWeights
IDataSet * SWeightPrecalculator::ProcessDataSet( IDataSet * InputData )
{
	cout << endl << "Calculating sWeights" << endl;

	//First fit to find number of signal events v number background
	//Make a sum PDF, with a new physics parameter for the ratio
	const string fractionName = "sWeightSignalFraction";
	NormalisedSumPDF * signalAndBackground = new NormalisedSumPDF( signalPDF, backgroundPDF, InputData->GetBoundary(), fractionName );

	//Make a new parameter set with this parameter in
	vector<string> allParameters = fitParameters->GetAllNames();
	vector<string>::iterator parameterIterator;
	allParameters.push_back(fractionName);
	ParameterSet * fractionFitParameters = new ParameterSet(allParameters);
	for ( parameterIterator = allParameters.begin(); parameterIterator != allParameters.end(); parameterIterator++ )
	{
		if  ( *parameterIterator == fractionName )
		{
			//Add the new parameter
			fractionFitParameters->SetPhysicsParameter( fractionName, 0.5, 0.0, 1.0, "Free", "Unitless" );
		}
		else
		{
			//Copy the parameter set
			PhysicsParameter temporaryParameter = *( fitParameters->GetPhysicsParameter( *parameterIterator ) );
			fractionFitParameters->SetPhysicsParameter( *parameterIterator, &temporaryParameter );
		}
	}

	//Fit
	vector< IPDF* > fitPDF;
	fitPDF.push_back(signalAndBackground);
	vector< IDataSet* > fitData;
	fitData.push_back(InputData);
	MinimiserConfiguration * minimiser = new MinimiserConfiguration("Minuit2");
	FitFunctionConfiguration * function = new FitFunctionConfiguration("NegativeLogLikelihood");
	FitResult * findFractionResult = FitAssembler::DoFit( minimiser, function, fractionFitParameters, fitPDF, fitData, vector< ConstraintFunction* >() );

	//Retrieve the correct fraction
	double signalFraction = findFractionResult->GetResultParameterSet()->GetResultParameter(fractionName)->GetValue();
	long numberSignalEvents = (long)floor( InputData->GetDataNumber() * signalFraction );
	long numberBackgroundEvents = (long)ceil( InputData->GetDataNumber() * ( 1 - signalFraction ) );

	//Reset the PDF parameters (just in case)
	signalPDF->SetPhysicsParameters(fitParameters);
	backgroundPDF->SetPhysicsParameters(fitParameters);

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
	MemoryDataSet * newDataSet = new MemoryDataSet( InputData->GetBoundary() );
	for ( int eventIndex = 0; eventIndex < InputData->GetDataNumber(); eventIndex++ )
	{
		//Calculate the sWeight
		double numerator = double( ( matrixElements.first * signalValues[eventIndex] ) + ( matrixElements.second * backgroundValues[eventIndex] ) );
		double denominator = double( ( double(numberSignalEvents) * signalValues[eventIndex] ) + ( double(numberBackgroundEvents) * backgroundValues[eventIndex] ) );
		
		//Make the new data point
		DataPoint * currentEvent = InputData->GetDataPoint(eventIndex);
		DataPoint * newEvent = new DataPoint(allObservables);
		for ( observableIterator = allObservables.begin(); observableIterator != allObservables.end(); observableIterator++ )
		{
			if ( *observableIterator == weightName )
			{
				//Add the sWeight as an Observable
				newEvent->SetObservable( weightName, numerator / denominator, 0.0, "Unitless" );
				newDataSet->AddDataPoint(newEvent);
			}
			else
			{
				//Copy data point contents
				Observable * copyObservable = currentEvent->GetObservable( *observableIterator );
				newEvent->SetObservable( *observableIterator, copyObservable );
			}
		}
	}

	//Output time information
	time_t timeNow;
        time(&timeNow);
	cout << "Finished calculating sWeights: " << ctime( &timeNow ) << endl;
	//exit(0);

	return newDataSet;
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
	RapidFitIntegrator * signalIntegrator = new RapidFitIntegrator(signalPDF);
	RapidFitIntegrator * backgroundIntegrator = new RapidFitIntegrator(backgroundPDF);

	//The matrix is a sum over all events
	for ( int eventIndex = 0; eventIndex < InputData->GetDataNumber(); eventIndex++ )
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
