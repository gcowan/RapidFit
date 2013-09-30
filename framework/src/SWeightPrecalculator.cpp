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
#include "ObservableRef.h"
//	System Headers
#include <math.h>
#include <stdlib.h>
#include <time.h>

SWeightPrecalculator::SWeightPrecalculator( FitResult* InputResult, string WeightName, unsigned int Inputconfig ) :
	inputResult(InputResult), signalPDF(NULL), backgroundPDF(NULL), weightName(WeightName), fractionName(), config( Inputconfig), useAlpha(false)
{
}

void SWeightPrecalculator::SetApplyAlphaCorrection( bool input )
{
	useAlpha = input;
}

void SWeightPrecalculator::ApplyAlphaCorrection( IDataSet* inputDataSet )
{
	double total_sWeight=0.;
	double total_sWeight_sq=0.;
	ObservableRef sWeightName("sWeight");
	ObservableRef sWeightsqName("sWeightSq");
	for( unsigned int i=0; i< (unsigned)inputDataSet->GetDataNumber(); ++i )
	{
		total_sWeight += inputDataSet->GetDataPoint( (int)i )->GetObservable( sWeightName )->GetValue();
		total_sWeight_sq += inputDataSet->GetDataPoint( (int)i )->GetObservable( sWeightsqName )->GetValue();
	}
	double alpha = total_sWeight / total_sWeight_sq;
	cout << total_sWeight << " / " << total_sWeight_sq << endl;
	cout << "Applying Weight: " << alpha << endl;

	for( unsigned int i=0; i< (unsigned)inputDataSet->GetDataNumber(); ++i )
	{
		DataPoint* thisPoint = inputDataSet->GetDataPoint( (int)i );
		double thisWeight = thisPoint->GetObservable( sWeightName )->GetValue();
		thisWeight *= alpha;
		Observable* sWeight_corrected = new Observable( sWeightName, thisWeight, "Weight-Unit" );
		Observable* sWeight_sq_corrected = new Observable( sWeightsqName, thisWeight*thisWeight, "WeightSq-Unit" );
		thisPoint->SetObservable( sWeightName, sWeight_corrected );
		thisPoint->SetObservable( sWeightsqName, sWeight_sq_corrected );
		delete sWeight_corrected;
		delete sWeight_sq_corrected;
	}
}

void SWeightPrecalculator::ConfigurePDFs( IPDF* InputPDF )
{
	if( InputPDF->GetName() == "NormalisedSumPDF" )
	{
		NormalisedSumPDF* inputNorm = (NormalisedSumPDF*)InputPDF;
		if( config == 1 )
		{
			if( signalPDF != NULL ) delete signalPDF;
			signalPDF = ClassLookUp::CopyPDF( inputNorm->GetFirstPDF() );
			if( backgroundPDF != NULL ) delete backgroundPDF;
			backgroundPDF = ClassLookUp::CopyPDF( inputNorm->GetSecondPDF() );
		}
		else if ( config == 2 )
		{
			if( signalPDF != NULL ) delete signalPDF;
			signalPDF = ClassLookUp::CopyPDF( inputNorm->GetSecondPDF() );
			if( backgroundPDF != NULL ) delete backgroundPDF;
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

	ParameterSet* inputSet = inputResult->GetResultParameterSet()->GetDummyParameterSet();
	signalPDF->UpdatePhysicsParameters( inputSet );
	backgroundPDF->UpdatePhysicsParameters( inputSet );
}

//Destructor
SWeightPrecalculator::~SWeightPrecalculator()
{
	if( signalPDF != NULL ) delete signalPDF;
	if( backgroundPDF != NULL ) delete backgroundPDF;
	// Can't copy FitResult, so no need to attempt to destroy it here!!!
}

//Calculate the sWeights
IDataSet * SWeightPrecalculator::ProcessDataSet( IDataSet * InputData, IPDF* InputPDF )
{
	this->ConfigurePDFs( InputPDF );

	cout << endl << "Calculating sWeights" << endl;

	//Retrieve the correct fraction
	double signalFraction = inputResult->GetResultParameterSet()->GetResultParameter( fractionName )->GetValue();
	double numberSignalEvents = /*(int)floor(*/ (double)(InputData->GetDataNumber()) * signalFraction;// );
	double numberBackgroundEvents = /*(int)ceil(*/ (double)(InputData->GetDataNumber()) * ( 1. - signalFraction );// );

	//Debug
	cout << "Signal events: " << numberSignalEvents << endl;
	cout << "Background events: " << numberBackgroundEvents << endl;

	cout << "SignalFraction: " << fractionName << "\t" << signalFraction << endl;
	cout << "Number of Events: " << InputData->GetDataNumber() << endl;

	double denom_2=0.; vector<double> numer_2;

	//Calculate the matrix elements to use in the sWeight
	vector<double> signalValues, backgroundValues;
	pair< double, double > matrixElements = CalculateMatrixElements( numberSignalEvents, numberBackgroundEvents, InputData, signalValues, backgroundValues, denom_2, numer_2 );

	//Now loop through input data, calculate sWeights and make a new DataSet
	vector<string> allObservables = InputData->GetBoundary()->GetAllNames();
	vector<string>::iterator observableIterator;
	allObservables.push_back(weightName);

	vector<DataPoint*> allPoints;
	vector<double> allValues, allValues2;

	string weightName2 = weightName+"Sq";

	double sum = 0.;
	double sum2 = 0.;
    double sum_weight_sq = 0.;
	for ( int eventIndex = 0; eventIndex < InputData->GetDataNumber(); ++eventIndex )
	{
		//Calculate the sWeight
		double numerator = double( ( matrixElements.first * signalValues[unsigned(eventIndex)] )
				+ ( matrixElements.second * backgroundValues[unsigned(eventIndex)] ) );

		double numerator2 =
		(numer_2[0]*signalValues[unsigned(eventIndex)])*(numer_2[0]*signalValues[unsigned(eventIndex)]) +
		(numer_2[1]*backgroundValues[unsigned(eventIndex)])*(numer_2[1]*backgroundValues[unsigned(eventIndex)])
		+ 2.*(numer_2[0]*signalValues[unsigned(eventIndex)])
		    *(numer_2[1]*backgroundValues[unsigned(eventIndex)]);
		numerator2/=denom_2;

		double denominator = double( ( double(numberSignalEvents) * signalValues[unsigned(eventIndex)] )
				+ ( double(numberBackgroundEvents) * backgroundValues[unsigned(eventIndex)] ) );

		double denominator2 =
		(double(numberSignalEvents)*signalValues[unsigned(eventIndex)])
		*(double(numberSignalEvents)*signalValues[unsigned(eventIndex)]) +
		(double(numberBackgroundEvents)*backgroundValues[unsigned(eventIndex)])
		*(double(numberBackgroundEvents)*backgroundValues[unsigned(eventIndex)])
		+ 2.*(double(numberSignalEvents)*signalValues[unsigned(eventIndex)])
		    *(double(numberBackgroundEvents)*backgroundValues[unsigned(eventIndex)]);

		//Make the new data point
		DataPoint * currentEvent = InputData->GetDataPoint(eventIndex);
		DataPoint * newEvent = new DataPoint( *currentEvent );

		//Add the sWeight as an Observable
		newEvent->AddObservable( weightName, numerator / denominator, "Unitless" );
		newEvent->AddObservable( weightName2, numerator2 / denominator2, "Unitless" );

                sum += numerator / denominator;
                sum2 += numerator2 / denominator2;

		allValues.push_back( numerator / denominator );
		allValues2.push_back( numerator2 / denominator2 );
		allPoints.push_back( newEvent );

        sum_weight_sq += (numerator / denominator)*(numerator / denominator);
	}

	cout << sum << " " << sum2 << endl;

    cout << "FOM " << sum*sum/sum_weight_sq << endl;

        PhaseSpaceBoundary* dataSetBoundary = new PhaseSpaceBoundary( *(InputData->GetBoundary()) );
	double min=0., max=0., min2=0., max2=0.;
	min = StatisticsFunctions::Minimum( allValues );
	max = StatisticsFunctions::Maximum( allValues );
	min2 = StatisticsFunctions::Minimum( allValues2 );
	max2 = StatisticsFunctions::Maximum( allValues2 );
	ObservableContinuousConstraint* weightConstraint = new ObservableContinuousConstraint( weightName, min, max, "Unitless" );
	ObservableContinuousConstraint* weightConstraint2 = new ObservableContinuousConstraint( weightName2, min2, max2, "Unitless" );
	dataSetBoundary->AddConstraint( weightName, weightConstraint );
	dataSetBoundary->AddConstraint( weightName2, weightConstraint2 );
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

	if( useAlpha ) ApplyAlphaCorrection( (IDataSet*)newDataSet );

	return (IDataSet*) newDataSet;
}

//A separate method for calculating the martix elements for the weight
pair< double, double > SWeightPrecalculator::CalculateMatrixElements( double NumberSignal, double NumberBackground, IDataSet * InputData,
		vector<double> & SignalValues, vector<double> & BackgroundValues, double& denom_2, vector<double>& numer_2 )
{
	//The three unique components of the matrix (one is used twice)
	double signalSignal = 0.0;
	double signalBackground = 0.0;
	double backgroundBackground = 0.0;

	//Store the PDF evaluations to save time
	vector<double> saveSignalValues, saveBackgroundValues;

	//Make the PDF integrators
	RapidFitIntegrator * signalIntegrator = new RapidFitIntegrator( signalPDF );
	signalIntegrator->ForceTestStatus( true );
	RapidFitIntegrator * backgroundIntegrator = new RapidFitIntegrator( backgroundPDF );
	backgroundIntegrator->ForceTestStatus( true );

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

	double den_2 = (signalSignal*backgroundBackground)*(signalSignal*backgroundBackground)
			+ (signalBackground*signalBackground)*(signalBackground*signalBackground)
			- 2.*(signalSignal*backgroundBackground)*(signalBackground*signalBackground);
	denom_2 = den_2;

	vector<double> num_2; num_2.push_back( backgroundBackground ); num_2.push_back( -signalBackground );
	numer_2 = num_2;

	//Return results
	SignalValues = saveSignalValues;
	BackgroundValues = saveBackgroundValues;
	delete signalIntegrator;
	delete backgroundIntegrator;
	return make_pair( returnSignalSignal, returnSignalBackground );
}

