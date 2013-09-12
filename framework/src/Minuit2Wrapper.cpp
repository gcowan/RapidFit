/**
  @class Minuit2Wrapper

  A wrapper for Minuit2, implementing IMinimiser

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	ROOT Headers
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MinosError.h"
#include "TMatrixDSym.h"
//	RapidFit Headers
#include "Minuit2Wrapper.h"
#include "ResultParameterSet.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <limits>
#include <ctime>

//const double MAXIMUM_MINIMISATION_STEPS = 100000.0;//800.0;
//const double FINAL_GRADIENT_TOLERANCE = 0.01;//;0.001;
const double STEP_SIZE = 0.01;
//const int MINUIT_QUALITY = 2;

//Default constructor
Minuit2Wrapper::Minuit2Wrapper() :
	function(NULL), RapidFunction(NULL), fitResult(NULL), contours(), maxSteps(), bestTolerance(), Options(), Quality(), debug( new DebugClass(false) ), nSigma(1), minimum(NULL)
{
}

//Destructor
Minuit2Wrapper::~Minuit2Wrapper()
{
	if( debug != NULL ) delete debug;
	if( minimum != NULL ) delete minimum;
}

void Minuit2Wrapper::SetSteps( int newSteps )
{
	maxSteps = newSteps;
}

void Minuit2Wrapper::SetTolerance( double newTolerance )
{
	bestTolerance = newTolerance;
}

void Minuit2Wrapper::SetOptions( vector<string> newOptions )
{
	Options = newOptions;
}

void Minuit2Wrapper::SetQuality( int newQuality )
{
	Quality = newQuality;
}

void Minuit2Wrapper::SetupFit( IFitFunction* NewFunction )
{
	//Make a wrapper for the function
	function = new Minuit2Function( NewFunction, nSigma );
	RapidFunction = NewFunction;
}

IFitFunction* Minuit2Wrapper::GetFitFunction()
{
	return RapidFunction;
}

void Minuit2Wrapper::FixParameters( vector<double> fix_values, vector<string> ParameterNames )
{
	(void) fix_values;
	(void) ParameterNames;
}

//Use Migrad to minimise the given function
void Minuit2Wrapper::Minimise()
{
	function->SetSigma(nSigma);

	cout << "Minuit2 Starting Fit" << endl;

	//Minimise the wrapped function
	MnMigrad mig( *function, *( function->GetMnUserParameters() ), (unsigned)Quality );//MINUIT_QUALITY );

	//Retrieve the result of the fit
	minimum = new FunctionMinimum( mig( (unsigned)maxSteps, bestTolerance ) );//(int)MAXIMUM_MINIMISATION_STEPS, FINAL_GRADIENT_TOLERANCE );

	//Work out the fit status - possibly dodgy
	int fitStatus=0;
	if ( !minimum->HasCovariance() )
	{
		fitStatus = 0;
	}
	else if ( !minimum->HasAccurateCovar() )
	{
		fitStatus = 1;
	}
	else if ( minimum->HasMadePosDefCovar() )
	{
		fitStatus = 2;
	}
	else
	{
		fitStatus = 3;
	}

	//Output time information
	time_t timeNow;
	time(&timeNow);

	cout << endl << "Minuit2 MnMigrad finished:\tStatus: " << fitStatus << "\t\t" << ctime( &timeNow ) << endl;

	// May also want to run Hesse before the minimisation to get better estimate
	// of the error matrix.
	string NoHesse("NoHesse");
	if( StringProcessing::VectorContains( &Options, &NoHesse ) == -1 )
	{
		cout << "Minuit2 Starting MnHesse!" << endl;
		//      Finally Now call HESSE to properly calculate the error matrix
		MnHesse hesse(1);
		hesse( *function, *minimum, 100000);

		const MnUserParameters * minimisedParameters = &(minimum->UserParameters());

		for( unsigned int i=0; i< RapidFunction->GetParameterSet()->GetAllFloatNames().size(); ++i )
		{
			cout << RapidFunction->GetParameterSet()->GetAllFloatNames()[i] << "\t±\t" << minimisedParameters->Error( RapidFunction->GetParameterSet()->GetAllFloatNames()[i] ) << endl;
		}

		//Work out the fit status - possibly dodgy
		if ( !minimum->HasCovariance() )
		{
			fitStatus = 0;
		}
		else if ( !minimum->HasAccurateCovar() )
		{
			fitStatus = 1;
		}
		else if ( minimum->HasMadePosDefCovar() )
		{
			fitStatus = 2;
		}
		else
		{
			fitStatus = 3;
		}

		time(&timeNow);
		cout << endl << "Minuit2 MnHesse finished:\tStatus: " << fitStatus << "\t\t" << ctime( &timeNow ) << endl;
	}

	vector<double> allMin, allMax;

	string MinosOption("MinosErrors");
	if( StringProcessing::VectorContains( &Options, &MinosOption ) != -1 )
	{	
		cout << "Minuit2 Starting MnMinos!" << endl;
		MnMinos minos( *function, *minimum, 100000 );
		for( unsigned int i=0; i< RapidFunction->GetParameterSet()->GetAllFloatNames().size(); ++i )
		{
			MinosError thisErr = minos.Minos( i );
			cout << RapidFunction->GetParameterSet()->GetAllFloatNames()[i] << "\t+\t" << thisErr.Upper() << "\t-\t" << thisErr.Lower() << endl;

			allMin.push_back( thisErr.Lower() ); allMax.push_back( thisErr.Upper() );
		}
		//Work out the fit status - possibly dodgy
		if ( !minimum->HasCovariance() )
		{
			fitStatus = 0;
		}
		else if ( !minimum->HasAccurateCovar() )
		{
			fitStatus = 1;
		}
		else if ( minimum->HasMadePosDefCovar() )
		{
			fitStatus = 2;
		}
		else
		{
			fitStatus = 3;
		}

		time(&timeNow);
		cout << endl << "Minuit2 MnMinos finished:\tStatus: " << fitStatus << "\t\t" << ctime( &timeNow ) << endl;
	}

	//if( minimum != NULL ) delete minimum;
	//minimum = new FunctionMinimum( mig( (unsigned)maxSteps, bestTolerance ) );//(int)MAXIMUM_MINIMISATION_STEPS, FINAL_GRADIENT_TOLERANCE );

	// Need to add in the running of Hesse and Minos here. Should be configurable.

	//Output time information
	time(&timeNow);

	if ( !minimum->HasCovariance() )
	{
		fitStatus = 0;
	}
	else if ( !minimum->HasAccurateCovar() )
	{
		fitStatus = 1;
	}
	else if ( minimum->HasMadePosDefCovar() )
	{
		fitStatus = 2;
	}
	else
	{
		fitStatus = 3;
	}

	cout << endl << "Minuit2 finished:\tStatus: " << fitStatus << "\t\t" << ctime( &timeNow ) << endl;

	//Make a set of the fitted parameters
	const MnUserParameters * minimisedParameters = &(minimum->UserParameters());
	ParameterSet * newParameters = RapidFunction->GetParameterSet();
	vector<string> allNames = newParameters->GetAllNames();
	ResultParameterSet * fittedParameters = new ResultParameterSet( allNames );
	for (unsigned int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
	{
		string parameterName = allNames[nameIndex];
		PhysicsParameter * oldParameter = newParameters->GetPhysicsParameter( parameterName );
		double parameterValue = minimisedParameters->Value( parameterName.c_str() );
		double parameterError = minimisedParameters->Error( parameterName.c_str() );

		fittedParameters->SetResultParameter( parameterName, parameterValue, oldParameter->GetOriginalValue(), parameterError,
				oldParameter->GetMinimum(), oldParameter->GetMaximum(),
				oldParameter->GetType(), oldParameter->GetUnit() );
	}

	// Get the covariance matrix. Stored as an n*(n+1)/2 vector
	const MnUserCovariance * covMatrix = &(minimum->UserCovariance());
	vector<double> covData = covMatrix->Data();

	//Make a location to store the contour plots
	vector< FunctionContour* > allContours;
	vector< pair< int, int > > allIndices;
	int sigmaMax = 2;
	for (unsigned int plotIndex = 0; plotIndex < contours.size(); ++plotIndex )
	{
		//Find the index (in Minuit) of the parameter names requested for the plot
		int xParameterIndex = StringProcessing::VectorContains( &allNames, &( contours[plotIndex].first ) );
		int yParameterIndex = StringProcessing::VectorContains( &allNames, &( contours[plotIndex].second ) );
		if ( xParameterIndex == -1 )
		{
			cerr << "Contour plotting skipped: x parameter \"" << contours[plotIndex].first << "\" not found" << endl;
		}
		else if ( newParameters->GetPhysicsParameter( contours[plotIndex].first )->GetType() == "Fixed" )
		{
			cerr << "Contour plotting skipped: x parameter \"" << contours[plotIndex].first << "\" is fixed" << endl;
		}
		else if ( yParameterIndex == -1 )
		{
			cerr << "Contour plotting skipped: y parameter \"" << contours[plotIndex].second << "\" not found" << endl;
		}
		else if ( newParameters->GetPhysicsParameter( contours[plotIndex].second )->GetType() == "Fixed" )
		{
			cerr << "Contour plotting skipped: y parameter \"" << contours[plotIndex].second << "\" is fixed" << endl;
		}
		else
		{
			allIndices.push_back( make_pair( xParameterIndex, yParameterIndex ) );
			FunctionContour * newContour = new FunctionContour( contours[plotIndex].first, contours[plotIndex].second, sigmaMax );
			allContours.push_back(newContour);
		}
	}

	//Make the contour plots
	const MnContours contoursFromMinuit = MnContours( *function, *minimum );
	for ( int sigma = 1; sigma <= sigmaMax; ++sigma )
	{
		//Set the sigma value for the contours
		function->SetSigma(sigma);

		for (unsigned int plotIndex = 0; plotIndex < allIndices.size(); ++plotIndex )
		{
			//If the parameters have valid indices, ask minuit to plot them
			int numberOfPoints = 40;
			//			int iErrf;
			//			double xCoordinates[numberOfPoints], yCoordinates[numberOfPoints];
			vector< pair< double, double> > oneContour = contoursFromMinuit( unsigned(allIndices[plotIndex].first), unsigned(allIndices[plotIndex].second), unsigned(numberOfPoints) );
			allContours[plotIndex]->SetPlot( sigma, oneContour );
		}
	}

	PhysicsBottle* newBottle = RapidFunction->GetPhysicsBottle();
	fitResult = new FitResult( minimum->Fval(), fittedParameters, fitStatus, newBottle, NULL, allContours );

	vector<string> floated = RapidFunction->GetParameterSet()->GetAllFloatNames();
	for( unsigned int i=0; i<allMax.size(); ++i )
	{
		ResultParameterSet* thisResultSet = fitResult->GetResultParameterSet();
		ResultParameter* thisResult = thisResultSet->GetResultParameter( floated[i] );
		thisResult->SetAssymErrors( allMax[i], fabs(allMin[i]), 0.5*(fabs(allMax[i])+fabs(allMin[i])) );
	}
}

//Return the result of minimisation
FitResult * Minuit2Wrapper::GetFitResult()
{
	return fitResult;
}

//Request contour plots
void Minuit2Wrapper::ContourPlots( vector< pair< string, string > > ContourParameters )
{
	contours = ContourParameters;
}


//	The following 3 functions are simply coded up to still have the class compile correctly, but these should be implemented in the future
void Minuit2Wrapper::CallHesse()
{
	MnHesse hesse(1);
	hesse( *function, *minimum, 100000);
}

RapidFitMatrix* Minuit2Wrapper::GetCovarianceMatrix()
{
	return NULL;
}

void Minuit2Wrapper::ApplyCovarianceMatrix( RapidFitMatrix* Input )
{
	(void)Input;
	return;
}

void Minuit2Wrapper::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
}

void Minuit2Wrapper::SetNSigma( int input )
{
	nSigma = input;
}

