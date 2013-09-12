// $Id: FumiliWrapper.cpp,v 1.2 2009/11/13 09:57:06 gcowan Exp $
/**
  @class FumiliWrapper

  A wrapper for Minuit2, implementing IMinimiser

  @author Greig A Cowan greig.cowan@cern.ch
  @date 2009-10-09
  */

//	ROOT Headers
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnFumiliMinimize.h"
#include "Minuit2/FumiliStandardMaximumLikelihoodFCN.h"
#include "TMatrixDSym.h"
//	RapidFit Headers
#include "FumiliWrapper.h"
#include "ResultParameterSet.h"
#include "PhysicsBottle.h"
#include "PhaseSpaceBoundary.h"
//	System Headers
#include <iostream>
#include <limits>

//const double MAXIMUM_MINIMISATION_STEPS = 800.0;
//const double FINAL_GRADIENT_TOLERANCE = 0.001;
const double STEP_SIZE = 0.01;
//const int MINUIT_QUALITY = 2;

//Default constructor
FumiliWrapper::FumiliWrapper() : function(NULL), RapidFunction(NULL), fitResult(NULL), contours(), maxSteps(), bestTolerance(), Options(), Quality(), debug(new DebugClass(false) ), nSigma(1)
{
}

//Destructor
FumiliWrapper::~FumiliWrapper()
{
	//delete minuit;
	if( debug != NULL ) delete debug;
}

void FumiliWrapper::SetSteps( int newSteps )
{
	maxSteps = newSteps;
}

void FumiliWrapper::SetTolerance( double newTolerance )
{
	bestTolerance = newTolerance;
}

void FumiliWrapper::SetOptions( vector<string> newOptions )
{
	Options = newOptions;
}

void FumiliWrapper::SetQuality( int newQuality )
{
	Quality = newQuality;
}

void FumiliWrapper::SetupFit( IFitFunction* NewFunction )
{
	// Instantiate the FumiliFunction (which is just a parametric function)
	// There is another constructor that takes a vector of parameters
	function = new FumiliFunction( NewFunction, nSigma );
	RapidFunction = NewFunction;
}

IFitFunction* FumiliWrapper::GetFitFunction()
{
	return RapidFunction;
}

void FumiliWrapper::FixParameters( vector<double> fix_values, vector<string> ParameterNames )
{
	(void) fix_values;
	(void) ParameterNames;
}

//Use Migrad to minimise the given function
void FumiliWrapper::Minimise()
{
	ParameterSet * newParameters = RapidFunction->GetParameterSet();
	vector<string> allNames = newParameters->GetAllNames();
	//	int numParams = allNames.size();

	/*
	// Fill a vector of doubles for each set of physics parameters that you
	// want to sample. What about case where Apara_sq + Aperp_sq > 1?
	vector< vector<double> > positions;
	PhysicsBottle* bottle = NewFunction->GetPhysicsBottle();
	ParameterSet* parameters = bottle->GetParameterSet();
	vector<string> names = parameters->GetAllNames();
	double nsteps = 1;
	for( int k = 0; k < names.size(); ++k)
	{
	for( int j = 0; j < nsteps; ++j)
	{
	vector<double> tempPos;
	for( int i = 0; i < names.size(); ++i )
	{
	double value;
	if ( i != k )
	{
	value = parameters->GetPhysicsParameter(names[i])->GetValue();
	}
	else
	{
	double min = parameters->GetPhysicsParameter(names[i])->GetMinimum();
	double max = parameters->GetPhysicsParameter(names[i])->GetMaximum();
	double step = (max - min)/nsteps;
	value = min + j * step;
	}
	tempPos.push_back(value);
	}
	positions.push_back(tempPos);
	}
	}
	*/

	// Fill a vector of doubles for each set of observables
	vector< vector<double> > positions;
	PhysicsBottle* bottle = RapidFunction->GetPhysicsBottle();
	PhaseSpaceBoundary* boundary = bottle->GetResultDataSet(0)->GetBoundary();
	vector<string> names = boundary->GetAllNames();

	vector<double> observableSteps;
	int nsteps = 10;

	// Could make this faster...
	for ( int step = 0; step < nsteps; ++step)
	{
		vector<double> tempPos;
		for( unsigned int observable = 0; observable < names.size(); ++observable )
		{
			if ( !boundary->GetConstraint(names[observable])->IsDiscrete() )
			{
				double min = boundary->GetConstraint(names[observable])->GetMinimum();
				double max = boundary->GetConstraint(names[observable])->GetMinimum();
				double delta = (max - min)/nsteps;
				double position = min + step*delta;
				tempPos.push_back( position );
			}
			else
			{
				double value = boundary->GetConstraint(names[observable])->CreateObservable()->GetValue();
				tempPos.push_back( value );
			}
		}
		positions.push_back(tempPos);
	}

	// Now, get the FumiliFCNBase function which will be passed to the Fumili minimiser
	FumiliStandardMaximumLikelihoodFCN fumFCN( *function, positions );

	// Setup the minimiser
	MnFumiliMinimize fumili( fumFCN, *( function->GetMnUserParameters() ), (unsigned)Quality);//MINUIT_QUALITY);

	// Do the minimisation
	FunctionMinimum minimum = fumili( (unsigned)maxSteps, bestTolerance );//(int)MAXIMUM_MINIMISATION_STEPS, FINAL_GRADIENT_TOLERANCE );

	// Once we have the FunctionMinimum, code same as in other Wrappers
	//Create the fit results
	const MnUserParameters * minimisedParameters = &minimum.UserParameters();
	ResultParameterSet * fittedParameters = new ResultParameterSet( allNames );
	for ( unsigned int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
	{
		string parameterName = allNames[nameIndex];
		PhysicsParameter * oldParameter = newParameters->GetPhysicsParameter( parameterName );
		double parameterValue = minimisedParameters->Value( parameterName.c_str() );
		double parameterError = minimisedParameters->Error( parameterName.c_str() );

		fittedParameters->SetResultParameter( parameterName, parameterValue, oldParameter->GetOriginalValue(), parameterError,
				-oldParameter->GetMinimum(), oldParameter->GetMaximum(),
				oldParameter->GetType(), oldParameter->GetUnit() );
	}

	int fitStatus;
	if ( !minimum.HasCovariance() )
	{
		fitStatus = 0;
	}
	else if ( !minimum.HasAccurateCovar() )
	{
		fitStatus = 1;
	}
	else if ( minimum.HasMadePosDefCovar() )
	{
		fitStatus = 2;
	}
	else
	{
		fitStatus = 3;
	}

	PhysicsBottle* newBottle = RapidFunction->GetPhysicsBottle();
	fitResult = new FitResult( minimum.Fval(), fittedParameters, fitStatus, newBottle  );
}

//Return the result of minimisation
FitResult * FumiliWrapper::GetFitResult()
{
	return fitResult;
}

//Request contour plots
void FumiliWrapper::ContourPlots( vector< pair< string, string > > ContourParameters )
{
	contours = ContourParameters;
}

//      The following 3 functions are simply coded up to still have the class compile correctly, but these should be implemented in the future
void FumiliWrapper::CallHesse()
{
	return;
}

RapidFitMatrix* FumiliWrapper::GetCovarianceMatrix()
{
	return NULL;
}

void FumiliWrapper::ApplyCovarianceMatrix( RapidFitMatrix* Input )
{
	(void)Input;
	return;
}

void FumiliWrapper::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
}

void FumiliWrapper::SetNSigma( int input )
{
	nSigma = input;
}

