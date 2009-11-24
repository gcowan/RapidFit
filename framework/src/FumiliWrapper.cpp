// $Id: FumiliWrapper.cpp,v 1.2 2009/11/13 09:57:06 gcowan Exp $
/**
        @class FumiliWrapper

        A wrapper for Minuit2, implementing IMinimiser

        @author Greig A Cowan greig.cowan@cern.ch
	@date 2009-10-09
*/

#include "FumiliWrapper.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnFumiliMinimize.h"
#include "Minuit2/FumiliStandardMaximumLikelihoodFCN.h"
#include <iostream>
#include <limits>
#include "ResultParameterSet.h"
#include "PhysicsBottle.h"
#include "PhaseSpaceBoundary.h"

const double MAXIMUM_MINIMISATION_STEPS = 800.0;
const double FINAL_GRADIENT_TOLERANCE = 0.001;
const double STEP_SIZE = 0.01;
const int MINUIT_QUALITY = 2;

//Default constructor
FumiliWrapper::FumiliWrapper()
{
}

//Destructor
FumiliWrapper::~FumiliWrapper()
{
	//delete minuit;
}

//Use Migrad to minimise the given function
void FumiliWrapper::Minimise( FitFunction * NewFunction )
{
	ParameterSet * newParameters = NewFunction->GetParameterSet();
        vector<string> allNames = newParameters->GetAllNames();
	int numParams = allNames.size();

	// Instantiate the FumiliFunction (which is just a parametric function)
	// There is another constructor that takes a vector of parameters
	function = new FumiliFunction( NewFunction );

	/*
	// Fill a vector of doubles for each set of physics parameters that you
	// want to sample. What about case where Apara_sq + Aperp_sq > 1? 
	vector< vector<double> > positions;
	PhysicsBottle* bottle = NewFunction->GetPhysicsBottle();
	ParameterSet* parameters = bottle->GetParameterSet();	
	vector<string> names = parameters->GetAllNames();
	double nsteps = 1;
	for( int k = 0; k < names.size(); k++)
	{
		for( int j = 0; j < nsteps; j++)
		{
			vector<double> tempPos;
			for( int i = 0; i < names.size(); i++ )
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

	// Fill a vector of doubles for each set of observables that you
        vector< vector<double> > positions;
	PhysicsBottle* bottle = NewFunction->GetPhysicsBottle();
	//PhaseSpaceBoundary* boundary = bottle->GetResultPDF(0)->GetPhaseSpaceBoundary();
        vector<string> names;// = boundary->GetAllNames();   
	int nsteps = 10;
	for ( int j = 0; j < nsteps; j++)
	{
        	vector<double> tempPos;
                for( int i = 0; i < names.size(); i++ )
                {
			double value = 0.;
			//double min = boundary->
			tempPos.push_back(value);
		}
		positions.push_back(tempPos);
	}


	// Now, get the FumiliFCNBase function which will be passed to the Fumili minimiser
	FumiliStandardMaximumLikelihoodFCN fumFCN( *function, positions );
	
	// Setup the minimiser
	MnFumiliMinimize fumili( fumFCN, *( function->GetMnUserParameters() ), MINUIT_QUALITY);

	// Do the minimisation
	FunctionMinimum minimum = fumili( (int)MAXIMUM_MINIMISATION_STEPS, FINAL_GRADIENT_TOLERANCE );

	// Once we have the FunctionMinimum, code same as in other Wrappers
	//Create the fit results
	const MnUserParameters * minimisedParameters = &minimum.UserParameters();
	ResultParameterSet * fittedParameters = new ResultParameterSet( allNames );
	for (int nameIndex = 0; nameIndex < allNames.size(); nameIndex++)
	{
		string parameterName = allNames[nameIndex];
		PhysicsParameter * oldParameter = newParameters->GetPhysicsParameter( parameterName );
		double parameterValue = minimisedParameters->Value( parameterName.c_str() );
		double parameterError = minimisedParameters->Error( parameterName.c_str() );

		fittedParameters->SetResultParameter( parameterName, parameterValue, oldParameter->GetOriginalValue(), parameterError,
			       -numeric_limits<double>::max(), numeric_limits<double>::max(),
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

	fitResult = new FitResult( minimum.Fval(), fittedParameters, fitStatus, *( NewFunction->GetPhysicsBottle() ) );
}

//Return the result of minimisation
FitResult * FumiliWrapper::GetFitResult()
{
	return fitResult;
}
