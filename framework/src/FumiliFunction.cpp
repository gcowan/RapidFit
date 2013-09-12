// $Id: FumiliFunction.cpp,v 1.2 2009/11/13 09:57:06 gcowan Exp $
/**
        @class FumiliFunction

        A wrapper making IPDFs work with the Minuit2 API
	using the Fumili minimisation algorithm. This should
	only be used with NLL or chi2 minimisations.

        @author Greig A Cowan greig.cowan@cern.ch
	@date 2009-10-09
*/

//	RapidFit Headers
#include "FumiliFunction.h"
//	System Headers
#include <iostream>
#include <cmath>

using namespace::std;

//Default constructor
//FumiliFunction::FumiliFunction() : ParametricFunction(1), function(), parameters()
//{
//}

//Constructor with correct argument
FumiliFunction::FumiliFunction( IFitFunction* NewFitFunction, int NSigma ) : ParametricFunction( int(NewFitFunction->GetParameterSet()->GetAllNames().size()) ), function(NewFitFunction), parameters(), nSigma(NSigma)
{
	// Need to change this constructor since we now pass the numParams and not the fit function
	// Not entirely sure what to do here.

	parameters = new MnUserParameters();

	//Populate the parameters
	ParameterSet * newParameters = function->GetParameterSet();
	vector<string> allNames = newParameters->GetAllNames();
	vector<string>::iterator nameIterator;
	for (nameIterator = allNames.begin(); nameIterator != allNames.end(); nameIterator++)
	{
		PhysicsParameter * newParameter = newParameters->GetPhysicsParameter( *nameIterator );
		if ( newParameter->GetType() == "Fixed")
		{
			//Fixed parameter
			parameters->Add( nameIterator->c_str(), newParameter->GetValue() );
		}
		else// if ( newParameter->GetType() == "Free" || newParameter->GetType() == "GaussianConstrained" )
		{
			//Free parameter with limits
			parameters->Add( nameIterator->c_str(), newParameter->GetValue(),
					0.1, //This is the parameter uncertainty - i.e. I think it's a limit on MiGrad convergence
					newParameter->GetMinimum(), newParameter->GetMaximum() );
		}
		/*
		else
		{
			cerr << "Unrecognised parameter type" << endl;
		}*/
	}
}

//Destructor
FumiliFunction::~FumiliFunction()
{
}

double FumiliFunction::operator()( const vector<double>& NewParameterValues) const
{
	//Make parameter set and pass to wrapped function
	ParameterSet * newParameters = function->GetParameterSet();
	vector<string> allNames = newParameters->GetAllNames();
	for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); nameIndex++)
	{
		PhysicsParameter * newParameter = newParameters->GetPhysicsParameter( allNames[nameIndex] );
		newParameter->SetValue( NewParameterValues[nameIndex] );
		newParameters->SetPhysicsParameter( allNames[nameIndex], newParameter );
		//cout << allNames[nameIndex] << " " << NewParameterValues[nameIndex] << endl;
	}
	function->SetParameterSet( newParameters );

	//Return function value
	double value = function->Evaluate();
	//cout << function->Evaluate() << endl;
	if ( value < 0. or std::isnan(value) )
	{
		return 1.;
	}
	else
	{
		return value;
	}
}

double FumiliFunction::Up() const
{
	return function->UpErrorValue( nSigma );
}

//Return the parameters to minimise with
MnUserParameters * FumiliFunction::GetMnUserParameters()
{
	return parameters;
}
