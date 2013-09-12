/**
        @class MinuitFunction

        A wrapper making IPDFs work with the Minuit2 API

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "Minuit2Function.h"
//	System Headers
#include <iostream>

//This is the parameter uncertainty - i.e. I think it's a limit on MiGrad convergence
const double ERROR_THINGY = 0.1;

//Default constructor
//Minuit2Function::Minuit2Function() : function(), parameters(), up(1.)
//{
//}

//Constructor with correct argument
Minuit2Function::Minuit2Function( IFitFunction * NewFitFunction, int nSigma ) : function(NewFitFunction), parameters(), up( NewFitFunction->UpErrorValue( nSigma ) )
{
	parameters = new MnUserParameters();

	//Populate the parameters
	ParameterSet * newParameters = function->GetParameterSet();
	vector<string> allNames = newParameters->GetAllNames();
	vector<string>::iterator nameIterator;
	for (nameIterator = allNames.begin(); nameIterator != allNames.end(); nameIterator++)
	{
		PhysicsParameter * newParameter = newParameters->GetPhysicsParameter( *nameIterator );
		if ( newParameter->GetType() == "Fixed" )
		{
			//Fixed parameter
			parameters->Add( nameIterator->c_str(), newParameter->GetValue() );
		}
		else if ( newParameter->GetType() == "Unbounded" || newParameter->GetType() == "GaussianConstrained" )
		{
			//Unbounded parameter
			parameters->Add( nameIterator->c_str(), newParameter->GetValue(), ERROR_THINGY );
		}
		else
		{
			//Free parameter with limits
			parameters->Add( nameIterator->c_str(), newParameter->GetValue(), ERROR_THINGY,
					newParameter->GetMinimum(), newParameter->GetMaximum() );
		}
	}
}

//Destructor
Minuit2Function::~Minuit2Function()
{
}

//Return the value to minimise, given the parameters passed
double Minuit2Function::operator()( const vector<double>& NewParameterValues ) const
{
	//Make parameter set and pass to wrapped function
	ParameterSet * temporaryParameters = function->GetParameterSet();

	try
	{
		temporaryParameters->SetPhysicsParameters(NewParameterValues);
		function->SetParameterSet(temporaryParameters);
		return function->Evaluate();
	}
	catch(...)
	{
		cerr << "Minuit2 does not provide the correct parameters" << endl;
		throw(-9999);
	}
}

//Set the up value for error calculation
void Minuit2Function::SetErrorDef( double thisUp )
{
	up = thisUp;
}
void Minuit2Function::SetSigma( int Sigma )
{
	up = function->UpErrorValue(Sigma);
}

//Return the up value for error calculation
double Minuit2Function::Up() const
{
	return up;
}
double Minuit2Function::ErrorDef() const
{
	return up;
}

//Return the parameters to minimise with
MnUserParameters * Minuit2Function::GetMnUserParameters()
{
	return parameters;
}

