/**
        @class MinuitFunction

        A wrapper making IPDFs work with the Minuit2 API

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "MinuitFunction.h"
#include <iostream>

//This is the parameter uncertainty - i.e. I think it's a limit on MiGrad convergence
const double ERROR_THINGY = 0.1;

//Default constructor
MinuitFunction::MinuitFunction()
{
}

//Constructor with correct argument
MinuitFunction::MinuitFunction( FitFunction* NewFitFunction ) : function(NewFitFunction)
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
		else if ( newParameter->GetType() == "Unbounded" )
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
MinuitFunction::~MinuitFunction()
{
}

//Return the value to minimise, given the parameters passed
double MinuitFunction::operator()( const vector<double>& NewParameterValues ) const
{
	//Make parameter set and pass to wrapped function
	ParameterSet * temporaryParameters = function->GetParameterSet();
	if ( temporaryParameters->SetPhysicsParameters(NewParameterValues) )
	{
		function->SetParameterSet(temporaryParameters);
		return function->Evaluate();
	}
	else
	{
		cerr << "Minuit does not provide the correct parameters" << endl;
		return 0.0;
	}
}

//Return the up value for error calculation
double MinuitFunction::Up() const
{
	return function->UpErrorValue();
}

//Return the parameters to minimise with
MnUserParameters * MinuitFunction::GetMnUserParameters()
{
	return parameters;
}
