/**
        @class MinuitFunction

        A wrapper making IPDFs work with the Minuit2 API

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef MINUIT_FUNCTION_H
#define MINUIT_FUNCTION_H

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/ParametricFunction.h"
#include "FitFunction.h"

using namespace ROOT::Minuit2;

class MinuitFunction : public FCNBase
//class MinuitFunction : public ParametricFunction
{
	public:
		MinuitFunction();
		MinuitFunction( FitFunction* );
		~MinuitFunction();

		MnUserParameters * GetMnUserParameters();

		//Interface functions
		virtual double operator()( const vector<double>& ) const;
		virtual double Up() const;
		virtual void SetErrorDef( double );
	
	protected:
		FitFunction * function;
		MnUserParameters * parameters;
		//int numParams;
};

#endif
