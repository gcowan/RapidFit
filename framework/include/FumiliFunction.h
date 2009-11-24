// $Id: FumiliFunction.h,v 1.2 2009/11/13 09:57:05 gcowan Exp $
/**
        @class FumiliFunction

        A wrapper making IPDFs work with the Minuit2 API
	using the Fumili minimisation method.

        @author Greig A Cowan greig.cowan@cern.ch
	@date 2009-10-09
*/

#ifndef FUMILI_FUNCTION_H
#define FUMILI_FUNCTION_H

#include "Minuit2/ParametricFunction.h"
#include "Minuit2/MnUserParameters.h"
#include "FitFunction.h"

using namespace ROOT::Minuit2;

class FumiliFunction : public ParametricFunction
{
	public:
		FumiliFunction();
		FumiliFunction( FitFunction* );
		//FumiliFunction( int );
		~FumiliFunction();

		MnUserParameters * GetMnUserParameters();

		//Interface functions
		// Need to have this operator here, even if I don't implement it in the src.	
		virtual double operator()( const vector<double>& ) const;
		//virtual double operator()( const vector<double>&, const vector<double>& ) const;
		virtual double Up() const;		

	protected:
		FitFunction * function;
		MnUserParameters * parameters;
};

#endif
