/**
        @class Minuit2Wrapper

        A wrapper for Minuit2, implementing IMinimiser

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef MINUIT_2_WRAPPER_H
#define MINUIT_2_WRAPPER_H

#include "IMinimiser.h"
#include "Minuit2/MnMigrad.h"
#include "MinuitFunction.h"
#include "FitResult.h"

using namespace ROOT::Minuit2;

class Minuit2Wrapper : public IMinimiser
{
	public:
		Minuit2Wrapper();
		~Minuit2Wrapper();

		//Interface functions
		virtual void Minimise( FitFunction* );
		virtual FitResult * GetFitResult();


	private:
		//MnMigrad minuit;
		MinuitFunction * function;
		FitResult * fitResult;
};

#endif
