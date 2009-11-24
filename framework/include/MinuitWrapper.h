/**
        @class MinuitWrapper

        A wrapper to integrate Minuit with RapidFit

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef MINUIT_WRAPPER_H
#define MINUIT_WRAPPER_H

#include "IMinimiser.h"
#include "TMinuit.h"
#include "FitFunction.h"
#include "FitResult.h"

void Function( int&, double*, double&, double*, int );

class MinuitWrapper : public IMinimiser
{
	public:
		MinuitWrapper();
		MinuitWrapper( int );
		~MinuitWrapper();

		//Interface functions
		virtual void Minimise( FitFunction* );
		virtual FitResult * GetFitResult();

	private:
		friend void Function( int&, double*, double&, double*, int );

		TMinuit * minuit;
		static FitFunction * function;
		//static vector<string> allNames;
		FitResult * fitResult;
};

#endif
