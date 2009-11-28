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
#include "Minuit2Function.h"
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
		virtual void ContourPlots( vector< pair< string, string > > );

	private:
		//MnMigrad minuit;
		Minuit2Function * function;
		FitResult * fitResult;
		vector< pair< string, string > > contours;
};

#endif
