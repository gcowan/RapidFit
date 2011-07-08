/**
        @class Minuit2Wrapper

        A wrapper for Minuit2, implementing IMinimiser

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef MINUIT_2_WRAPPER_H
#define MINUIT_2_WRAPPER_H

//	ROOT Headers
#include "Minuit2/MnMigrad.h"
#include "Minuit2Function.h"
//	RapidFit Headers
#include "IMinimiser.h"
#include "FitResult.h"

using namespace ROOT::Minuit2;

class Minuit2Wrapper : public IMinimiser
{
	public:
		Minuit2Wrapper();
		~Minuit2Wrapper();

		//Interface functions
		virtual void SetOutputLevel( int ){};
		virtual void Minimise( FitFunction* );
		virtual FitResult * GetFitResult();
		virtual void ContourPlots( vector< pair< string, string > > );

		virtual void SetSteps( int );
                virtual void SetTolerance( double );
                virtual void SetOptions( vector<string> );
		virtual void SetQuality( int );

	private:
		//	Uncopyable!
		Minuit2Wrapper ( const Minuit2Wrapper& );
		Minuit2Wrapper& operator = ( const Minuit2Wrapper& );

		//MnMigrad minuit;
		Minuit2Function * function;
		FitResult * fitResult;
		vector< pair< string, string > > contours;
                int maxSteps;
                double bestTolerance;
                vector<string> Options;
		int Quality;
};

#endif
