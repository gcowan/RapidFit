/**
        @class MinuitWrapper

        A wrapper to integrate Minuit with RapidFit

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef MINUIT_WRAPPER_H
#define MINUIT_WRAPPER_H

//	ROOT Headers
#include "TMinuit.h"
//	RapidFit Headers
#include "IMinimiser.h"
#include "FitFunction.h"
#include "FitResult.h"

void Function( int&, double*, double&, double*, int );

class MinuitWrapper : public IMinimiser
{
	public:
		MinuitWrapper();
		MinuitWrapper( int, int=0 );
		~MinuitWrapper();

		//Interface functions
		void SetOutputLevel( int );
		virtual void Minimise( FitFunction* );
		virtual FitResult * GetFitResult();
		virtual void ContourPlots( vector< pair< string, string > > );

		virtual void SetSteps( int );
		virtual void SetTolerance( double );
		virtual void SetOptions( vector<string> );
		virtual void SetQuality( int );

	private:
		//	Uncopyable!
		MinuitWrapper ( const MinuitWrapper& );
		MinuitWrapper& operator = ( const MinuitWrapper& );

		static void Function( int&, double*, double&, double*, int );
		static FitFunction * function;

		TMinuit * minuit;
		FitResult * fitResult;
		vector< pair< string, string > > contours;
		int print_verbosity;
		int maxSteps;
		double bestTolerance;
		vector<string> Options;
		int Quality;
};

#endif
