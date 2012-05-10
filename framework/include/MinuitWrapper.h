/**
        @class MinuitWrapper

        A wrapper to integrate Minuit with RapidFit

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef MINUIT_WRAPPER_H
#define MINUIT_WRAPPER_H

///	ROOT Headers
#include "TMinuit.h"
#include "TMatrixDSym.h"
///	RapidFit Headers
#include "IMinimiser.h"
#include "FitFunction.h"
#include "FitResult.h"
#include "FunctionContour.h"
#include "ParameterSet.h"
///	System Headers
#include <vector>
#include <string>

using namespace::std;

void Function( Int_t&, Double_t*, Double_t&, Double_t*, Int_t);

class MinuitWrapper : public IMinimiser
{
	public:
		MinuitWrapper();
		MinuitWrapper( int, int=0 );
		~MinuitWrapper();

		//Interface functions
		void SetOutputLevel( int );
                virtual void SetupFit( FitFunction* );
		virtual void FixParameters( vector<double>, vector<string> );
		virtual void Minimise();
		virtual FitResult * GetFitResult();
		virtual void ContourPlots( vector< pair< string, string > > );

		virtual void SetSteps( int );
		virtual void SetTolerance( double );
		virtual void SetOptions( vector<string> );
		virtual void SetQuality( int );
		virtual FitFunction* GetFitFunction();

		ResultParameterSet* GetResultParameters( vector<string> allNames, ParameterSet* );
		TMatrixDSym* GetCovarianceMatrix();
		vector<double> oldGetCovarianceMatrix( int numParams );
		virtual vector<FunctionContour*> ConstructContours( vector<string>, ParameterSet* );

		void CallHesse();

		void ApplyCovarianceMatrix( TMatrixDSym* Input );

	private:
		//	Uncopyable!
		MinuitWrapper ( const MinuitWrapper& );
		MinuitWrapper& operator = ( const MinuitWrapper& );

		static void Function( Int_t&, Double_t*, Double_t&, Double_t*, Int_t );
		static FitFunction * function;

		TMinuit * minuit;
		FitResult * fitResult;
		vector< pair< string, string > > contours;
		int print_verbosity;
		int maxSteps;
		double bestTolerance;
		vector<string> Options;
		int Quality;
		ParameterSet* test;
};

#endif

