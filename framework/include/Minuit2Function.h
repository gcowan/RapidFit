/**
        @class Minuit2Function

        A wrapper making IPDFs work with the Minuit2 API

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef MINUIT2_FUNCTION_H
#define MINUIT2_FUNCTION_H

//	ROOT Headers
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/ParametricFunction.h"
//	RapidFit Headers
#include "IFitFunction.h"

using namespace ROOT::Minuit2;

class Minuit2Function : public FCNBase
{
	public:
		//Minuit2Function();
		Minuit2Function( IFitFunction*, int );
		~Minuit2Function();

		void SetSigma(int);
		MnUserParameters * GetMnUserParameters();

		//Interface functions
		virtual double operator()( const vector<double>& ) const;
		virtual double Up() const;
		virtual void SetErrorDef( double );
		virtual double ErrorDef() const;
	
	private:
		//	Uncopyable!
		Minuit2Function ( const Minuit2Function& );
		Minuit2Function& operator = ( const Minuit2Function& );

		IFitFunction * function;
		MnUserParameters * parameters;
		double up;
};

#endif

