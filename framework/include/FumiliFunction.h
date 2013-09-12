// $Id: FumiliFunction.h,v 1.2 2009/11/13 09:57:05 gcowan Exp $
/*!
 * @class FumiliFunction
 *
 * @brief A wrapper making IPDFs work with the Minuit2 API
 *        using the Fumili minimisation method.
 *
 * @author Greig A Cowan greig.cowan@cern.ch
 */

#pragma once
#ifndef FUMILI_FUNCTION_H
#define FUMILI_FUNCTION_H

///	ROOT Headers
#include "Minuit2/ParametricFunction.h"
#include "Minuit2/MnUserParameters.h"
///	RapidFit Headers
#include "IFitFunction.h"

using namespace ROOT::Minuit2;

class FumiliFunction : public ParametricFunction
{
	public:
		//FumiliFunction();
		FumiliFunction( IFitFunction*, int );
		//FumiliFunction( int );
		~FumiliFunction();

		MnUserParameters * GetMnUserParameters();

		//Interface functions
		// Need to have this operator here, even if I don't implement it in the src.	
		virtual double operator()( const vector<double>& ) const;
		//virtual double operator()( const vector<double>&, const vector<double>& ) const;
		virtual double Up() const;		

	protected:
		//	Uncopyable!
		FumiliFunction ( const FumiliFunction& );
		FumiliFunction& operator = ( const FumiliFunction& );

		IFitFunction * function;
		MnUserParameters * parameters;
		int nSigma;
};

#endif

