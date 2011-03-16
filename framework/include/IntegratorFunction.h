/**
	@class IntegratorFunction

	A wrapper to make IPDF usable by the numerical integrator

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-8
*/

#ifndef INTEGRATOR_FUNCTION_H
#define INTEGRATOR_FUNCTION_H

#include "Math/IFunction.h"
#include "TFoamIntegrand.h"
#include "IPDF.h"

using namespace ROOT::Math;

class IntegratorFunction : public IBaseFunctionMultiDim, public IBaseFunctionOneDim, public TFoamIntegrand
{
	public:
		IntegratorFunction();
		IntegratorFunction( IPDF*, DataPoint*, vector<string>, vector<string> );
		IntegratorFunction( IPDF*, DataPoint*, vector<string>, vector<string>, vector<double>, vector<double> );

		IPDF * GetWrappedFunction();

		//Interface functions
		virtual ~IntegratorFunction();
		virtual IntegratorFunction * Clone() const;
		virtual unsigned int NDim() const;
		virtual double operator()( double x ) const;
		virtual double operator()( const double * x ) const;
		virtual Double_t Density(Int_t ndim, Double_t *);
		virtual IBaseFunctionMultiDim & operator=( const IntegratorFunction& );

	private:
		//Interface function
		virtual double DoEval( double x ) const;
		virtual double DoEval( const double * x ) const;

		IPDF * wrappedFunction;
		DataPoint * currentPoint;
		vector<string> doIntegrate, dontIntegrate;
		vector<double> minima, ranges;
		vector<int> cache_positions;
};

#endif
