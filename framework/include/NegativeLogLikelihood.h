/**
        @class NegativeLogLikelihood

        A fit function with evaulate methods for an NLL calculation

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef NEGATIVE_LOG_LIKELIHOOD_H
#define NEGATIVE_LOG_LIKELIHOOD_H

//	RapidFit Headers
#include "FitFunction.h"

class NegativeLogLikelihood : public FitFunction
{
	public:
		NegativeLogLikelihood();
		~NegativeLogLikelihood();

		virtual double UpErrorValue(int);

	protected:
		virtual double EvaluateDataSet( IPDF*, IDataSet*, int );
};
#endif

