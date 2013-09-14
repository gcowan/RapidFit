/**
        @class NegativeLogLikelihoodThreadedNew

        A fit function with a multithreaded method for an NLL calculation

        @author Robert Currie rcurrie@cern.ch
	@date 2011-10
*/

#pragma once
#ifndef NEGATIVE_LOG_LIKELIHOOD_THREADED_NEW_H
#define NEGATIVE_LOG_LIKELIHOOD_THREADED_NEW_H

//	RapidFit Headers
#include "FitFunction.h"
#include "Threading.h"

class NegativeLogLikelihoodThreadedNew : public FitFunction
{
	public:
		NegativeLogLikelihoodThreadedNew();
		~NegativeLogLikelihoodThreadedNew();

		virtual double UpErrorValue(int);

	protected:
		virtual double EvaluateDataSet( IPDF*, IDataSet*, int );

};

#endif

