/**
        @class NegativeLogLikelihoodThreaded

        A fit function with a multithreaded method for an NLL calculation

        @author Robert Currie rcurrie@cern.ch
	@date 2011-10
*/

#pragma once
#ifndef NEGATIVE_LOG_LIKELIHOOD_THREADED_H
#define NEGATIVE_LOG_LIKELIHOOD_THREADED_H

//	RapidFit Headers
#include "FitFunction.h"
#include "Threading.h"

class NegativeLogLikelihoodThreaded : public FitFunction
{
	public:
		NegativeLogLikelihoodThreaded();
		~NegativeLogLikelihoodThreaded();

		virtual double UpErrorValue(int);

	protected:
		virtual double EvaluateDataSet( IPDF*, IDataSet*, int );

	private:
		#ifndef __CINT__
			//	CINT behaves badly with this attribute
			//	and,
			//	g++ complains that this is a good place for it...
			//	let's keep em happy
			//	
			static void* ThreadWork( void* ) __attribute__ ((noreturn));
		#else
			static void* ThreadWork( void* );
		#endif
};

#endif

