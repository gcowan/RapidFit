
#pragma once
#ifndef _MULTI_THREADED_FUNCTIONS_H_
#define _MULTI_THREADED_FUNCTIONS_H_

#include "IPDF.h"
#include "IDataSet.h"
#include "DataPoint.h"
#include "Threading.h"
#include "ThreadingConfig.h"

class MultiThreadedFunctions
{
	public:

		static vector<double> ParallelEvaulate( IPDF* thisFunction, IDataSet* thesePoints, ThreadingConfig* threadingInfo );

		static vector<double> ParallelEvaulate( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, ThreadingConfig* threadingInfo );

		static vector<double> ParallelIntegrate( IPDF* thisFunction, IDataSet* thesePoints, PhaseSpaceBoundary* thisBoundary, ThreadingConfig* threadingInfo );

		static vector<double> ParallelIntegrate( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, vector<PhaseSpaceBoundary*> thisBoundary, ThreadingConfig* threadingInfo );
	private:

		MultiThreadedFunctions();
		~MultiThreadedFunctions();

		static vector<double> ParallelEvaluate_pthreads( IPDF* thisFunction, IDataSet* thesePoints, unsigned int nThreads );

		static vector<double> ParallelEvaluate_pthreads( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, unsigned int nThreads );

                #ifndef __CINT__
                        //      CINT behaves badly with this attribute
                        //      and,
                        //      g++ complains that this is a good place for it...
                        //      let's keep em happy
                        //
			static void* Evaluate_pthread( void *input_data ) __attribute__ ((noreturn));
			static void* Integrate_pthread( void *input_data ) __attribute__ ((noreturn));
		#else
			static void* Evaluate_pthread( void *input_data );
			static void* Integrate_pthread( void *input_data );
		#endif

		static vector<double> ParallelIntegrate_pthreads( IPDF* thisFunction, IDataSet* thesePoints, PhaseSpaceBoundary* thisBoundary, unsigned int nThreads );

		static vector<double> ParallelIntegrate_pthreads( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, vector<PhaseSpaceBoundary*> thisBoundary, unsigned int nThreads );

};

#endif

