
#pragma once
#ifndef _MULTI_THREADED_FUNCTIONS_H_
#define _MULTI_THREADED_FUNCTIONS_H_

#include "IPDF.h"
#include "IDataSet.h"
#include "DataPoint.h"
#include "Threading.h"
#include "ThreadingConfig.h"

#include <stdio.h>
#include <pthread.h>
#include <vector>

#ifdef __CINT__
#undef __GNUC__
#define _SYS__SELECT_H_
struct pthread_t;
#undef __SYS__SELECT_H_
#define __GNUC__
#endif

using namespace::std;

class MultiThreadedFunctions
{
	public:

		static vector<double>* ParallelEvaulate( IPDF* thisFunction, IDataSet* thesePoints, ThreadingConfig* threadingInfo );

		static vector<double>* ParallelEvaulate( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, ThreadingConfig* threadingInfo );

		static vector<double>* ParallelIntegrate( IPDF* thisFunction, IDataSet* thesePoints, PhaseSpaceBoundary* thisBoundary, ThreadingConfig* threadingInfo );

		static vector<double>* ParallelIntegrate( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, vector<PhaseSpaceBoundary*> thisBoundary, ThreadingConfig* threadingInfo );
	private:

		MultiThreadedFunctions();
		~MultiThreadedFunctions();

		static vector<double>* ParallelEvaluate_pthreads( IPDF* thisFunction, IDataSet* thesePoints, unsigned int nThreads, ComponentRef* thisRef=NULL );

		static vector<double>* ParallelEvaluate_pthreads( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, unsigned int nThreads, ComponentRef* thisRef=NULL );

                #ifndef __CINT__
                        //      CINT behaves badly with this attribute
                        //      and,
                        //      g++ complains that this is a good place for it...
                        //      let's keep em happy
                        //
			static void* Evaluate_pthread( void *input_data ) __attribute__ ((noreturn));
			static void* EvaluateComponent_pthread( void *input_data ) __attribute__ ((noreturn));
			static void* Integrate_pthread( void *input_data ) __attribute__ ((noreturn));
		#else
			static void* Evaluate_pthread( void *input_data );
			static void* EvaluateComponent_pthread( void *input_data );
			static void* Integrate_pthread( void *input_data );
		#endif

		static Fitting_Thread* GetFittingThreadData( unsigned int nThreads );

		static Fitting_Thread* stored_tread_data;

		static unsigned int stored_thread_data_n;

		static pthread_t* GetFittingThreads( unsigned int nThreads );

		static pthread_t* Thread;

		static unsigned int Threads_n;

		static vector<double>* ParallelIntegrate_pthreads( IPDF* thisFunction, IDataSet* thesePoints, PhaseSpaceBoundary* thisBoundary, unsigned int nThreads );

		static vector<double>* ParallelIntegrate_pthreads( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, vector<PhaseSpaceBoundary*> thisBoundary, unsigned int nThreads );

};

#endif

