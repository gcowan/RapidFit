
#include "IPDF.h"
#include "IDataSet.h"
#include "DataPoint.h"
#include "Threading.h"
#include "MultiThreadedFunctions.h"
#include "ClassLookUp.h"
#include "MemoryDataSet.h"

#include <string>
#include <float.h>

using namespace::std;

vector<double> MultiThreadedFunctions::ParallelEvaulate( IPDF* thisFunction, IDataSet* thesePoints, ThreadingConfig* threadingInfo )
{
	vector<double> results;

	if( threadingInfo->MultiThreadingInstance == "pthreads" )
	{
		results = ParallelEvaluate_pthreads( thisFunction, thesePoints, threadingInfo->numThreads );
	}

	return results;
}

vector<double> MultiThreadedFunctions::ParallelEvaulate( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, ThreadingConfig* threadingInfo )
{
	vector<double> results;

	if( threadingInfo->MultiThreadingInstance == "pthreads" )
	{
		results = ParallelEvaluate_pthreads( thisFunction, thesePoints, threadingInfo->numThreads );
	}

	return results;
}

vector<double> MultiThreadedFunctions::ParallelIntegrate( IPDF* thisFunction, IDataSet* thesePoints, PhaseSpaceBoundary* thisBoundary, ThreadingConfig* threadingInfo )
{
	vector<double> results;

	if( threadingInfo->MultiThreadingInstance == "pthreads" )
	{
		results = ParallelIntegrate_pthreads( thisFunction, thesePoints, thisBoundary, threadingInfo->numThreads );
	}

	return results;
}

vector<double> MultiThreadedFunctions::ParallelIntegrate( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, vector<PhaseSpaceBoundary*> thisBoundary, ThreadingConfig* threadingInfo )
{
	vector<double> results;

	if( threadingInfo->MultiThreadingInstance == "pthreads" )
	{
		results = ParallelIntegrate_pthreads( thisFunction, thesePoints, thisBoundary, threadingInfo->numThreads );
	}

	return results;
}

vector<double> MultiThreadedFunctions::ParallelEvaluate_pthreads( IPDF* thisFunction, IDataSet* thesePoints, unsigned int nThreads )
{
	vector<IDataSet*> payLoad;
	vector<vector<DataPoint*> > datasets_data = Threading::divideData( thesePoints, nThreads );

	for( unsigned int i=0; datasets_data.size(); ++i )
	{
		payLoad.push_back( (IDataSet*) new MemoryDataSet( thesePoints->GetBoundary(), datasets_data[i] ) );
	}

	vector<IPDF*> functions;
	for( unsigned int i=0; i< nThreads; ++i )
	{
		functions.push_back( ClassLookUp::CopyPDF( thisFunction ) );
	}

	vector<double> returnable = MultiThreadedFunctions::ParallelEvaluate_pthreads( functions, payLoad, nThreads );

	for( unsigned int i=0; i< payLoad.size(); ++i ) if( payLoad[i] != NULL ) delete payLoad[i];
	for( unsigned int i=0; i< functions.size(); ++i ) if( functions[i] != NULL ) delete functions[i];

	return returnable;
}

vector<double> MultiThreadedFunctions::ParallelEvaluate_pthreads( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, unsigned int nThreads )
{
        if( ( thisFunction.size() != thesePoints.size() ) || ( ( thesePoints.size() != nThreads ) || ( thisFunction.size() != nThreads ) ) )
        {
                cerr << "Catastrophy, exiting...." << endl;
                exit(87356);
        }

	//      1 thread per core
	pthread_t* Thread = new pthread_t[ nThreads ];
	pthread_attr_t attrib;

	//      Threads HAVE to be joinable
	//      We CANNOT _AND_SHOULD_NOT_ ***EVER*** return information to Minuit without the results from ALL threads successfully returned
	//      Not all pthread implementations are required to obey this as default when constructing the thread _SO_BE_EXPLICIT_
	pthread_attr_init(&attrib);
	pthread_attr_setdetachstate(&attrib, PTHREAD_CREATE_JOINABLE);

	Fitting_Thread* fit_thread_data = new Fitting_Thread[ (unsigned) nThreads ];

	for( unsigned int i=0; i< nThreads; ++i )
	{
		fit_thread_data[i].fittingPDF = thisFunction[i];
		fit_thread_data[i].dataSet = thesePoints[i];
		fit_thread_data[i].dataPoint_Result = vector<double>();
	}

	//cout << "Creating Threads" << endl;

	//      Create the Threads and set them to be joinable
	for( unsigned int threadnum=0; threadnum< nThreads ; ++threadnum )
	{
		int status = pthread_create( &Thread[threadnum], &attrib, MultiThreadedFunctions::Evaluate_pthread, (void *) &fit_thread_data[threadnum] );
		if( status )
		{
			cerr << "ERROR:\tfrom pthread_create()\t" << status << "\t...Exiting\n" << endl;
			exit(-1);
		}
	}

	//cout << "Joining Threads!!" << endl;

	//      Join the Threads
	for( unsigned int threadnum=0; threadnum< nThreads ; ++threadnum )
	{
		int status = pthread_join( Thread[threadnum], NULL);
		if( status )
		{
			cerr << "Error Joining a Thread:\t" << threadnum << "\t:\t" << status << "\t...Exiting\n" << endl;
		}
	}

	//      Do some cleaning Up
	pthread_attr_destroy(&attrib);

	vector<double> final_output;

	for( unsigned int i=0; i< nThreads; ++i )
	{
		for( unsigned j=0; j< fit_thread_data[i].dataPoint_Result.size(); ++j )
		{
			final_output.push_back( fit_thread_data[i].dataPoint_Result[j] );
		}
	}

	delete fit_thread_data;
	delete Thread;

	return final_output;
}

void* MultiThreadedFunctions::Evaluate_pthread( void *input_data )
{
	struct Fitting_Thread *thread_input = (struct Fitting_Thread*) input_data;

	vector<double> values;
	double value=0;
	for( int i=0; i< thread_input->dataSet->GetDataNumber(); ++i )
	{
		//pthread_mutex_t* debug_lock = thread_input->fittingPDF->DebugMutex();
		try
		{
			value = thread_input->fittingPDF->Evaluate( thread_input->dataSet->GetDataPoint(i) );
		}
		catch( ... )
		{
			value = DBL_MAX;
		}

		thread_input->dataPoint_Result.push_back( value );
	}

	//      Finished evaluating this thread
	pthread_exit( NULL );
}

vector<double> MultiThreadedFunctions::ParallelIntegrate_pthreads( IPDF* thisFunction, IDataSet* thesePoints, PhaseSpaceBoundary* thisBoundary, unsigned int nThreads )
{
	vector<IDataSet*> payLoad;
	vector<vector<DataPoint*> > datasets_data = Threading::divideData( thesePoints, nThreads );

	for( unsigned int i=0; datasets_data.size(); ++i )
	{
		payLoad.push_back( (IDataSet*) new MemoryDataSet( thesePoints->GetBoundary(), datasets_data[i] ) );
	}

	vector<IPDF*> functions;
	for( unsigned int i=0; i< nThreads; ++i )
	{
		functions.push_back( ClassLookUp::CopyPDF( thisFunction ) );
	}

	vector<PhaseSpaceBoundary*> boundaries;
	for( unsigned int i=0; i< nThreads; ++i )
	{
		boundaries.push_back( new PhaseSpaceBoundary( *thisBoundary ) );
	}

	vector<double> returnable = MultiThreadedFunctions::ParallelIntegrate_pthreads( functions, payLoad, boundaries, nThreads );

	for( unsigned int i=0; i< payLoad.size(); ++i ) if( payLoad[i] != NULL ) delete payLoad[i];
	for( unsigned int i=0; i< functions.size(); ++i ) if( functions[i] != NULL ) delete functions[i];

	return returnable;
}

void* MultiThreadedFunctions::Integrate_pthread( void *input_data )
{
	struct Fitting_Thread *thread_input = (struct Fitting_Thread*) input_data;

	vector<double> values;
	double value=0;
	for( int i=0; i< thread_input->dataSet->GetDataNumber(); ++i )
	{
		//pthread_mutex_t* debug_lock = thread_input->fittingPDF->DebugMutex();
		try
		{
			value = thread_input->fittingPDF->Integral( thread_input->dataSet->GetDataPoint(i), thread_input->FitBoundary );
		}
		catch( ... )
		{
			value = DBL_MAX;
		}

		thread_input->dataPoint_Result.push_back( value );
	}

	//      Finished evaluating this thread
	pthread_exit( NULL );
}

vector<double> MultiThreadedFunctions::ParallelIntegrate_pthreads( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, vector<PhaseSpaceBoundary*> theseBoundarys, unsigned int nThreads )
{

	if(        ( ( ( thisFunction.size() != thesePoints.size() ) || ( thesePoints.size() != theseBoundarys.size() ) ) || ( thisFunction.size() != theseBoundarys.size() ) )
		|| ( ( ( thisFunction.size() != nThreads ) || ( thesePoints.size() != nThreads ) ) || ( theseBoundarys.size() != nThreads ) )      )
	{
		cerr << "Catastrophy, exiting...." << endl;
		exit(87356);
	}

	//      1 thread per core
	pthread_t* Thread = new pthread_t[ nThreads ];
	pthread_attr_t attrib;

	//      Threads HAVE to be joinable
	//      We CANNOT _AND_SHOULD_NOT_ ***EVER*** return information to Minuit without the results from ALL threads successfully returned
	//      Not all pthread implementations are required to obey this as default when constructing the thread _SO_BE_EXPLICIT_
	pthread_attr_init(&attrib);
	pthread_attr_setdetachstate(&attrib, PTHREAD_CREATE_JOINABLE);

	Fitting_Thread* fit_thread_data = new Fitting_Thread[ nThreads ];

	for( unsigned int i=0; i< nThreads; ++i )
	{
		fit_thread_data[i].fittingPDF = thisFunction[i];
		fit_thread_data[i].dataSet = thesePoints[i];
		fit_thread_data[i].FitBoundary = theseBoundarys[i];
		fit_thread_data[i].dataPoint_Result = vector<double>();
	}

	//cout << "Creating Threads" << endl;

	//      Create the Threads and set them to be joinable
	for( unsigned int threadnum=0; threadnum< nThreads ; ++threadnum )
	{
		int status = pthread_create( &Thread[threadnum], &attrib, MultiThreadedFunctions::Integrate_pthread, (void *) &fit_thread_data[threadnum] );
		if( status )
		{
			cerr << "ERROR:\tfrom pthread_create()\t" << status << "\t...Exiting\n" << endl;
			exit(-1);
		}
	}

	//cout << "Joining Threads!!" << endl;

	//      Join the Threads
	for( unsigned int threadnum=0; threadnum< nThreads ; ++threadnum )
	{
		int status = pthread_join( Thread[threadnum], NULL);
		if( status )
		{
			cerr << "Error Joining a Thread:\t" << threadnum << "\t:\t" << status << "\t...Exiting\n" << endl;
		}
	}

	//      Do some cleaning Up
	pthread_attr_destroy(&attrib);

	vector<double> final_output;

	for( unsigned int i=0; i< nThreads; ++i )
	{
		for( unsigned j=0; j< fit_thread_data[i].dataPoint_Result.size(); ++j )
		{
			final_output.push_back( fit_thread_data[i].dataPoint_Result[j] );
		}
	}

	delete fit_thread_data;
	delete Thread;

	return final_output;
}

