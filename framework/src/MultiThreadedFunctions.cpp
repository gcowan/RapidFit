
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

vector<double>* MultiThreadedFunctions::ParallelEvaulate( IPDF* thisFunction, IDataSet* thesePoints, ThreadingConfig* threadingInfo )
{
	if( threadingInfo->MultiThreadingInstance == "pthreads" )
	{
		return MultiThreadedFunctions::ParallelEvaluate_pthreads( thisFunction, thesePoints, threadingInfo->numThreads, threadingInfo->wantedComponent );
	}
	else
	{
		return new vector<double>();
	}
}

vector<double>* MultiThreadedFunctions::ParallelEvaulate( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, ThreadingConfig* threadingInfo )
{
	if( threadingInfo->MultiThreadingInstance == "pthreads" )
	{
		return MultiThreadedFunctions::ParallelEvaluate_pthreads( thisFunction, thesePoints, threadingInfo->numThreads, threadingInfo->wantedComponent );
	}
	else
	{
		return new vector<double>();
	}
}

vector<double>* MultiThreadedFunctions::ParallelIntegrate( IPDF* thisFunction, IDataSet* thesePoints, PhaseSpaceBoundary* thisBoundary, ThreadingConfig* threadingInfo )
{
	if( threadingInfo->MultiThreadingInstance == "pthreads" )
	{
		return MultiThreadedFunctions::ParallelIntegrate_pthreads( thisFunction, thesePoints, thisBoundary, threadingInfo->numThreads );
	}
	else
	{
		return new vector<double>();
	}
}

vector<double>* MultiThreadedFunctions::ParallelIntegrate( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, vector<PhaseSpaceBoundary*> thisBoundary, ThreadingConfig* threadingInfo )
{
	if( threadingInfo->MultiThreadingInstance == "pthreads" )
	{
		return MultiThreadedFunctions::ParallelIntegrate_pthreads( thisFunction, thesePoints, thisBoundary, threadingInfo->numThreads );
	}
	else
	{
		return new vector<double>();
	}
}

vector<double>* MultiThreadedFunctions::ParallelEvaluate_pthreads( IPDF* thisFunction, IDataSet* thesePoints, unsigned int nThreads, ComponentRef* thisComponent )
{
	vector<IDataSet*> payLoad;
	vector<vector<DataPoint*> > datasets_data = Threading::divideData( thesePoints, nThreads );


	for( unsigned int i=0; i < nThreads; ++i )
	{
		payLoad.push_back( (IDataSet*) new MemoryDataSet( thesePoints->GetBoundary(), datasets_data[i] ) );
	}


	vector<IPDF*> functions;
	for( unsigned int i=0; i < nThreads; ++i )
	{
		functions.push_back( ClassLookUp::CopyPDF( thisFunction ) );
	}


	vector<double>* returnable = MultiThreadedFunctions::ParallelEvaluate_pthreads( functions, payLoad, nThreads, thisComponent );


	while( !payLoad.empty() )
	{
		if( payLoad.back() != NULL ) delete payLoad.back();

		payLoad.pop_back();
	}
	while( !functions.empty() )
	{
		if( functions.back() != NULL ) delete functions.back();
		functions.pop_back();
	}

	//for( unsigned int i=0; i< payLoad.size(); ++i ) if( payLoad[i] != NULL ) delete payLoad[i];
	//for( unsigned int i=0; i< functions.size(); ++i ) if( functions[i] != NULL ) delete functions[i];

	return returnable;
}

vector<double>* MultiThreadedFunctions::ParallelEvaluate_pthreads( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, unsigned int nThreads, ComponentRef* thisComponent )
{
	if( ( thisFunction.size() != thesePoints.size() ) || ( ( thesePoints.size() != nThreads ) || ( thisFunction.size() != nThreads ) ) )
	{
		cerr << "Catastrophy, exiting...." << endl;
		cout << "nThreads: " << nThreads << endl;
		cout << "thisFunction.size(): " << thisFunction.size() << endl;
		cout << "thesePoints.size(): " << thesePoints.size() << endl;
		exit(87356);
	}

	//      1 thread per core
	pthread_t* Thread = MultiThreadedFunctions::GetFittingThreads( nThreads );//new pthread_t[ nThreads ];
	pthread_attr_t attrib;

	//      Threads HAVE to be joinable
	//      We CANNOT _AND_SHOULD_NOT_ ***EVER*** return information to Minuit without the results from ALL threads successfully returned
	//      Not all pthread implementations are required to obey this as default when constructing the thread _SO_BE_EXPLICIT_
	pthread_attr_init(&attrib);
	pthread_attr_setdetachstate(&attrib, PTHREAD_CREATE_JOINABLE);

	Fitting_Thread* fit_thread_data = MultiThreadedFunctions::GetFittingThreadData( nThreads );//new Fitting_Thread[ (unsigned) nThreads ];

	for( unsigned int i=0; i < nThreads; ++i )
	{
		fit_thread_data[i].fittingPDF = thisFunction[i];
		fit_thread_data[i].dataSet = thesePoints[i];
		fit_thread_data[i].dataPoint_Result.clear();
		if( thisComponent != NULL ) fit_thread_data[i].thisComponent = new ComponentRef( *thisComponent );
	}

	//cout << "Creating Threads" << endl;

	//      Create the Threads and set them to be joinable
	for( unsigned int threadnum=0; threadnum< nThreads ; ++threadnum )
	{
		int status=0;
		if( thisComponent == NULL ) status = pthread_create( &Thread[threadnum], &attrib, MultiThreadedFunctions::Evaluate_pthread, (void *) &fit_thread_data[threadnum] );
		else status = pthread_create( &Thread[threadnum], &attrib, MultiThreadedFunctions::EvaluateComponent_pthread, (void *) &fit_thread_data[threadnum] );
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

	unsigned int size=0;
	for( unsigned int i=0; i< thesePoints.size(); ++i ) size+=thesePoints[i]->GetDataNumber();

	vector<double>* final_output = new vector<double>( size, 0. );

	unsigned int k=0;
	for( unsigned int i=0; i < nThreads; ++i )
	{
		for( unsigned j=0; j< fit_thread_data[i].dataPoint_Result.size(); ++j )
		{
			(*final_output)[k] = fit_thread_data[i].dataPoint_Result[j];
			++k;
		}
	}

	//delete[] fit_thread_data;
	//delete[] Thread;

	return final_output;
}

void* MultiThreadedFunctions::Evaluate_pthread( void *input_data )
{
	struct Fitting_Thread *thread_input = (struct Fitting_Thread*) input_data;

	double value=0;
	for( unsigned int i=0; i < (unsigned)thread_input->dataSet->GetDataNumber(); ++i )
	{
		//pthread_mutex_t* debug_lock = thread_input->fittingPDF->DebugMutex();
		try
		{
			value = thread_input->fittingPDF->Evaluate( thread_input->dataSet->GetDataPoint( i ) );
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

void* MultiThreadedFunctions::EvaluateComponent_pthread( void *input_data )
{
        struct Fitting_Thread *thread_input = (struct Fitting_Thread*) input_data;
 
        double value=0;
        for( unsigned int i=0; i < (unsigned)thread_input->dataSet->GetDataNumber(); ++i )
        {
                //pthread_mutex_t* debug_lock = thread_input->fittingPDF->DebugMutex();
                try
                {
                        value = thread_input->fittingPDF->EvaluateComponent( thread_input->dataSet->GetDataPoint( i ), thread_input->thisComponent );
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

vector<double>* MultiThreadedFunctions::ParallelIntegrate_pthreads( IPDF* thisFunction, IDataSet* thesePoints, PhaseSpaceBoundary* thisBoundary, unsigned int nThreads )
{
	vector<IDataSet*> payLoad;
	vector<vector<DataPoint*> > datasets_data = Threading::divideData( thesePoints, nThreads );

	for( unsigned int i=0; i < nThreads; ++i )
	{
		payLoad.push_back( (IDataSet*) new MemoryDataSet( thesePoints->GetBoundary(), datasets_data[i] ) );
	}

	double tempVal = thisFunction->Integral( payLoad[0]->GetDataPoint(0), thisBoundary );

	(void) tempVal;

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

	vector<double>* returnable = MultiThreadedFunctions::ParallelIntegrate_pthreads( functions, payLoad, boundaries, nThreads );

	for( unsigned int i=0; i< payLoad.size(); ++i ) if( payLoad[i] != NULL ) delete payLoad[i];
	for( unsigned int i=0; i< functions.size(); ++i ) if( functions[i] != NULL ) delete functions[i];

	return returnable;
}

void* MultiThreadedFunctions::Integrate_pthread( void *input_data )
{
	struct Fitting_Thread *thread_input = (struct Fitting_Thread*) input_data;

	double value=0;
	for( unsigned int i=0; i< (unsigned) thread_input->dataSet->GetDataNumber(); ++i )
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

vector<double>* MultiThreadedFunctions::ParallelIntegrate_pthreads( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, vector<PhaseSpaceBoundary*> theseBoundarys, unsigned int nThreads )
{

	if(        ( ( ( thisFunction.size() != thesePoints.size() ) || ( thesePoints.size() != theseBoundarys.size() ) ) || ( thisFunction.size() != theseBoundarys.size() ) )
			|| ( ( ( thisFunction.size() != nThreads ) || ( thesePoints.size() != nThreads ) ) || ( theseBoundarys.size() != nThreads ) )      )
	{
		cerr << "Catastrophy, exiting...." << endl;
		cout << "nThreads: " << nThreads << endl;
		cout << "thisFunction.size(): " << thisFunction.size() << endl;
		cout << "thesePoints.size(): " << thesePoints.size() << endl;
		cout << "theseBoundarys.size(): " << theseBoundarys.size() << endl;
		exit(87356);
	}

	for( unsigned int i=0; i< thisFunction.size(); ++i )
	{
		double tempVal = thisFunction[i]->Integral( thesePoints[i]->GetDataPoint(0), theseBoundarys[i] );
		(void) tempVal;
	}

	//      1 thread per core
	pthread_t* Thread = MultiThreadedFunctions::GetFittingThreads( nThreads );//new pthread_t[ nThreads ];
	pthread_attr_t attrib;

	//      Threads HAVE to be joinable
	//      We CANNOT _AND_SHOULD_NOT_ ***EVER*** return information to Minuit without the results from ALL threads successfully returned
	//      Not all pthread implementations are required to obey this as default when constructing the thread _SO_BE_EXPLICIT_
	pthread_attr_init(&attrib);
	pthread_attr_setdetachstate(&attrib, PTHREAD_CREATE_JOINABLE);

	Fitting_Thread* fit_thread_data = MultiThreadedFunctions::GetFittingThreadData( nThreads );//new Fitting_Thread[ nThreads ];

	for( unsigned int i=0; i< nThreads; ++i )
	{
		fit_thread_data[i].fittingPDF = thisFunction[i];
		fit_thread_data[i].dataSet = thesePoints[i];
		fit_thread_data[i].FitBoundary = theseBoundarys[i];
		fit_thread_data[i].dataPoint_Result.clear();
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

	unsigned int size=0;
	for( unsigned int i=0; i< thesePoints.size(); ++i ) size+=thesePoints[i]->GetDataNumber();

	vector<double>* final_output = new vector<double>( size, 0. );

	unsigned int k=0;
	for( unsigned int i=0; i< nThreads; ++i )
	{
		for( unsigned j=0; j< fit_thread_data[i].dataPoint_Result.size(); ++j )
		{
			(*final_output)[k] = fit_thread_data[i].dataPoint_Result[j];
			++k;
		}
	}

	//delete[] fit_thread_data;
	//delete[] Thread;
	return final_output;
}

Fitting_Thread* MultiThreadedFunctions::GetFittingThreadData( unsigned int nThreads )
{
	if( stored_tread_data == NULL )
	{
		stored_thread_data_n = 0;

		stored_tread_data = new Fitting_Thread[ nThreads ];
		stored_thread_data_n = nThreads;
		return stored_tread_data;
	}
	else
	{
		if( stored_thread_data_n == nThreads )
		{
			return stored_tread_data;
		}
		else
		{
			delete[] stored_tread_data;
			stored_tread_data = new Fitting_Thread[ nThreads ];
			stored_thread_data_n = nThreads;
			return stored_tread_data;
		}
	}
}

Fitting_Thread* MultiThreadedFunctions::stored_tread_data = NULL;

unsigned int MultiThreadedFunctions::stored_thread_data_n = 0;

pthread_t* MultiThreadedFunctions::GetFittingThreads( unsigned int nThreads )
{
	if( Thread == NULL )
	{
		Threads_n = 0;

		Thread = new pthread_t[ nThreads ];
		Threads_n = nThreads;
		return Thread;
	}
	else
	{
		if( Threads_n == nThreads )
		{
			return Thread;
		}
		else
		{
			delete[] Thread;
			Thread = new pthread_t[ nThreads ];
			Threads_n = nThreads;
			return Thread;
		}
	}
}

pthread_t* MultiThreadedFunctions::Thread = NULL;

unsigned int MultiThreadedFunctions::Threads_n = 0;

