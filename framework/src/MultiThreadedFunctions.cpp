
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

vector<double>* MultiThreadedFunctions::ParallelEvaluate( IPDF* thisFunction, IDataSet* thesePoints, ThreadingConfig* threadingInfo )
{
	if( threadingInfo == NULL ) return MultiThreadedFunctions::ParallelEvaluate_pthreads( thisFunction, thesePoints, 4, NULL );
	if( threadingInfo->MultiThreadingInstance == "pthreads" )
	{
		return MultiThreadedFunctions::ParallelEvaluate_pthreads( thisFunction, thesePoints, threadingInfo->numThreads, threadingInfo->wantedComponent );
	}
	else
	{
		return new vector<double>();
	}
}

vector<double>* MultiThreadedFunctions::ParallelEvaluate( vector<IPDF*> thisFunction, vector<IDataSet*> thesePoints, ThreadingConfig* threadingInfo )
{
	if( threadingInfo == NULL ) return MultiThreadedFunctions::ParallelEvaluate_pthreads( thisFunction, thesePoints, 4, NULL );
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
	if( threadingInfo == NULL ) return MultiThreadedFunctions::ParallelIntegrate_pthreads( thisFunction, thesePoints, thisBoundary, 4 );
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
	if( threadingInfo == NULL ) return MultiThreadedFunctions::ParallelIntegrate_pthreads( thisFunction, thesePoints, thisBoundary, 4 );
	if( threadingInfo->MultiThreadingInstance == "pthreads" )
	{
		return MultiThreadedFunctions::ParallelIntegrate_pthreads( thisFunction, thesePoints, thisBoundary, threadingInfo->numThreads );
	}
	else
	{
		return new vector<double>();
	}
}

vector<IPDF*> MultiThreadedFunctions::StoredFunctions = vector<IPDF*>();

vector<IPDF*> MultiThreadedFunctions::GetFunctions( IPDF* thisFunction, unsigned int nThreads )
{
	if( !StoredFunctions.empty() )
	{
		if( thisFunction->GetLabel() == StoredFunctions[0]->GetLabel() && StoredFunctions.size() == nThreads )
		{
			for( unsigned int i=0; i< nThreads; ++i )
			{
				StoredFunctions[i]->UpdatePhysicsParameters( thisFunction->GetPhysicsParameters() );
			}
			return StoredFunctions;
		}

		while( !StoredFunctions.empty() )
		{
			if( StoredFunctions.back() != NULL ) delete StoredFunctions.back();
			StoredFunctions.pop_back();
		}
	}

	for( unsigned int i=0; i < nThreads; ++i )
	{
		StoredFunctions.push_back( ClassLookUp::CopyPDF( thisFunction ) );
	}

	return StoredFunctions;
}

vector<double>* MultiThreadedFunctions::ParallelEvaluate_pthreads( IPDF* thisFunction, IDataSet* thesePoints, unsigned int nThreads, ComponentRef* thisComponent )
{
	vector<IDataSet*> payLoad;
	vector<vector<DataPoint*> > datasets_data = Threading::divideData( thesePoints, (int)nThreads );


	for( unsigned int i=0; i < nThreads; ++i )
	{
		payLoad.push_back( (IDataSet*) new MemoryDataSet( thesePoints->GetBoundary(), datasets_data[i] ) );
	}


	vector<IPDF*> functions = MultiThreadedFunctions::GetFunctions( thisFunction, nThreads );

	vector<double>* returnable = MultiThreadedFunctions::ParallelEvaluate_pthreads( functions, payLoad, nThreads, thisComponent );


	while( !payLoad.empty() )
	{
		if( payLoad.back() != NULL ) delete payLoad.back();

		payLoad.pop_back();
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
		cerr << "nThreads: " << nThreads << endl;
		cerr << "thisFunction.size(): " << thisFunction.size() << endl;
		cerr << "thesePoints.size(): " << thesePoints.size() << endl;
		exit(87356);
	}

	//      1 thread per core
	pthread_t* local_Thread = MultiThreadedFunctions::GetFittingThreads( nThreads );//new pthread_t[ nThreads ];
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
		//thesePoints[i]->GetDataPoint(0)->Print();
		fit_thread_data[i].dataPoint_Result.clear();
		if( thisComponent != NULL )
		{
			if( fit_thread_data[i].thisComponent != NULL )
			{
				//if( fit_thread_data[i].thisComponent->getComponentNumber() == thisComponent->getComponentNumber() )
				//{
					delete fit_thread_data[i].thisComponent;
					fit_thread_data[i].thisComponent = new ComponentRef( *thisComponent );
				//}
			}
			else fit_thread_data[i].thisComponent = new ComponentRef( *thisComponent );
		}
		else fit_thread_data[i].thisComponent = NULL;
	}

	//cout << "Creating Threads" << endl;

	//      Create the Threads and set them to be joinable
	for( unsigned int threadnum=0; threadnum< nThreads ; ++threadnum )
	{
		int status=0;
		if( thisComponent == NULL ) status = pthread_create( &local_Thread[threadnum], &attrib, MultiThreadedFunctions::Evaluate_pthread, (void *) &(fit_thread_data[threadnum]) );
		else status = pthread_create( &local_Thread[threadnum], &attrib, MultiThreadedFunctions::EvaluateComponent_pthread, (void *) &(fit_thread_data[threadnum]) );

		//cout << "status: " << status << endl;
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
		int status = pthread_join( local_Thread[threadnum], NULL);
		if( status )
		{
			cerr << "Error Joining a Thread:\t" << threadnum << "\t:\t" << status << "\t...Exiting\n" << endl;
		}
	}

	//      Do some cleaning Up
	pthread_attr_destroy(&attrib);

	unsigned int size=0;
	for( unsigned int i=0; i< thesePoints.size(); ++i ) size+=(unsigned int)thesePoints[i]->GetDataNumber();

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

	//for( unsigned int i=0; i < nThreads; ++i )
	//{
	//	if( fit_thread_data[i].thisComponent != NULL ) delete fit_thread_data[i].thisComponent;
	//}

	//delete[] fit_thread_data;
	//delete[] local_Thread;

	return final_output;
}

void* MultiThreadedFunctions::Evaluate_pthread( void *input_data )
{
	struct Fitting_Thread *thread_input = (struct Fitting_Thread*) input_data;

	double value=0;
	IDataSet* myDataSet = thread_input->dataSet;
	unsigned int number = (unsigned) myDataSet->GetDataNumber();
	for( unsigned int i=0; i < number; ++i )
	{
		//pthread_mutex_t* debug_lock = thread_input->fittingPDF->DebugMutex();
		//cout << i << endl;
		//myDataSet->GetDataPoint( i )->Print();
		try
		{
			value = thread_input->fittingPDF->Evaluate( myDataSet->GetDataPoint( (int)i ) );
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
	IDataSet* myDataSet = thread_input->dataSet;
	unsigned int number = (unsigned) myDataSet->GetDataNumber();
	for( unsigned int i=0; i < number; ++i )
	{
		//pthread_mutex_t* debug_lock = thread_input->fittingPDF->DebugMutex();
		try
		{
			value = thread_input->fittingPDF->EvaluateComponent( myDataSet->GetDataPoint( (int)i ), thread_input->thisComponent );
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
	vector<vector<DataPoint*> > datasets_data = Threading::divideData( thesePoints, (int)nThreads );

	for( unsigned int i=0; i < nThreads; ++i )
	{
		payLoad.push_back( (IDataSet*) new MemoryDataSet( thesePoints->GetBoundary(), datasets_data[i] ) );
	}

	double tempVal = thisFunction->Integral( payLoad[0]->GetDataPoint(0), thisBoundary ); (void) tempVal;

	vector<IPDF*> functions = MultiThreadedFunctions::GetFunctions( thisFunction, nThreads );

	vector<PhaseSpaceBoundary*> boundaries;
	for( unsigned int i=0; i< nThreads; ++i )
	{
		boundaries.push_back( new PhaseSpaceBoundary( *thisBoundary ) );
	}

	vector<double>* returnable = MultiThreadedFunctions::ParallelIntegrate_pthreads( functions, payLoad, boundaries, nThreads );

	//for( unsigned int i=0; i< payLoad.size(); ++i ) if( payLoad[i] != NULL ) delete payLoad[i];
	//for( unsigned int i=0; i< functions.size(); ++i ) if( functions[i] != NULL ) delete functions[i];
	for( unsigned int i=0; i< boundaries.size(); ++i ) if( boundaries[i] != NULL ) delete boundaries[i];

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
			value = thread_input->fittingPDF->Integral( thread_input->dataSet->GetDataPoint((int)i), thread_input->FitBoundary );
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
		cerr << "nThreads: " << nThreads << endl;
		cerr << "thisFunction.size(): " << thisFunction.size() << endl;
		cerr << "thesePoints.size(): " << thesePoints.size() << endl;
		cerr << "theseBoundarys.size(): " << theseBoundarys.size() << endl;
		exit(87356);
	}

	for( unsigned int i=0; i< thisFunction.size(); ++i )
	{
		double tempVal = thisFunction[i]->Integral( thesePoints[i]->GetDataPoint(0), theseBoundarys[i] );
		(void) tempVal;
	}

	//      1 thread per core
	pthread_t* local_Thread = MultiThreadedFunctions::GetFittingThreads( nThreads );//new pthread_t[ nThreads ];
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
		int status = pthread_create( &local_Thread[threadnum], &attrib, MultiThreadedFunctions::Integrate_pthread, (void *) &(fit_thread_data[threadnum]) );
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
		int status = pthread_join( local_Thread[threadnum], NULL);
		if( status )
		{
			cerr << "Error Joining a Thread:\t" << threadnum << "\t:\t" << status << "\t...Exiting\n" << endl;
		}
	}

	//      Do some cleaning Up
	pthread_attr_destroy(&attrib);

	unsigned int size=0;
	for( unsigned int i=0; i< thesePoints.size(); ++i ) size+=(unsigned int)thesePoints[i]->GetDataNumber();

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
	//delete[] local_Thread;
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

