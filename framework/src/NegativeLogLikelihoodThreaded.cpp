/**
  @class NegativeLogLikelihoodThreaded

  A fit function with evaulate methods for an NLL calculation using mutliple cores (shamelessly copied from the NLL class)

  There are 2 ways of multithreading RapidFit to make it more efficient/faster

  1)
  Make Minuit or your minimisation engine faster and call the fit function (PDF) in thread safe ways
  (GOOD LUCK!!!)

  2)
  Behind the scenes of the minimisation engine evaluate the dataset using multiple threads at each point
  (carries a slight overhead due to thread initialization/cleaning up,
  but this is << than the *n speed improvement that is obtained by running on mutliple cores)

  @author Rob Currie rcurrie@cern.ch
  @date 2011
  */

//	ROOT Headers
#include "RooMath.h"
//	RapidFit Headers
#include "NegativeLogLikelihoodThreaded.h"
#include "ClassLookUp.h"
#include "IPDF.h"
//	System Headers
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <pthread.h>
#include <float.h>

using namespace::std;

pthread_mutex_t eval_lock;

//Default constructor
NegativeLogLikelihoodThreaded::NegativeLogLikelihoodThreaded() : FitFunction()
{
	Name="NegativeLogLikelihoodThreaded";
}

//Destructor
NegativeLogLikelihoodThreaded::~NegativeLogLikelihoodThreaded()
{
}

//Return the negative log likelihood for a PDF/DataSet result
double NegativeLogLikelihoodThreaded::EvaluateDataSet( IPDF * FittingPDF, IDataSet * TotalDataSet, int number )
{
	(void) FittingPDF;

	if( TotalDataSet->GetDataNumber() == 0 ) return 0.;

	if( Threads <= 0 )
	{
		cerr<< "Bad Number of Threads: " << Threads << " check your XML!!!" << endl << endl;
		exit(-125);
	}

	//	1 thread per core
	pthread_t* Thread = new pthread_t[ (unsigned)Threads ];
	pthread_attr_t attrib;

	//	Threads HAVE to be joinable
	//	We CANNOT _AND_SHOULD_NOT_ ***EVER*** return information to Minuit without the results from ALL threads successfully returned
	//	Not all pthread implementations are required to obey this as default when constructing the thread _SO_BE_EXPLICIT_
	pthread_attr_init(&attrib);
	pthread_attr_setdetachstate(&attrib, PTHREAD_CREATE_JOINABLE);

	//cout << "Setup Threads: " << Threads << endl;
	ObservableRef weightObservableRef( weightObservableName );

	/*
	   for( int i=0; i< StoredDataSubSet.size(); ++i )
	   {
	   if( i == number )
	   {
	   int tot=0;
	   for( int j=0; j< StoredDataSubSet[i].size(); ++j )
	   {
	   cout << "+" << StoredDataSubSet[i][j].size() << endl;
	   tot+=StoredDataSubSet[i][j].size();
	   }
	   cout << "=" << tot << endl;
	   }
	   }
	   */

	//	Initialize the Fitting_Thread objects which contain the objects to be passed to each thread
	for( unsigned int threadnum=0; threadnum< (unsigned)Threads; ++threadnum )
	{
		fit_thread_data[threadnum].dataSubSet = StoredDataSubSet[(unsigned)number][threadnum];
		fit_thread_data[threadnum].fittingPDF = stored_pdfs[((unsigned)number)*(unsigned)Threads + threadnum];
		fit_thread_data[threadnum].fittingPDF->SetDebugMutex( &eval_lock, false );
		fit_thread_data[threadnum].useWeights = useWeights;					//	Defined in the fitfunction baseclass
		fit_thread_data[threadnum].FitBoundary = StoredBoundary[(unsigned)Threads*((unsigned)number)+threadnum];
		fit_thread_data[threadnum].dataPoint_Result = vector<double>();
		fit_thread_data[threadnum].weightsSquared = weightsSquared;
	}

	//cout << "Creating Threads" << endl;

	//	Create the Threads and set them to be joinable
	for( unsigned int threadnum=0; threadnum< (unsigned)Threads ; ++threadnum )
	{
		int status = pthread_create(&Thread[threadnum], &attrib, this->ThreadWork, (void *) &fit_thread_data[threadnum] );
		if( status )
		{
			cerr << "ERROR:\tfrom pthread_create()\t" << status << "\t...Exiting\n" << endl;
			exit(-1);
		}
	}

	//cout << "Joining Threads!!" << endl;

	//	Join the Threads
	for( unsigned int threadnum=0; threadnum< (unsigned)Threads ; ++threadnum )
	{
		int status = pthread_join( Thread[threadnum], NULL);
		if( status )
		{
			cerr << "Error Joining a Thread:\t" << threadnum << "\t:\t" << status << "\t...Exiting\n" << endl;
		}
	}

	//      Do some cleaning Up
	pthread_attr_destroy(&attrib);

	//cout << "Leaving Threads" << endl;

	double total=0;

	for( unsigned int threadnum=0; threadnum< (unsigned)Threads; ++threadnum )
	{
		for( unsigned int point_num=0; point_num< fit_thread_data[threadnum].dataPoint_Result.size(); ++point_num )
		{
			if( fit_thread_data[threadnum].dataPoint_Result[ point_num ] >= DBL_MAX )
			{
				return DBL_MAX;
			}
			total+= fit_thread_data[threadnum].dataPoint_Result[ point_num ];
		}
		vector<double> empty;
		fit_thread_data[threadnum].dataPoint_Result.swap( empty );
	}

	//delete [] fit_thread_data;
	delete [] Thread;

	//cout << total << endl;
	//exit(0);

	//cout << total << endl;
	return -total;
}

void* NegativeLogLikelihoodThreaded::ThreadWork( void *input_data )
{
	struct Fitting_Thread *thread_input = (struct Fitting_Thread*) input_data;

	double value=0, weight=0, integral=0, result=0;
	int num=0;
	//bool isnorm = thread_input->fittingPDF->GetName()=="NormalisedSum";
	for( vector<DataPoint*>::iterator data_i=thread_input->dataSubSet.begin(); data_i != thread_input->dataSubSet.end(); ++data_i, ++num )
	{

		pthread_mutex_t* debug_lock = thread_input->fittingPDF->DebugMutex();
		try
		{
			value = thread_input->fittingPDF->Evaluate( *data_i );
		}
		catch( ... )
		{
			value = DBL_MAX;
		}

		try
		{
			integral = thread_input->fittingPDF->Integral( *data_i, thread_input->FitBoundary );
		}
		catch( ... )
		{
			integral = DBL_MAX;
		}

		/*
		   if(num==5) exit(0);
		   pthread_mutex_lock( debug_lock );
		   cout << "hello" << endl;
		   pthread_mutex_unlock( debug_lock );
		   */

		//cout << value << "\t" << integral << endl;

		if( std::isnan(value) == true )
		{
			pthread_mutex_lock( debug_lock );
			thread_input->dataPoint_Result.push_back( DBL_MAX );
			cout << endl << "PDF is nan" << endl;
			(*data_i)->Print();
			pthread_mutex_unlock( debug_lock );
			break;
		}
		if( std::isnan(integral) == true )
		{
			pthread_mutex_lock( debug_lock );
			thread_input->dataPoint_Result.push_back( DBL_MAX );
			cout << endl << "Integral is nan" << endl;
			(*data_i)->Print();
			pthread_mutex_unlock( debug_lock );
			break;
		}
		if( value <= 0 )
		{
			pthread_mutex_lock( debug_lock );
			thread_input->dataPoint_Result.push_back( DBL_MAX );
			cout << endl << "Value is <=0 " << value << endl;
			(*data_i)->Print();
			pthread_mutex_unlock( debug_lock );
			break;
		}
		if( integral <= 0 )
		{
			pthread_mutex_lock( debug_lock );
			thread_input->dataPoint_Result.push_back( DBL_MAX );
			cout << endl << "Integral is <= 0 " << integral << endl;
			(*data_i)->Print();
			pthread_mutex_unlock( debug_lock );
			break;
		}

		if( value >= DBL_MAX || integral >= DBL_MAX )
		{
			pthread_mutex_lock( debug_lock );
			thread_input->dataPoint_Result.push_back( DBL_MAX );
			cerr << endl << "Caught invalid value from PDF: " << endl;
			cerr << "Val: " << value << "\tNorm: " << integral << endl;
			(*data_i)->Print();
			pthread_mutex_unlock( debug_lock );
			break;
		}

		//if( value / integral > 1. )
		//{
		//	cout << "SERIOUSE: " << value / integral << endl;
		//}

		//	Result of evaluating the DataPoint
		result = log( value / integral );

		//	If we have a weighted dataset then weight the result (if not don't perform a *1.)
		if( thread_input->useWeights == true )
		{
			weight = (*data_i)->GetEventWeight();
			//pthread_mutex_lock( &eval_lock );
			result *= weight;
			if( thread_input->weightsSquared )
			{
				result *= weight;
				if( weight < 0 ) result *= -1.;
			}
			//pthread_mutex_unlock( &eval_lock );
		}

		//	Push back the result from evaluating this datapoint
		thread_input->dataPoint_Result.push_back( result );
	}

	//	Finished evaluating this thread
	pthread_exit( NULL );
}

//Return the up value for error calculations
double NegativeLogLikelihoodThreaded::UpErrorValue( int Sigma )
{
	if ( Sigma == 1 )
	{
		return 0.5;
	}
	else if ( Sigma == 2 )
	{
		return 2.0;
	}
	else if ( Sigma == 3 )
	{
		return 4.5;
	}
	else
	{
		cerr << "I don't know UP for NLL sigma > 3" << endl;
		exit(1);
	}
}

