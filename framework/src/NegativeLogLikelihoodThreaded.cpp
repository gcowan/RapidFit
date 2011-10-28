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
//	System Headers
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <pthread.h>

pthread_mutex_t eval_lock;

//Default constructor
NegativeLogLikelihoodThreaded::NegativeLogLikelihoodThreaded()
{
}

//Destructor
NegativeLogLikelihoodThreaded::~NegativeLogLikelihoodThreaded()
{
}

//Return the negative log likelihood for a PDF/DataSet result
double NegativeLogLikelihoodThreaded::EvaluateDataSet( IPDF * FittingPDF, IDataSet * TotalDataSet, RapidFitIntegrator * ResultIntegrator )
{
	//Initialise the integral caching
	ResultIntegrator->UpdateIntegralCache( TotalDataSet->GetBoundary() );

	//	Get the number of cores on the compile machine
	unsigned int num_cores = (unsigned)Threading::numCores();

	if( this->GetThreads() > 0 )	{	num_cores = (unsigned)this->GetThreads();	}

	//	1 thread per core
	pthread_t Thread[ num_cores ];
	pthread_attr_t attrib;


	//	Threads HAVE to be joinable
	//	We CANNOT _AND_SHOULD_NOT_ ***EVER*** return information to Minuit without the results from ALL threads successfully returned
	//	Not all pthread implementations are required to obey this as default when constructing the thread _SO_BE_EXPLICIT_
	pthread_attr_init(&attrib);
	pthread_attr_setdetachstate(&attrib, PTHREAD_CREATE_JOINABLE);

	//	Create simple data subsets. We no longer care about the handles that IDataSet takes care of
	vector<vector<DataPoint*> > DataSet = Threading::divideData( TotalDataSet, (int)num_cores );
	struct Fitting_Thread* fit_thread_data = new Fitting_Thread[ num_cores ];

	//	Initialize the Fitting_Thread objects which contain the objects to be passed to each thread
	vector<double> new_results;
	for( unsigned int threadnum=0; threadnum< num_cores; ++threadnum )
	{
		fit_thread_data[threadnum].dataSubSet = DataSet[threadnum];
		fit_thread_data[threadnum].fittingPDF = ClassLookUp::CopyPDF( FittingPDF );
		fit_thread_data[threadnum].useWeights = useWeights;				//	Defined in the fitfunction baseclass
		fit_thread_data[threadnum].weightName = weightObservableName;			//	Defined in the fitfunction baseclass
		fit_thread_data[threadnum].dataPoint_Result = new_results;
		fit_thread_data[threadnum].FitBoundary = TotalDataSet->GetBoundary();

		//	Required to give the integral cache code a sharp kick!
		ResultIntegrator->Integral( DataSet[threadnum][0], TotalDataSet->GetBoundary(), true );

		fit_thread_data[threadnum].ResultIntegrator = new RapidFitIntegrator( *(ResultIntegrator) );
	}

	//	Create the Threads and set them to be joinable
	for( unsigned int threadnum=0; threadnum<num_cores ; ++threadnum )
	{
		int status = pthread_create(&Thread[threadnum], &attrib, this->ThreadWork, (void *) &fit_thread_data[threadnum] ); 
		if (status)
		{
			cerr << "ERROR:\tfrom pthread_create()\t" << status << "\t...Exiting\n" << endl;
			exit(-1);
		}
	}

	//	Do some cleaning Up
	pthread_attr_destroy(&attrib);

	//cout << "Joining Threads!!" << endl;

	//	Join the Threads
	for( unsigned int threadnum=0; threadnum<num_cores ; ++threadnum )
	{
		int status = pthread_join( Thread[threadnum], NULL);
		if( status )
		{
			cerr << "Error Joining a Thread:\t" << threadnum << "\t:\t" << status << "\t...Exiting\n" << endl;
		}
	}

	//cout << "Leaving Threads" << endl;

	double total=0;

	for( unsigned int threadnum=0; threadnum<num_cores; ++threadnum )
	{
		for( unsigned int point_num=0; point_num< fit_thread_data[threadnum].dataPoint_Result.size(); ++point_num )
		{
			total+= fit_thread_data[threadnum].dataPoint_Result[ point_num ];
		}
	}

	//	Cleanup the objects we created specifically for threading
	for( unsigned int i=0; i< num_cores; ++i )
	{
		delete fit_thread_data[i].fittingPDF;
		delete fit_thread_data[i].ResultIntegrator;
	}

	//cout << total << endl;
	return -total;

}

void* NegativeLogLikelihoodThreaded::ThreadWork( void *input_data )
{
	struct Fitting_Thread *thread_input = (struct Fitting_Thread*) input_data;

	double value=0, integral=0, result=0;
	for( unsigned int i=0; i< thread_input->dataSubSet.size(); ++i )
	{
		//	Evaluate the PDF at this point
		value = thread_input->fittingPDF->Evaluate( thread_input->dataSubSet[i] );

		//	Get the PDF Integral for this point (usually not dependent on the datapoint)
		integral = thread_input->ResultIntegrator->Integral( thread_input->dataSubSet[i], thread_input->FitBoundary, false );

		//	Result of evaluating the DataPoint
		result = log( value / integral );

		//	If we have a weighted dataset then weight the result (if not don't perform a *1)
		if( thread_input->useWeights ) { result *= thread_input->dataSubSet[i]->GetObservable(thread_input->weightName)->GetValue(); }

		//	Push back the result from evaluating this datapoint
		thread_input->dataPoint_Result.push_back( result );
	}

	//	Finished evaluating this thread
	pthread_exit( (void*) input_data);
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
		cerr << "I don't know UP for NLL sigma > 2" << endl;
		exit(1);
	}
}
