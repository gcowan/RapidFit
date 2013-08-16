/**
  @class NegativeLogLikelihoodNumerical

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
#include "NegativeLogLikelihoodNumerical.h"
#include "ClassLookUp.h"
//	System Headers
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <pthread.h>
#include <float.h>
#include "TFile.h"
#include "TNtuple.h"

using namespace::std;

pthread_mutex_t _n_eval_lock;
pthread_mutex_t _n_int_lock;

//Default constructor
NegativeLogLikelihoodNumerical::NegativeLogLikelihoodNumerical() : FitFunction()
{
	Name="NegativeLogLikelihoodNumerical";
}

//Destructor
NegativeLogLikelihoodNumerical::~NegativeLogLikelihoodNumerical()
{
}

//Return the negative log likelihood for a PDF/DataSet result
double NegativeLogLikelihoodNumerical::EvaluateDataSet( IPDF * FittingPDF, IDataSet * TotalDataSet, int number )
{
	(void) FittingPDF;

	//	Have to provide a datapoint even though one is _not_ expected to be explicitly used for this fittting function
	double this_integral = FittingPDF->GetPDFIntegrator()->Integral( TotalDataSet->GetDataPoint(0), TotalDataSet->GetBoundary() );

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


	//	Initialize the Fitting_Thread objects which contain the objects to be passed to each thread
	for( unsigned int threadnum=0; threadnum< (unsigned)Threads; ++threadnum )
	{
		fit_thread_data[threadnum].dataSubSet = StoredDataSubSet[(unsigned)number][threadnum];
		fit_thread_data[threadnum].fittingPDF = stored_pdfs[((unsigned)number)*(unsigned)Threads + threadnum];
		fit_thread_data[threadnum].useWeights = useWeights;					//	Defined in the fitfunction baseclass
		fit_thread_data[threadnum].FitBoundary = StoredBoundary[(unsigned)Threads*((unsigned)number)+threadnum];
		//fit_thread_data[threadnum].ResultIntegrator = StoredIntegrals[(unsigned)Threads*((unsigned)number)+threadnum];
		fit_thread_data[threadnum].stored_integral = this_integral;
		fit_thread_data[threadnum].dataPoint_Result = vector<double>();
		fit_thread_data[threadnum].weightsSquared = weightsSquared;
	}

	//cout << "Create Threads" << endl;

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
		while( !fit_thread_data[threadnum].dataPoint_Result.empty() ) fit_thread_data[threadnum].dataPoint_Result.pop_back();
	}
	//cout << -total << endl;
	//exit(100);
	//delete [] fit_thread_data;
	delete [] Thread;

	//cout << total << endl;
	return -total;
}

void* NegativeLogLikelihoodNumerical::ThreadWork( void *input_data )
{
	struct Fitting_Thread *thread_input = (struct Fitting_Thread*) input_data;

	double value=0, weight=0, integral=0, result=0;
	int num=0;
	//TFile * file = TFile::Open("integral.root", "UPDATE");
	//TNtuple * ntuple = new TNtuple("test_good","test", "value");
	//TNtuple * ntuple = new TNtuple("integral","test", "integral");
	//bool isnorm = thread_input->fittingPDF->GetName()=="NormalisedSum";
	for( vector<DataPoint*>::iterator data_i=thread_input->dataSubSet.begin(); data_i != thread_input->dataSubSet.end(); ++data_i, ++num )
	{
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
			integral = thread_input->stored_integral;// ResultIntegrator->Integral( *data_i, thread_input->FitBoundary );
		}
		catch( ... )
		{
			integral = DBL_MAX;
		}

		if( std::isnan(value) == true )
		{
			thread_input->dataPoint_Result.push_back( DBL_MAX );
			cout << endl << "PDF is nan" << endl;
			(*data_i)->Print();
			break;
		}
		if( std::isnan(integral) == true )
		{
			thread_input->dataPoint_Result.push_back( DBL_MAX );
			cout << endl << "Integral is nan" << endl;
			(*data_i)->Print();
			break;
		}
		if( value <= 0 )
		{
			thread_input->dataPoint_Result.push_back( DBL_MAX );
			cout << endl << "Value is <=0 " << value << endl;
			(*data_i)->Print();
			break;
		}
		if( integral <= 0 )
		{
			thread_input->dataPoint_Result.push_back( DBL_MAX );
			cout << endl << "Integral is <= 0 " << integral << endl;
			(*data_i)->Print();
			break;
		}

		if( value >= DBL_MAX || integral >= DBL_MAX )
		{
			thread_input->dataPoint_Result.push_back( DBL_MAX );
			cerr << endl << "Caught invalid value from PDF" << endl;
			(*data_i)->Print();
			break;
		}

		//	Result of evaluating the DataPoint
		result = log( value / integral );
		//result = integral;
		//result = value;
		//cout << value << " " << integral << endl;

		//	If we have a weighted dataset then weight the result (if not don't perform a *1.)
		if( thread_input->useWeights == true )
		{
			weight = (*data_i)->GetEventWeight();
			//pthread_mutex_lock( &_n_eval_lock );
			result *= weight;
			if( thread_input->weightsSquared ) result *= weight;
			//pthread_mutex_unlock( &_n_eval_lock );
		}

		//	Push back the result from evaluating this datapoint
		thread_input->dataPoint_Result.push_back( result );
	}
	//ntuple->Fill(integral);
	//file->Write();
	//file->Close();
	//	Finished evaluating this thread
	pthread_exit( NULL );
}

//Return the up value for error calculations
double NegativeLogLikelihoodNumerical::UpErrorValue( int Sigma )
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

