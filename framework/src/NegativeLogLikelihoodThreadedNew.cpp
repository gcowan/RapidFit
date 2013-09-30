/**
  @class NegativeLogLikelihoodThreadedNew

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
#include "NegativeLogLikelihoodThreadedNew.h"
#include "MultiThreadedFunctions.h"
#include "ThreadingConfig.h"
#include "IPDF.h"
#include "MemoryDataSet.h"
//	System Headers
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <pthread.h>
#include <float.h>

using namespace::std;

//Default constructor
NegativeLogLikelihoodThreadedNew::NegativeLogLikelihoodThreadedNew() : FitFunction()
{
	Name="NegativeLogLikelihoodThreadedNew";
}

//Destructor
NegativeLogLikelihoodThreadedNew::~NegativeLogLikelihoodThreadedNew()
{
}

//Return the negative log likelihood for a PDF/DataSet result
double NegativeLogLikelihoodThreadedNew::EvaluateDataSet( IPDF * FittingPDF, IDataSet * TotalDataSet, int number )
{
	if( TotalDataSet->GetDataNumber() == 0 ) return 0.;

	if( Threads <= 0 )
	{
		cerr<< "Bad Number of Threads: " << Threads << " check your XML!!!" << endl << endl;
		exit(-125);
	}

	ThreadingConfig* thisConfig = new ThreadingConfig();

	thisConfig->MultiThreadingInstance = "pthreads";
	thisConfig->numThreads=(unsigned)Threads;
	thisConfig->wantedComponent = NULL;

	//cout << "Breaking into Threads" << endl;
	vector<double>* values = MultiThreadedFunctions::ParallelEvaulate( FittingPDF, TotalDataSet, thisConfig ); 
	//vector<double>* values = MultiThreadedFunctions::ParallelEvaulate( stored_pdfs, stored_datasets[number], thisConfig );

	//cout << "Continuing Threads" << endl;
	vector<double>* integrals = MultiThreadedFunctions::ParallelIntegrate( FittingPDF, TotalDataSet, TotalDataSet->GetBoundary(), thisConfig );
	//vector<double>* integrals = MultiThreadedFunctions::ParallelIntegrate( stored_pdfs, stored_datasets[number], StoredBoundary, thisConfig );


	//cout << "Collecting Results" << endl;
	//cout << TotalDataSet->GetDataNumber() << "\t\t" << values->size() << "\t" << integrals->size() << endl;

	double total=0;

	if( TotalDataSet->GetWeightsWereUsed() )
	{
		for( unsigned int i=0; i< values->size(); ++i )
		{
			double weight = TotalDataSet->GetDataPoint( i )->GetEventWeight();
			total += log( (*values)[i] / (*integrals)[i] ) * weight;
		}
	}
	else
	{
                for( unsigned int i=0; i< values->size(); ++i )
                {
                        total += log( (*values)[i] / (*integrals)[i] );
                }
	}

	delete values;
	delete integrals;
	delete thisConfig;

	//cout << total << endl;
	//exit(0);
	return -total;
}

//Return the up value for error calculations
double NegativeLogLikelihoodThreadedNew::UpErrorValue( int Sigma )
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

