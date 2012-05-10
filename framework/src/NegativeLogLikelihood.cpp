/**
        @class NegativeLogLikelihood

        A fit function with evaulate methods for an NLL calculation

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "NegativeLogLikelihood.h"
//	System Headers
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>

//Default constructor
NegativeLogLikelihood::NegativeLogLikelihood()
{
}

//Destructor
NegativeLogLikelihood::~NegativeLogLikelihood()
{
}

//Return the negative log likelihood for a PDF/DataSet result
double NegativeLogLikelihood::EvaluateDataSet( IPDF * TestPDF, IDataSet * TestDataSet, RapidFitIntegrator * ResultIntegrator, int number )
{
	(void)number;
	//Initialise the integral caching
	//ResultIntegrator->UpdateIntegralCache( TestDataSet->GetBoundary() );

	//Loop over all data points
	double total = 0.0;
	double integral = 0.0;
	double weight = 1.0;
	double value = 0.0;
	DataPoint* temporaryDataPoint=NULL;
	bool flag = false;

	for (int dataIndex = 0; dataIndex < TestDataSet->GetDataNumber(); ++dataIndex)
	{
		temporaryDataPoint = TestDataSet->GetDataPoint(dataIndex);
		value = TestPDF->Evaluate(temporaryDataPoint);

		//Idiot check
		//if ( value < 0 || isnan(value) )
		//{
			//cerr << "PDF evaluates to " << value << endl;
			//	Quickest way to train the fitter not to go wherever this was in phase-space
			//if( TestPDF->GetMCCacheStatus() )	TestPDF->SetMCCacheStatus( false );
			//total-=1;
		//}
		flag = ( (value < 0) || isnan(value) );
		
		//Find out the integral
		integral = ResultIntegrator->Integral( temporaryDataPoint, TestDataSet->GetBoundary() );
		
		//Get the weight for this DataPoint (event)
		weight = 1.0;
		if (useWeights)
		{
			weight = temporaryDataPoint->GetObservable(weightObservableName)->GetValue();
		}

		if( weightsSquared ) weight*=weight;
		total += weight * log( value / integral );
	}

	if( flag ) cerr << "PDF evaluates to " << value << endl;

	//Return negative log likelihood
	return -/*1.0 **/ total;
}

//Return the up value for error calculations
double NegativeLogLikelihood::UpErrorValue( int Sigma )
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
