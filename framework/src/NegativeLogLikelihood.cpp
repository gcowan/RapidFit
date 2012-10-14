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
NegativeLogLikelihood::NegativeLogLikelihood() : FitFunction()
{
	Name="NegativeLogLikelihood";
}

//Destructor
NegativeLogLikelihood::~NegativeLogLikelihood()
{
}

//Return the negative log likelihood for a PDF/DataSet result
double NegativeLogLikelihood::EvaluateDataSet( IPDF * TestPDF, IDataSet * TestDataSet, int number )
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
	//bool flag = false;

	for (int dataIndex = 0; dataIndex < TestDataSet->GetDataNumber(); ++dataIndex)
	{
		temporaryDataPoint = TestDataSet->GetDataPoint(dataIndex);
		value = TestPDF->Evaluate(temporaryDataPoint);

		if( debug != NULL )	if( debug->DebugThisClass( "NegativeLogLikelihood" ) ) cout << "V: " << value << endl;
		//Idiot check
		//if ( value < 0 || isnan(value) )
		//{
			//cerr << "PDF evaluates to " << value << endl;
			//	Quickest way to train the fitter not to go wherever this was in phase-space
			//if( TestPDF->GetMCCacheStatus() )	TestPDF->SetMCCacheStatus( false );
			//total-=1;
		//}
		//flag = ( (value < 0) || isnan(value) );
		
		//Find out the integral
		integral = TestPDF->Integral( temporaryDataPoint, TestDataSet->GetBoundary() );
		
		if( debug != NULL )	if( debug->DebugThisClass( "NegativeLogLikelihood" ) ) cout << "I: " << integral << endl;

		//Get the weight for this DataPoint (event)
		weight = 1.0;
		if( useWeights )
		{
			weight = temporaryDataPoint->GetEventWeight();
		}

		double pointValue= log( value / integral );

		if( useWeights ) pointValue *= weight;
		if( useWeights && weightsSquared ) pointValue *= weight;

		total+=pointValue;

		//cout << total << " " << value << " " << integral << endl;
	}

	if( false ) cerr << "PDF evaluates to " << value << endl;

	//Return negative log likelihood
	return -total;
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
