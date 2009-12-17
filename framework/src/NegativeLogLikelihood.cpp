/**
        @class NegativeLogLikelihood

        A fit function with evaulate methods for an NLL calculation

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#include "NegativeLogLikelihood.h"
#include <stdlib.h>
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
double NegativeLogLikelihood::EvaluateDataSet( IPDF * TestPDF, IDataSet * TestDataSet, RapidFitIntegrator * ResultIntegrator )
{
	//Loop over all data points
	double total = 0.0;
	for (int dataIndex = 0; dataIndex < TestDataSet->GetDataNumber(); dataIndex++)
	{
		DataPoint * temporaryDataPoint = TestDataSet->GetDataPoint(dataIndex);
		double value = TestPDF->Evaluate(temporaryDataPoint);

		//Idiot check
		if ( value < 0 || isnan(value) )
		{
			cerr << "PDF evaluates to " << value << endl;
		}
		
		//Find out the integral
		double integral = ResultIntegrator->Integral( temporaryDataPoint, TestDataSet->GetBoundary() );
		
		//Get the weight for this DataPoint (event)
		double weight = 1.0;	
		if (useWeights)
		{
			weight = temporaryDataPoint->GetObservable(weightObservableName)->GetValue();
		}
		total += weight * log( value / integral );
	}
	
	//Return negative log likelihood
	return -1.0 * total;
}

//Return the negative log likelihood
double NegativeLogLikelihood::EvaluateParameterSet( ParameterSet * TestParameterSet, vector<string> InterestingParameters )
{
	//Loop over all parameters
	double total = 0.0;
	for (nameIterator = InterestingParameters.begin(); nameIterator != InterestingParameters.end(); nameIterator++)
	{
		PhysicsParameter * testParameter = TestParameterSet->GetPhysicsParameter( *nameIterator );

		double difference = testParameter->GetValue() - testParameter->GetOriginalValue();
		double error = testParameter->GetMaximum() - testParameter->GetMinimum();
		total += log( difference * difference / error );
	}

	//Return negative log likelihood
	return -1.0 * total;
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
	else
	{
		cerr << "I don't know UP for NLL sigma > 2" << endl;
		exit(1);
	}
}
