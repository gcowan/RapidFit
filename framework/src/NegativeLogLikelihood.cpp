/**
        @class NegativeLogLikelihood

        A fit function with evaulate methods for an NLL calculation

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#include "NegativeLogLikelihood.h"
#include <cmath>
//#include "Rtypes.h"
#include <iostream>

//Default constructor
NegativeLogLikelihood::NegativeLogLikelihood() : up(0.5)
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
		
		//Find out the integral
		//double integral = TestPDF->Integral( temporaryDataPoint, TestDataSet->GetBoundary() );
		double integral = ResultIntegrator->Integral( temporaryDataPoint, TestDataSet->GetBoundary() );
		
		//Get the weight for this DataPoint (event)
		Observable * weightObs = temporaryDataPoint->GetObservable("eventWeight");
		double weight = 1.;	
		if ( weightObs->GetUnit() != "NameNotFoundError")
		{
			double weight = weightObs->GetValue();
		}
		total += weight * log( value / integral );
		
		//delete temporaryDataPoint;
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

//Set the up value for error calculations
void NegativeLogLikelihood::SetUpErrorValue(double value)
{
	up = value;
}

//Return the up value for error calculations
double NegativeLogLikelihood::UpErrorValue()
{
	return up;
}

