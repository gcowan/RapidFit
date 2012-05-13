/**
	@class AcceptReject

        Class for generating toy data from a PDF.
	Can inherti from this to implement preselection for a particular PDF.

        @author Benjamin M Wynne bwynne@cern.ch
        @date 2009-10-02
*/

//	RapidFit Headers
#include "AcceptReject.h"
#include "PhaseSpaceBoundary.h"
//	System Headers
#include <iostream>
#include <math.h>
#include <float.h>

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

//Constructor with correct argument
AcceptReject::AcceptReject( PhaseSpaceBoundary * NewBoundary, IPDF * NewPDF ) : generationFunction(NewPDF),
	generationBoundary(NewBoundary), dataNumber(0), newDataSet(), rootRandom(), moreThanMaximum(0.01), numberAttempts(0)
{
	newDataSet = new MemoryDataSet(generationBoundary);
	rootRandom = NewPDF->GetRandomFunction();
}

//Destructor
AcceptReject::~AcceptReject()
{
	delete newDataSet;
	//delete rootRandom;
	//delete generationFunction;
	//delete generationBoundary;
}

//Use accept/reject method to create data
int AcceptReject::GenerateData( int DataAmount )
{
	int numberAccepted = 0;
	vector<string> allNames = generationBoundary->GetAllNames();
	vector<string>::iterator nameIterator;

	//Keep trying until required amount of data is generated
	while (numberAccepted < DataAmount)
	{
		++numberAttempts;

		//Create a new point in N-space
		DataPoint * testDataPoint = new DataPoint(allNames);
		for ( nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
		{
			//Create a new observable value
			Observable * newObservable = generationBoundary->GetConstraint(*nameIterator)->CreateObservable(rootRandom);
			testDataPoint->SetObservable( *nameIterator, newObservable );
			delete newObservable;
		}

		//Apply preselection of test values
		double testValue = moreThanMaximum * rootRandom->Rndm();
		if ( Preselection( testDataPoint, testValue ) )
		{
			//Accept/reject
			double functionValue = generationFunction->Evaluate(testDataPoint);
			if ( fabs(functionValue - 0.0) < DOUBLE_TOLERANCE )
			{
				//Will get stuck in infinite loop
				cerr << "Function value zero" << endl;
				delete testDataPoint;
				return dataNumber;
			}

			if (functionValue > moreThanMaximum)
			{
				//Restart accept/reject
				cout << "Function value " << functionValue << " is more than expected maximum " << moreThanMaximum << ": restarting data generation" << endl;
				dataNumber = 0;
				numberAccepted = 0;
				newDataSet->Clear();
				moreThanMaximum *= 2.0;
				delete testDataPoint;
				return GenerateData(DataAmount);
			}
			else
			{
				if (testValue < functionValue)
				{
					//Accept
					newDataSet->AddDataPoint(testDataPoint);
					++numberAccepted;
					delete testDataPoint;
				}
				else
				{
					//Reject
					delete testDataPoint;
				}
			}
		}
		else
		{
			//Reject
			delete testDataPoint;
		}
	}

	//Return data generation statistics
	dataNumber += numberAccepted;
	cout << "Data generation: " << numberAccepted << " accepted from " << numberAttempts << endl;
	numberAttempts = 0;
	return dataNumber;
}

//Return data set
IDataSet * AcceptReject::GetDataSet() const
{
	return newDataSet;
}

//Overload in child functions to speed data generation for complex functions
bool AcceptReject::Preselection( DataPoint * TestDataPoint, double TestValue )
{
	(void)TestDataPoint;
	(void)TestValue;
	return true;
}
