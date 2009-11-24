/**
        @class InvalidObject

        A generic return type indicating the requested class could not be instantiated

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#include "InvalidObject.h"
#include <iostream>

InvalidObject::InvalidObject() : errorMessage( "Invalid object created" )
{
}

InvalidObject::InvalidObject( string ErrorMessage ) : errorMessage( ErrorMessage )
{
}

InvalidObject::~InvalidObject()
{
}

//Indicate whether the function has been set up correctly
bool InvalidObject::IsValid()
{
	cerr << errorMessage << endl;
	return false;
}

//Set the function parameters
bool InvalidObject::SetPhysicsParameters( ParameterSet * NewParameters )
{
	cerr << errorMessage << endl;
	return false;
}

//Return the integral of the function over the given boundary
double InvalidObject::Integral( DataPoint* NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	cerr << errorMessage << endl;
	return 1.0;
}

//Return the function value at the given point
double InvalidObject::Evaluate( DataPoint * NewPoint )
{
	cerr << errorMessage << endl;
	return 0.0;
}

//Return a prototype data point
vector<string> InvalidObject::GetPrototypeDataPoint()
{
	cerr << errorMessage << endl;
	vector<string> empty;
	return empty;
}

//Return a prototype set of physics parameters
vector<string> InvalidObject::GetPrototypeParameterSet()
{
	cerr << errorMessage << endl;
	vector<string> empty;
	return empty;
}

//Return a list of parameters not to be integrated
vector<string> InvalidObject::GetDoNotIntegrateList()
{
	cerr << errorMessage << endl;
        vector<string> empty;
        return empty;
}

void InvalidObject::Minimise( FitFunction * NewFunction )
{
	cerr << errorMessage << endl;
}

FitResult * InvalidObject::GetFitResult()
{
	cerr << errorMessage << endl;
	return NULL;
}

double InvalidObject::UpErrorValue()
{
	cerr << errorMessage << endl;
	return 1.0;
}

double InvalidObject::EvaluateDataSet( IPDF * NewPDF, IDataSet * NewDataSet )
{
	cerr << errorMessage << endl;
	return 1.0;
}

double InvalidObject::EvaluateParameterSet( ParameterSet * NewParameterSet )
{
	cerr << errorMessage << endl;
	return 1.0;
}

DataPoint * InvalidObject::GetDataPoint( int DataIndex )
{
	cerr << errorMessage << endl;
	return NULL;
}

void InvalidObject::AddDataPoint( DataPoint * dp )
{
        cerr << errorMessage << endl;
}


int InvalidObject::GetDataNumber()
{
	cerr << errorMessage << endl;
	return 0;
}

PhaseSpaceBoundary * InvalidObject::GetBoundary()
{
	cerr << errorMessage << endl;
	return NULL;
}

int InvalidObject::GenerateData(int)
{
	cerr << errorMessage << endl;
	return 0;
}

IDataSet * InvalidObject::GetDataSet()
{
	cerr << errorMessage << endl;
	return NULL;
}
