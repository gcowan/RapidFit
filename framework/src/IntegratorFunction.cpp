/**
  @class IntegratorFunction

  A wrapper to make IPDF usable by the numerical integrator

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-8
  */

#include "StringProcessing.h"
#include "IntegratorFunction.h"
#include <iostream>

using namespace std;

//Default constructor
IntegratorFunction::IntegratorFunction()
{
}

//Constructor with correct argument
IntegratorFunction::IntegratorFunction( IPDF * InputFunction, DataPoint * InputPoint, vector<string> IntegrateThese, vector<string> DontIntegrateThese )
	: wrappedFunction(InputFunction), currentPoint(InputPoint), doIntegrate(IntegrateThese), dontIntegrate(DontIntegrateThese)
{
}

//Constructor with additional information needed for the Foam coordinate transform
IntegratorFunction::IntegratorFunction( IPDF * InputFunction, DataPoint * InputPoint, vector<string> IntegrateThese, vector<string> DontIntegrateThese, vector<double> InputMinima, vector<double> InputRanges )
	: wrappedFunction(InputFunction), currentPoint(InputPoint), doIntegrate(IntegrateThese), dontIntegrate(DontIntegrateThese), minima(InputMinima), ranges(InputRanges)
{
}

//Destructor
IntegratorFunction::~IntegratorFunction()
{
}

//Return the IPDF inside the wrapper
IPDF * IntegratorFunction::GetWrappedFunction()
{
	return wrappedFunction;
}

//Copy the wrapper
IntegratorFunction * IntegratorFunction::Clone() const
{
	return new IntegratorFunction( wrappedFunction, currentPoint, doIntegrate, dontIntegrate );
}

//Return the number of dimensions
unsigned int IntegratorFunction::NDim() const
{
	return int(doIntegrate.size());
}

//Return the function value at x
double IntegratorFunction::operator()( const double * x ) const
{
	return DoEval(x);
}
double IntegratorFunction::operator()( double x ) const
{
	return DoEval(x);
}

//Assignment operator
IBaseFunctionMultiDim & IntegratorFunction::operator=( const IntegratorFunction & NewFunction )
{
	wrappedFunction = NewFunction.wrappedFunction;
	return *this;
}

//Return the function value at x
double IntegratorFunction::DoEval( const double * x ) const
{
	//Make a new data point
	DataPoint * newDataPoint = new DataPoint( currentPoint->GetAllNames() );//WrappedFunction->GetPrototypeDataPoint() );

	//Load the array into the data point
	for (unsigned int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
	{
		Observable * currentObservable = currentPoint->GetObservable( doIntegrate[observableIndex] );
		newDataPoint->SetObservable( doIntegrate[observableIndex], x[observableIndex], currentObservable->GetError(), currentObservable->GetUnit() );
	}

	//Load values of other observables
	for (unsigned int observableIndex = 0; observableIndex < dontIntegrate.size(); observableIndex++ )
	{
		Observable * currentObservable = currentPoint->GetObservable( dontIntegrate[observableIndex] );
		newDataPoint->SetObservable( dontIntegrate[observableIndex], currentObservable->GetValue(), currentObservable->GetError(), currentObservable->GetUnit() );
	}

	//Evaluate
	double result = wrappedFunction->Evaluate(newDataPoint);
	delete newDataPoint;
	return result;
}
double IntegratorFunction::DoEval( double x ) const
{
	if ( doIntegrate.size() == 1 )
	{
		double xArray[1];
		xArray[0] = x;
		return DoEval(xArray);
	}
	else
	{
		cerr << "One dimensional evaluation has been called for multidimensional function" << endl;
		return 0.0;
	}
}
Double_t IntegratorFunction::Density( Int_t ndim, Double_t * xArray )
{
	if ( ndim == int(doIntegrate.size()) )
	{
		//Coordinate transform
		double* transformedArray = new double[ndim];
		for ( int observableIndex = 0; observableIndex < ndim; observableIndex++ )
		{
			transformedArray[observableIndex] = minima[observableIndex] + ( ranges[observableIndex] * xArray[observableIndex] );
		}

		return DoEval(transformedArray);
	}
	else
	{
		cerr << "TFoamIntegrand problem - dimension number mismatch" << endl;
		return 0.0;
	}
}
