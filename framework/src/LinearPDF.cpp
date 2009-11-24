/**
        @class LinearPDF

        An example PDF implementing IPDF directly, without BasePDF

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#include "LinearPDF.h"
#include <iostream>

//Constructor
LinearPDF::LinearPDF() : interceptName( "cIntercept" ), gradientName( "mGradient" )
{
	MakePrototypes();
}

//Constructor allowing you to specify the parameter names
LinearPDF::LinearPDF(string CInterceptName, string MGradientName) : interceptName( CInterceptName ), gradientName( MGradientName )
{
	MakePrototypes();
}

//Make the data point and parameter set
void LinearPDF::MakePrototypes()
{
	//Make the DataPoint prototype
	observableName = "xValue";
	allObservables.push_back( observableName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( interceptName );
	parameterNames.push_back( gradientName );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
LinearPDF::~LinearPDF()
{
}

//Do the integration
double LinearPDF::Normalisation(PhaseSpaceBoundary * NewBoundary)
{
	double xMinimum = NewBoundary->GetConstraint( observableName )->GetMinimum();
	double xMaximum = NewBoundary->GetConstraint( observableName )->GetMaximum();
	double mGradient = allParameters.GetPhysicsParameter( gradientName )->GetValue();
	double cIntercept = allParameters.GetPhysicsParameter( interceptName )->GetValue();

	double yOne = (mGradient * xMinimum) + cIntercept;
	double yTwo = (mGradient * xMaximum) + cIntercept;
	return (xMaximum - xMinimum) * (yOne + yTwo) / 2.0;
}

//Calculate the function value
double LinearPDF::Evaluate(DataPoint * NewDataPoint)
{
	double xValue = NewDataPoint->GetObservable( observableName )->GetValue();
	double mGradient = allParameters.GetPhysicsParameter( gradientName )->GetValue();
	double cIntercept = allParameters.GetPhysicsParameter( interceptName )->GetValue();

	return  (mGradient * xValue) + cIntercept;
}
