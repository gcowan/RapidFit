/**
        @class FitResult

        Container for all results from a minimisation

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#include "FitResult.h"

//Default constructor
FitResult::FitResult()
{
}

//Constructor with correct arguments
FitResult::FitResult( double MinimumValue, ResultParameterSet * FittedParameters, int FitStatus, PhysicsBottle FittedBottle ) : minimumValue( MinimumValue ),
	fittedParameters( FittedParameters ), fitStatus( FitStatus ), fittedBottle( FittedBottle )
{
}

//Constructor with correct arguments, including covariance Matrix
FitResult::FitResult( double MinimumValue, ResultParameterSet * FittedParameters, int FitStatus, PhysicsBottle FittedBottle, vector<double> covMatrix ) :
	minimumValue( MinimumValue ), fittedParameters( FittedParameters ), fitStatus( FitStatus ), fittedBottle( FittedBottle )
	, covarianceMatrix( covMatrix )
{
}

//Destructor
FitResult::~FitResult()
{
}

vector<double> FitResult::GetCovarianceMatrix()
{
	return covarianceMatrix;
}

double FitResult::GetMinimumValue()
{
	return minimumValue;
}

ResultParameterSet * FitResult::GetResultParameterSet()
{
	return fittedParameters;
}

int FitResult::GetFitStatus()
{
	return fitStatus;
}

PhysicsBottle * FitResult::GetPhysicsBottle()
{
	return &fittedBottle;
}
