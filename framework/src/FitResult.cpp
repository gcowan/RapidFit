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
FitResult::FitResult( double MinimumValue, ResultParameterSet * FittedParameters, int FitStatus, PhysicsBottle FittedBottle ) : minimumValue( MinimumValue ), fittedParameters( FittedParameters ), fitStatus( FitStatus ), fittedBottle( FittedBottle )
{
}

//Constructor with correct arguments, including covariance Matrix
FitResult::FitResult( double MinimumValue, ResultParameterSet * FittedParameters, int FitStatus, PhysicsBottle FittedBottle, vector<double> covMatrix ) : minimumValue( MinimumValue ), fittedParameters( FittedParameters ), covarianceMatrix( covMatrix ), fitStatus( FitStatus ), fittedBottle( FittedBottle )
{
}

//Constructor with correct arguments, including covariance Matrix and contours
FitResult::FitResult( double MinimumValue, ResultParameterSet * FittedParameters, int FitStatus, PhysicsBottle FittedBottle, vector<double> covMatrix, vector< FunctionContour* > ContourPlots ) : minimumValue( MinimumValue ), fittedParameters( FittedParameters ), covarianceMatrix( covMatrix ), contours( ContourPlots ), fitStatus( FitStatus ), fittedBottle( FittedBottle )
{
}


//Destructor
FitResult::~FitResult()
{
}

vector< FunctionContour* > FitResult::GetContours()
{
	return contours;
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
