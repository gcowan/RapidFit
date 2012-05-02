/**
        @class FitResult

        Container for all results from a minimisation

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "FitResult.h"

//Constructor with correct arguments
FitResult::FitResult( double MinimumValue, ResultParameterSet * FittedParameters, int FitStatus, PhysicsBottle* FittedBottle ) : minimumValue( MinimumValue ), fittedParameters( new ResultParameterSet(*FittedParameters) ), covarianceMatrix(), contours(), fitStatus( FitStatus ), fittedBottle( new PhysicsBottle(*FittedBottle) )
{
}

//Constructor with correct arguments, including covariance Matrix
FitResult::FitResult( double MinimumValue, ResultParameterSet * FittedParameters, int FitStatus, PhysicsBottle* FittedBottle, vector<double> covMatrix ) : minimumValue( MinimumValue ), fittedParameters( new ResultParameterSet(*FittedParameters) ), covarianceMatrix( covMatrix ), contours(), fitStatus( FitStatus ), fittedBottle( new PhysicsBottle(*FittedBottle) )
{
}

//Constructor with correct arguments, including covariance Matrix and contours
FitResult::FitResult( double MinimumValue, ResultParameterSet * FittedParameters, int FitStatus, PhysicsBottle* FittedBottle, vector<double> covMatrix, vector< FunctionContour* > ContourPlots ) : minimumValue( MinimumValue ), fittedParameters( new ResultParameterSet(*FittedParameters) ), covarianceMatrix( covMatrix ), contours( ContourPlots ), fitStatus( FitStatus ), fittedBottle( new PhysicsBottle(*FittedBottle) )
{
}

//FitResult::FitResult( double MinimumValue, ResultParameterSet* FittedParameters, int FitStatus ) : minimumValue( MinimumValue ), fittedParameters( new ResultParameterSet(*FittedParameters) ), covarianceMatrix(), contours(), fitStatus( FitStatus ), fittedBottle(NULL)
//{
//}

FitResult::FitResult( const FitResult& input ) :
minimumValue(input.minimumValue), fittedParameters( new ResultParameterSet(*input.fittedParameters) ), covarianceMatrix(input.covarianceMatrix), contours(input.contours), fitStatus(input.fitStatus), fittedBottle(new PhysicsBottle(*input.fittedBottle))
{
}

//Destructor
FitResult::~FitResult()
{
	if( fittedBottle != NULL )
		delete fittedBottle;
	if( fittedParameters != NULL )
		delete fittedParameters;
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

void FitResult::ForceFitStatus( int input_status )
{
	fitStatus = input_status;
}

PhysicsBottle* FitResult::GetPhysicsBottle()
{
	return fittedBottle;
}

void FitResult::Print() const
{
	fittedParameters->Print();
}

