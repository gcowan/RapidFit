/*!
 * @class FitResult
 *
 * @brief Container for all results from a minimisation
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

///	ROOT Headers
#include "TMatrixDSym.h"
///	RapidFit Headers
#include "FitResult.h"
#include "ResultParameterSet.h"
///	System Headers
#include <vector>

using namespace::std;

/*!
 * Constructor with correct arguments, including covariance Matrix and contours
 */
FitResult::FitResult( double MinimumValue, ResultParameterSet * FittedParameters, int FitStatus, PhysicsBottle* FittedBottle,
		TMatrixDSym* covarianceMatrix, vector< FunctionContour* > ContourPlots ) :
	minimumValue( MinimumValue ), fittedParameters( new ResultParameterSet(*FittedParameters) ), covarianceMatrix( covarianceMatrix ),
	contours( ContourPlots ), fitStatus( FitStatus ), fittedBottle( new PhysicsBottle(*FittedBottle) )
{
}

FitResult::FitResult( const FitResult& input ) :
	minimumValue(input.minimumValue), fittedParameters( new ResultParameterSet(*input.fittedParameters) ),
	covarianceMatrix(input.covarianceMatrix==NULL?NULL:new TMatrixDSym(*input.covarianceMatrix)), contours(input.contours),
	fitStatus(input.fitStatus), fittedBottle(new PhysicsBottle(*input.fittedBottle))
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

TMatrixDSym* FitResult::GetCovarianceMatrix()
{
	return covarianceMatrix;
}

void FitResult::ApplyCovarianceMatrix( TMatrixDSym* Input )
{
	if( covarianceMatrix != NULL ) delete covarianceMatrix;
	covarianceMatrix = Input;
	fittedParameters->ApplyCovarianceMatrix( Input );
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

