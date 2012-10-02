/*!
 * @class FitResult
 *
 * @brief Container for all results from a minimisation
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

///	RapidFit Headers
#include "FitResult.h"
#include "ResultParameterSet.h"
#include "RapidFitMatrix.h"
///	System Headers
#include <vector>

using namespace::std;

/*!
 * Constructor with correct arguments, including covariance Matrix and contours
 */
FitResult::FitResult( double MinimumValue, ResultParameterSet * FittedParameters, int FitStatus, PhysicsBottle* FittedBottle,
		RapidFitMatrix* CovarianceMatrix, vector< FunctionContour* > ContourPlots ) :
	minimumValue( MinimumValue ), fittedParameters( new ResultParameterSet(*FittedParameters) ), covarianceMatrix( CovarianceMatrix ),
	contours( ContourPlots ), fitStatus( FitStatus ), fittedBottle( new PhysicsBottle(*FittedBottle) )
{
}

FitResult::FitResult( const FitResult& input ) :
	minimumValue(input.minimumValue), fittedParameters( new ResultParameterSet(*input.fittedParameters) ),
	covarianceMatrix(input.covarianceMatrix==NULL?NULL:new RapidFitMatrix(*input.covarianceMatrix)), contours(input.contours),
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

RapidFitMatrix* FitResult::GetCovarianceMatrix()
{
	return covarianceMatrix;
}

void FitResult::ApplyCovarianceMatrix( RapidFitMatrix* Input )
{
	cout << "Applying FitResult:  ";
	for( unsigned int i=0; i< (unsigned)Input->theseParameters.size(); ++i ) cout << Input->theseParameters[i] << "\t";
	cout << endl;
	if( covarianceMatrix != NULL ) delete covarianceMatrix;
	covarianceMatrix = Input;
	fittedParameters->ApplyCovarianceMatrix( Input );
}

double FitResult::GetMinimumValue()
{
	return minimumValue;
}

void FitResult::SetMinimumValue( double min )
{
	minimumValue = min;
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

