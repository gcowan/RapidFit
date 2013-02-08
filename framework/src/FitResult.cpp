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
FitResult::FitResult( const double MinimumValue, const ResultParameterSet * FittedParameters, const int FitStatus, const PhysicsBottle* FittedBottle,
		const RapidFitMatrix* CovarianceMatrix, const vector< FunctionContour* > ContourPlots ) :
	minimumValue( MinimumValue ), fittedParameters( new ResultParameterSet(*FittedParameters) ), covarianceMatrix( (CovarianceMatrix==NULL)?(new RapidFitMatrix()):(new RapidFitMatrix(*CovarianceMatrix)) ),
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

vector< FunctionContour* > FitResult::GetContours() const
{
	return contours;
}

RapidFitMatrix* FitResult::GetCovarianceMatrix() const
{
	return covarianceMatrix;
}

void FitResult::ApplyCovarianceMatrix( const RapidFitMatrix* Input )
{
	cout << "Applying FitResult:  ";
	for( unsigned int i=0; i< (unsigned)Input->theseParameters.size(); ++i ) cout << Input->theseParameters[i] << "\t";
	cout << endl;
	if( covarianceMatrix != NULL ) delete covarianceMatrix;
	covarianceMatrix = new RapidFitMatrix(*Input);
	fittedParameters->ApplyCovarianceMatrix( covarianceMatrix );
}

double FitResult::GetMinimumValue() const
{
	return minimumValue;
}

void FitResult::SetMinimumValue( const double min )
{
	minimumValue = min;
}

ResultParameterSet * FitResult::GetResultParameterSet() const
{
	return fittedParameters;
}

int FitResult::GetFitStatus() const
{
	return fitStatus;
}

void FitResult::ForceFitStatus( const int input_status )
{
	fitStatus = input_status;
}

PhysicsBottle* FitResult::GetPhysicsBottle() const
{
	return fittedBottle;
}

void FitResult::Print() const
{
	fittedParameters->Print();
}

