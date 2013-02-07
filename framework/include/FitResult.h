/*!
 * @class FitResult
 *
 * Container for all results from a minimisation
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef FIT_RESULT_H
#define FIT_RESULT_H

///	RapidFit Headers
#include "ResultParameterSet.h"
#include "PhysicsBottle.h"
#include "FunctionContour.h"
#include "RapidFitMatrix.h"
///	System Headers
#include <vector>

using namespace::std;

class FitResult
{
	public:
		/*!
		 * Constructor containing the ResultParameterSet and the PhysicsBottle as well as the Covarience Matrix and Contours
		 */
		FitResult( const double, const ResultParameterSet*, const int, const PhysicsBottle*, const RapidFitMatrix* =NULL, const vector<FunctionContour*> =vector<FunctionContour*>());

		/*!
		 * Correct Copy Constructor
		 */
		FitResult ( const FitResult& );

		/*!
		 * Destructor
		 */
		~FitResult();

		/*!
		 * Returns the Minimum total function value i.e. the NLL for this Result
		 */
		double GetMinimumValue() const;
		void SetMinimumValue( double );

		/*!
		 * Return a copy of the Covarience Matrix as provided by the Minimiser
		 */
		RapidFitMatrix* GetCovarianceMatrix() const;

		void ApplyCovarianceMatrix( const RapidFitMatrix* Input );

		/*!
		 * Get the Contours extracted using Minuit's internal Functions
		 */
		vector<FunctionContour*> GetContours() const;

		/*!
		 * Get the Result Parameter Set from this FitReult
		 */
		ResultParameterSet* GetResultParameterSet() const;

		/*!
		 * Get the Final FitStatus from the Minimiser
		 */
		int GetFitStatus() const;

		/*!
		 * Forcibly alter the FitStatus
		 */
		void ForceFitStatus( const int );

		/*!
		 * Get the Pointer to the Physics Bottle from this FitResult
		 */
		PhysicsBottle* GetPhysicsBottle() const;

		/*!
		 * Output Useful Debugging Info
		 */
		void Print() const;

	private:
		/*!
		 * Don't Copy the class this way!
		 */
		FitResult& operator = ( const FitResult& );

		double minimumValue;			/*!	The Final value from the Function After Minimisation		*/
		ResultParameterSet * fittedParameters;	/*!	Pointer to the Result ParameterSet provided at Construction	*/
		RapidFitMatrix* covarianceMatrix;	/*!	The Covariance Matrix as constructed from the Minimiser		*/
		vector< FunctionContour* > contours;	/*!	Contours Provided directly from the Minimiser			*/
		int fitStatus;				/*!	Final Minimiser Status value					*/
		PhysicsBottle* fittedBottle;		/*!	Pointer to the Physics Bottle defined at Construction		*/
};

#endif

