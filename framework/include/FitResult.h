/**
        @class FitResult

        Container for all results from a minimisation

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef FIT_RESULT_H
#define FIT_RESULT_H

#include "ResultParameterSet.h"
#include "PhysicsBottle.h"
#include "FunctionContour.h"
#include <vector>

class FitResult
{
	public:
		FitResult();
		FitResult( double, ResultParameterSet*, int, PhysicsBottle );
		FitResult( double, ResultParameterSet*, int, PhysicsBottle, vector<double> );
		FitResult( double, ResultParameterSet*, int, PhysicsBottle, vector<double>, vector< FunctionContour* >);
		~FitResult();

		double GetMinimumValue();
		vector<double> GetCovarianceMatrix();
		vector< FunctionContour* > GetContours();
		ResultParameterSet * GetResultParameterSet();
		int GetFitStatus();
		void ForceFitStatus( int );
		PhysicsBottle * GetPhysicsBottle();

	private:
		double minimumValue;
		ResultParameterSet * fittedParameters;
		vector<double> covarianceMatrix;
		vector< FunctionContour* > contours;
		int fitStatus;
		PhysicsBottle fittedBottle;
};

#endif
