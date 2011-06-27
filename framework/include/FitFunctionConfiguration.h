/**
	@class FitFunctionConfiguration

	Container that stores all information related to FitFunction configuration, and returns an appropriate instance of a FitFunction

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-11-27
*/

#ifndef FIT_FUNCTION_CONFIGURATION_H
#define FIT_FUNCTION_CONFIGURATION_H

//	ROOT Headers
#include "TString.h"
//	RapidFit Headers
#include "FitFunction.h"
#include "PhysicsBottle.h"

class FitFunctionConfiguration
{
	public:
		FitFunctionConfiguration();
		FitFunctionConfiguration( string );
		FitFunctionConfiguration( string, string );
		~FitFunctionConfiguration();

		FitFunction * GetFitFunction( PhysicsBottle* );
	
		bool GetWeightsWereUsed() ;
		string GetWeightName() ;

		void SetupTrace( TString );
	private:
		string functionName, weightName;
		bool hasWeight;
		bool wantTrace;
		TString TraceFileName;
};

#endif
