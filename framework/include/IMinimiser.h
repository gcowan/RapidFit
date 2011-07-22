/**
        @interface IMinimiser

        Interface for all function minimisers

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef I_MINIMISER_H
#define I_MINIMISER_H

//	RapidFit Headers
#include "FitFunction.h"
#include "FitResult.h"
#include "PhysicsParameter.h"

class IMinimiser
{
	public:
		virtual ~IMinimiser(){};
		virtual void SetupFit( FitFunction* ) = 0;
		virtual void FixParameters( vector<double>, vector<string> ) = 0;
		virtual void Minimise() = 0;
		virtual FitResult * GetFitResult() = 0;
		virtual void ContourPlots( vector< pair< string, string > > ) = 0;
		virtual void SetOutputLevel( int ) = 0;
		virtual void SetSteps( int ) = 0;
		virtual void SetTolerance( double ) = 0;
		virtual void SetOptions( vector<string> ) = 0;
		virtual void SetQuality( int ) = 0;
		virtual FitFunction* GetFitFunction() = 0;
};

#endif
