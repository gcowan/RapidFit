
#pragma once
#ifndef RAPIDFITINTEGRATOR_CONFIG_H
#define RAPIDFITINTEGRATOR_CONFIG_H

#define __DEFAULT_RAPIDFIT_FIXEDINTEGRATIONPOINTS 100000
#define __DEFAULT_RAPIDFIT_USEGSL false
#define __DEFAULT_RAPIDFIT_MAXINTEGRALSTEPS 1000000
#define __DEFAULT_RAPIDFIT_INTABSTOL 1E-9
#define __DEFAULT_RAPIDFIT_INTRELTOL 1E-9

#include "Threading.h"

using namespace::std;

class RapidFitIntegratorConfig
{
	public:
		RapidFitIntegratorConfig() :
			FixedIntegrationPoints( __DEFAULT_RAPIDFIT_FIXEDINTEGRATIONPOINTS ), useGSLIntegrator( __DEFAULT_RAPIDFIT_USEGSL ),
			MaxIntegrationSteps( __DEFAULT_RAPIDFIT_MAXINTEGRALSTEPS ), IntegrationAbsTolerance( __DEFAULT_RAPIDFIT_INTABSTOL ),
			IntegrationRelTolerance( __DEFAULT_RAPIDFIT_INTRELTOL ), numThreads( (unsigned)Threading::numCores() )
		{
		}

		unsigned int FixedIntegrationPoints;
		bool useGSLIntegrator;
		unsigned int MaxIntegrationSteps;
		double IntegrationAbsTolerance;
		double IntegrationRelTolerance;
		unsigned int numThreads;
};

#endif

