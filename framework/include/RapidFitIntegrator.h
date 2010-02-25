/**
	@class RapidFitIntegrator

	A numerical integrator to be used for projecting the PDF, and as a fallback if the PDF does not provide its own normalisation
	This class uses two integrator classes provided by root: AdaptiveIntegratorMultiDim and GaussLegendreIntegrator.
	Both of these assume the function to be integrated is well behaved.

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-8
*/

#ifndef RAPIDFIT_INTEGRATOR_H
#define RAPIDFIT_INTEGRATOR_H

#include "IntegratorFunction.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/Integrator.h"
#include "FoamIntegrator.h"
#include "IDataSet.h"

using namespace ROOT::Math;

class RapidFitIntegrator
{
	public:
		RapidFitIntegrator();
		RapidFitIntegrator( IPDF*, bool ForceNumerical = false );
		~RapidFitIntegrator();

		double Integral( DataPoint*, PhaseSpaceBoundary*, bool UseCache = false );
		double ProjectObservable( DataPoint*, PhaseSpaceBoundary*, string );
		double GetRatioOfIntegrals();
		IPDF * GetPDF();
		void UpdateIntegralCache( PhaseSpaceBoundary* );
	
	private:
		double DoNumericalIntegral( DataPoint*, PhaseSpaceBoundary*, vector<string> );
		double GetCachedIntegral( DataPoint* );
		void SetUpIntegralCache( PhaseSpaceBoundary* );

		double ratioOfIntegrals;
		double cumulativeError;
		double numberCalls;
		bool testFast;
		FoamIntegrator * fastIntegrator;
		//BenIntegrator * fastIntegrator;
		
		IPDF * functionToWrap;
		AdaptiveIntegratorMultiDim * multiDimensionIntegrator;
		IntegratorOneDim * oneDimensionIntegrator;
		bool functionCanIntegrate, functionCanProject, haveTestedIntegral, forceNumerical, cacheSetUp;

		vector<string> discreteNames, continuousNames;
		vector< vector<double> > discreteValues, discreteCombinations;
		vector<double> cachedIntegrals;
};

#endif
