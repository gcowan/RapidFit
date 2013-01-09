/*!
 * @namespace GoodnessOfFit
 *
 * @brief Where all of teh Code relating to GoodNess of Fit is
 *
 * @author Greig A Cowan greig.cowan@cern.ch
 */

#pragma once
#ifndef RAPIDFIT_GOODNESSOFFIT
#define RAPIDFIT_GOODNESSOFFIT

//	ROOT Headers
#include "TH1D.h"
//	RapidFit Headers
#include "IDataSet.h"
#include "IPDF.h"
#include "PhaseSpaceBoundary.h"
#include "XMLConfigReader.h"
#include "MinimiserConfiguration.h"
#include "FitFunctionConfiguration.h"

namespace GoodnessOfFit
{
	double gofLoop( XMLConfigReader * xmlFile, MinimiserConfiguration * theMinimiser, FitFunctionConfiguration * theFunction, ParameterSet* argumentParameterSet, vector<string> CommandLineParam, int nData );
        double fitDataCalculatePvalue( XMLConfigReader * xmlFile, MinimiserConfiguration * theMinimiser, FitFunctionConfiguration * theFunction, ParameterSet* argumentParameterSet, FitResult * result);
        void generateFitAndCalculatePvalue( XMLConfigReader * xmlFile, ParameterSet* parSet, MinimiserConfiguration * theMinimiser, FitFunctionConfiguration * theFunction, ParameterSet* argumentParameterSet, int nData, int repeats, vector<double> * pvalues);
	double getPvalue( double datavalue, vector<double> distribution );
	double pValueFromPoint2PointDissimilarity( IDataSet * data, IDataSet * mc );
	double calculateTstatistic( IDataSet * data, IDataSet * mc);
	vector<double> permutation( IDataSet * data, IDataSet * mc, int nPerm );
	double permutationCore( IDataSet * data, IDataSet * mc, int iteration );
	double sumEvents( IDataSet * data );
	double sumDataMCEvents( IDataSet * data, IDataSet * mcData );
	void plotUstatistic( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, string plot );
	void calculateUstatistic( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, TH1D * distances);
        void calculateUstatisticNum( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, TH1D * distances);
	double getWeight( DataPoint * x );
	double getDistance( DataPoint * x, DataPoint * y );
	vector<double> getDistances( DataPoint * x, DataPoint * y );
        void copyPhaseSpaceBoundary( PhaseSpaceBoundary * newBoundary, PhaseSpaceBoundary * oldBoundary );
	void updatePhaseSpaceBoundary( DataPoint * x, PhaseSpaceBoundary * newPhase, PhaseSpaceBoundary * oldPhase, vector<double> distances );
        double diffmod2pi( double input );
}

#endif

