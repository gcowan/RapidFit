/*!
 * @namespace JackKnife
 *
 * @brief JackKnifing to the max 
 *
 * @author Greig A Cowan greig.cowan@cern.ch
 */

#pragma once
#ifndef RAPIDFIT_JACKKNIFE
#define RAPIDFIT_JACKKNIFE

//	ROOT Headers
#include "TH1D.h"
//	RapidFit Headers
#include "IDataSet.h"
#include "IPDF.h"
#include "PhaseSpaceBoundary.h"
#include "XMLConfigReader.h"
#include "MinimiserConfiguration.h"
#include "FitFunctionConfiguration.h"

namespace JackKnife {
	void jackknife( XMLConfigReader * xmlFile, MinimiserConfiguration * theMinimiser, 
		FitFunctionConfiguration * theFunction, ParameterSet* argumentParameterSet, vector<string> CommandLineParam, int start, int stop);
	void plotUstatistic( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, string plot );
}

#endif

