
#pragma once
#ifndef _RapidFitConfig_H
#define _RapidFitConfig_H

#include "MinimiserConfiguration.h"
#include "FitFunctionConfiguration.h"
#include "PrecalculatorConfig.h"
#include "IPrecalculator.h"
#include "FitResult.h"
#include "FitResultVector.h"
#include "ConstraintFunction.h"
#include "OutputConfiguration.h"
#include "ParameterSet.h"
#include "XMLConfigReader.h"
#include "MCStudy.h"
#include "DebugClass.h"
#include "ResultFormatter.h"

#include <string>
#include <vector>

using namespace::std;

class RapidFitConfiguration
{
	public:
		RapidFitConfiguration();

		~RapidFitConfiguration();

		//Variables to store command line arguments
		int numberRepeats;
		unsigned int Nuisencemodel;
		int jobNum;
		int nData;
		int jackStartNum;
		int jackStopNum;
		string configFileName;
		vector<string> parameterTemplates;
		MinimiserConfiguration * theMinimiser;
		FitFunctionConfiguration * theFunction;
		PrecalculatorConfig* calcConfig;
		IPrecalculator* calculator;
		string saveOneDataSetFileName;
		string observableName;
		string plotFileName;
		string pullFileName;
		string LLscanFileName;
		string LLcontourFileName;
		string FCOutputFile;
		vector< PDFWithData* > pdfsAndData;
		vector<int> RuntimeSeed;
		vector<string> Scan_X;
		vector<string> Contour_X;
		vector<string> Contour_Y;
		vector<pair<string, string> >* XMLOverrideList;
		vector<string> CommandLineParamvector;
		string MCStepSize;
		string MCStartEntry;
		int OutputLevel;
		int OutputLevel2;

		//Flags for which arguments have been received
		bool numberRepeatsFlag;
		bool configFileNameFlag;
		bool parameterTemplateFlag;
		bool theMinimiserFlag;
		bool theFunctionFlag;
		bool saveOneDataSetFlag;
		bool saveOneFoamDataSetFlag;
		bool testIntegratorFlag;
		bool testComponentPlotFlag;
		bool observableNameFlag;
		bool doPlottingFlag;
		bool doPullsFlag;
		bool doLLscanFlag;
		bool doLLcontourFlag;
		bool testRapidIntegratorFlag;
		bool calculateFitFractionsFlag;
		bool calculateAcceptanceWeights;
		bool calculateAcceptanceCoefficients;
		bool calculateAcceptanceWeightsWithSwave;
		bool calculatePerEventAcceptance;
		bool defineContourFlag;
		bool defineScanFlag;
		bool doFC_Flag;
		bool UUID_Flag;
		bool BurnToROOTFlag;
		bool MCStudyFlag;
		bool Force_Continue_Flag;
		bool FC_LL_PART_Flag;
		bool GOF_Flag;
		bool StartAtCenterFlag;
		bool WeightDataSet;
		bool OutputLevelSet;
		bool saveFitXML;
		bool generateToyXML;
		bool makeTemplateXML;
		bool JackKnife_Flag;
		bool fixedTotalToys;
		bool saveAllToys;
		bool BuildConstraints;
		bool disableLatexOutput;

		string currentArgument;
		FitResultVector* _2DResultForFC;
		FitResultVector* GlobalFitResult;
		FitResult * GlobalResult;
		vector<ConstraintFunction*> XMLConstraints;
		OutputConfiguration* makeOutput;
		ParameterSet* argumentParameterSet;
		XMLConfigReader* xmlFile;
		vector<FitResultVector*> SoloContourResults;

		vector<string> templatePDFs;

		vector<string> runtimeArgs;

		DebugClass* debug;
	private:
		/*!
		 * @brief DO NOT COPY THE CONFIGURATION OBJECT!!!!!
		 */
		RapidFitConfiguration& operator = ( const RapidFitConfiguration& );
		RapidFitConfiguration( const RapidFitConfiguration& );
};

#endif

