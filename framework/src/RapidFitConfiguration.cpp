
#include "RapidFitConfiguration.h"

#include <string>
#include <vector>

using namespace::std;

RapidFitConfiguration::RapidFitConfiguration() :
numberRepeats(),
	Nuisencemodel(),
	jobNum(),
	nData(),
	jackStartNum(),
	jackStopNum(),
	configFileName(),
	parameterTemplates(),
	theMinimiser(),
	theFunction(),
	calcConfig(),
	calculator(),
	saveOneDataSetFileName(),
	observableName(),
	plotFileName(),
	pullFileName(),
	LLscanFileName(),
	LLcontourFileName(),
	FCOutputFile(),
	pdfsAndData(),
	RuntimeSeed(),
	Scan_X(),
	Contour_X(),
	Contour_Y(),
	XMLOverrideList(),
	CommandLineParamvector(),
	MCStepSize(),
	MCStartEntry(),
	OutputLevel(),
	OutputLevel2(),
	numberRepeatsFlag(),
	configFileNameFlag(),
	parameterTemplateFlag(),
	theMinimiserFlag(),
	theFunctionFlag(),
	saveOneDataSetFlag(),
	saveOneFoamDataSetFlag(),
	testIntegratorFlag(),
	testComponentPlotFlag(),
	observableNameFlag(),
	doPlottingFlag(),
	doPullsFlag(),
	doLLscanFlag(),
	doLLcontourFlag(),
	testRapidIntegratorFlag(),
	calculateFitFractionsFlag(),
	calculateAcceptanceWeights(),
	calculateAcceptanceCoefficients(),
	calculateAcceptanceWeightsWithSwave(),
	calculatePerEventAcceptance(),
	defineContourFlag(),
	defineScanFlag(),
	doFC_Flag(),
	UUID_Flag(),
	BurnToROOTFlag(),
	MCStudyFlag(),
	Force_Continue_Flag(),
	FC_LL_PART_Flag(),
	GOF_Flag(),
	StartAtCenterFlag(),
	WeightDataSet(),
	OutputLevelSet(),
	saveFitXML(),
	generateToyXML(),
	makeTemplateXML(),
	JackKnife_Flag(),
	fixedTotalToys(),
	saveAllToys(),
	currentArgument(),
	_2DResultForFC(),
	GlobalFitResult(),
	GlobalResult(),
	XMLConstraints(),
	makeOutput(),
	argumentParameterSet(),
	xmlFile(),
	SoloContourResults(),
	templatePDFs(),
	debug(),
	disableLatexOutput(),
	BuildConstraints(),
	runtimeArgs()
{
		//Variables to store command line arguments
		numberRepeats = 0;
		Nuisencemodel=2;
		jobNum = 0;
		nData = 0;
		configFileName = "";
		parameterTemplates = vector<string>();
		theMinimiser=NULL;
		theFunction=NULL;
		calcConfig = NULL;
		calculator = NULL;
		saveOneDataSetFileName = "";
		observableName = "time";
		plotFileName = "FitPlots.root";
		pullFileName = "PullPlots.root";
		LLscanFileName = "LLscanPlots.root";
		LLcontourFileName = "LLcontourPlots";
		FCOutputFile = "FCscanOutputData.root";
		pdfsAndData = vector<PDFWithData*>();
		RuntimeSeed = vector<int>();
		Scan_X = vector<string>();
		Contour_X = vector<string>();
		Contour_Y = vector<string>();
		XMLOverrideList = new vector<pair<string,string> >();
		CommandLineParamvector = vector<string>();
		MCStepSize="";
		MCStartEntry="";
		OutputLevel=0;
		OutputLevel2=-1;

		jackStartNum = 0;
		jackStopNum = 0;

		//Flags for which arguments have been received
		numberRepeatsFlag = false;
		configFileNameFlag = false;
		parameterTemplateFlag = false;
		theMinimiserFlag = false;
		theFunctionFlag = false;
		saveOneDataSetFlag = false;
		saveOneFoamDataSetFlag = false;
		testIntegratorFlag = false;
		testComponentPlotFlag = false;
		observableNameFlag = false;
		doPlottingFlag = false;
		doPullsFlag = false;
		doLLscanFlag = false;
		doLLcontourFlag = false;
		testRapidIntegratorFlag = false;
		calculateFitFractionsFlag = false;
		calculateAcceptanceWeights = false;
		calculateAcceptanceCoefficients = false;
		calculateAcceptanceWeightsWithSwave = false;
		calculatePerEventAcceptance = false;
		defineContourFlag = false;
		defineScanFlag = false;
		doFC_Flag = false;
		UUID_Flag = false;
		BurnToROOTFlag = false;
		MCStudyFlag=false;
		Force_Continue_Flag=false;
		FC_LL_PART_Flag=false;
		GOF_Flag=false;
		StartAtCenterFlag=true;
		WeightDataSet=false;
		OutputLevelSet=false;
		saveFitXML=false;
		generateToyXML=false;
		makeTemplateXML=false;
		JackKnife_Flag=false;
		fixedTotalToys = false;
		saveAllToys = false;
		disableLatexOutput = false;

		_2DResultForFC = NULL;
		GlobalFitResult = NULL;
		GlobalResult = NULL;
		XMLConstraints = vector<ConstraintFunction*>();
		makeOutput = NULL;
		argumentParameterSet = NULL;
		xmlFile = NULL;
		SoloContourResults = vector<FitResultVector*>();

		templatePDFs = vector<string>();
		debug = new DebugClass(false);
		debug->SetClassNames(vector<string>(1,"default"));

		BuildConstraints = false;

		runtimeArgs = vector<string>();
}

RapidFitConfiguration::~RapidFitConfiguration()
{
		if( theMinimiser != NULL ) delete theMinimiser;
		if( theFunction != NULL ) delete theFunction;
		if( calcConfig != NULL ) delete calcConfig;
		if( calculator != NULL ) delete calculator;

		while( !pdfsAndData.empty() )
		{
			if( pdfsAndData.back() != NULL ) delete pdfsAndData.back();
			pdfsAndData.pop_back();
		}

		if( _2DResultForFC != NULL ) delete _2DResultForFC;
		if( GlobalFitResult != NULL ) delete GlobalFitResult;
		if( GlobalResult != NULL ) delete GlobalResult;

		while( !XMLConstraints.empty() )
		{
			if( XMLConstraints.back() != NULL ) delete XMLConstraints.back();
			XMLConstraints.pop_back();
		}

		if( makeOutput != NULL ) delete makeOutput;
		if( argumentParameterSet != NULL ) delete argumentParameterSet;
		if( xmlFile != NULL ) delete xmlFile;

		while( !SoloContourResults.empty() )
		{
			if( SoloContourResults.back() != NULL ) delete SoloContourResults.back();
			SoloContourResults.pop_back();
		}
		ResultFormatter::CleanUp();
		RapidFitIntegrator::clearGSLIntegrationPoints();
}

