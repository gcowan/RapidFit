/**
  @file main.cpp

  Entry point for RapidFit

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 **/

//  Root Headers
#include "TString.h"
//  RapidFit Headers
#include "Mathematics.h"
#include "FitAssembler.h"
#include "ToyStudy.h"
#include "XMLConfigReader.h"
#include "ResultFormatter.h"
#include "InputParsing.h"
#include "RapidFitIntegrator.h"
#include "Plotter.h"
#include "MakeFoam.h"
#include "PerEventAngularAcceptance.h"
#include "OutputConfiguration.h"
#include "FitResultVector.h"
#include "main.h"
#include "DataSetConfiguration.h"
#include "IDataSet.h"
#include "StringProcessing.h"
#include "MCStudy.h"
#include "GoodnessOfFit.h"
#include "EdStyle.h"
#include "MemoryDataSet.h"
#include "PhaseSpaceBoundary.h"
#include "VectorScan.h"
//  System Headers
#include <string>
#include <vector>
#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

//	Standard opening into RapidFit when not in CINT/pyROOT
#ifndef __CINT__
int main( int argc, char * argv[] )
{
	return RapidFit( argc, argv );
}
#endif

//	The 'meat' of the Code
int RapidFit( int argc, char * argv[] )
{
	//Welcome blurb
	time_t timeNow;
	time(&timeNow);

	cout << endl << "RapidFit" << endl;
	cout << "SVN Rev:\t" << STR(SVN_REV) << endl;
	string mod="M";
	size_t found = string( STR(SVN_REV) ).find( mod );

	if( found != string::npos )
	{
		cout << "!!YOU HAVE LOCAL CODE MODIFICATIONS!!" << endl;
	}
	cout << "Build Info:\t" << STR(BUILD_INFO) << endl;
	cout << "Starting time: " << ctime( &timeNow ) << endl << endl;

	//Variables to store command line arguments
	int numberRepeats = 0;
	int Nuisencemodel=2;
        int jobNum = 0;
        int nData = 0;
	string configFileName = "";
	vector<string> parameterTemplates;
	MinimiserConfiguration * theMinimiser=NULL;
	FitFunctionConfiguration * theFunction=NULL;
	string saveOneDataSetFileName = "";
	string plotFileName = "FitPlots.root";
	string pullFileName = "PullPlots.root";
	string LLscanFileName = "LLscanPlots.root";
	string LLcontourFileName = "LLcontourPlots";
	string FCOutputFile = "FCscanOutputData.root";
	vector< PDFWithData* > pdfsAndData;
	vector<int> RuntimeSeed;
	vector<string> Scan_X;
	vector<string> Contour_X;
	vector<string> Contour_Y;
	vector<pair<string, string> >* XMLOverrideList = new vector<pair<string,string> >;
	vector<string> CommandLineParam;
	string MCStepSize;
	string MCStartEntry;
	int OutputLevel=0;
	int OutputLevel2=-1;

	//Flags for which arguments have been received
	bool numberRepeatsFlag = false;
	bool configFileNameFlag = false;
	bool parameterTemplateFlag = false;
	bool theMinimiserFlag = false;
	bool theFunctionFlag = false;
	bool saveOneDataSetFlag = false;
	bool testIntegratorFlag = false;
	bool testPlotFlag = false;
	bool doPlottingFlag = false;
	bool doPullsFlag = false;
	bool doLLscanFlag = false;
	bool doLLcontourFlag = false;
	bool donewLLscanFlag = false;
	bool donewLLcontourFlag = false;
	bool testRapidIntegratorFlag = false;
	bool calculateAcceptanceWeights = false;
	bool calculateAcceptanceWeightsWithSwave = false;
	bool calculatePerEventAcceptance = false;
	bool defineContourFlag = false;
	bool defineScanFlag = false;
	bool doFC_Flag = false;
	bool UUID_Flag = false;
	bool FC_Debug_Flag = false;
	bool BurnToROOTFlag = false;
	bool MCStudyFlag=false;
	bool Force_Scan_Flag=false;
	bool FC_LL_PART_Flag=false;
	bool GOF_Flag=false;
	bool StartAtCenterFlag=true;

	//	This should do something, it doesn't anymore... what was it?	-	Rob Currie
	(void) testRapidIntegratorFlag;


	//	Some parameters with large scope

	//	These have to be initialized early before the use of goto statements in the main file
	//	I re-wrote main.cpp file to make it clearer what was going on and we were over indenting code by the end of the file
	//	the FC code was over 7th level indented and hard to read/maintain

	string currentArgument;
	FitResultVector* _2DResultForFC=NULL;
	FitResultVector* GlobalFitResult=NULL;
	FitResult * GlobalResult = NULL;
	FitResultVector* FC_Input=NULL;
	vector< ConstraintFunction* > XMLConstraints;
	OutputConfiguration * makeOutput = NULL;
	vector<ParameterSet*> argumentParameterSet;
	XMLConfigReader* xmlFile=NULL;
	MCStudy* newMCStudy = NULL;
	vector<FitResultVector*> SoloContourResults;
	TString FC_Input_File;
	int FC_INPUT_MIN=0, FC_INPUT_MAX=0;
	TString FC_INPUT_PARAM1, FC_INPUT_PARAM2;

	int BAD_COMMAND_LINE_ARG = -39;
	int BAD_XML = -40;
	int NO_XML = -41;



	if ( argc == 1 )
	{
		cout << endl << "This is the Edinburgh, RapidFit Fitter!" << endl;
		cout << endl <<"for the priority of command line arguments run:" << endl <<endl<<"\t"<< argv[0] << " --about" << endl;
		cout << endl <<"run: " << argv[0] << " --help\t for more information" << endl << endl;
		goto exit_RapidFit;
	}

	//Parse command line arguments

	for ( int argumentIndex = 1; argumentIndex < argc; ++argumentIndex )
	{
		currentArgument = argv[argumentIndex];

		if( currentArgument == "--about" )
		{
			cout << endl <<"RapidFit is a fitter which is configured by a control XML to be provided at runtime" << endl << endl;

			cout << "RapidFit has many arguments, some which have higher priority than others." << endl;
			cout << endl <<"The order RapidFit will execute different arguments is:" << endl <<endl;

			cout << "\t 1)	save One Data Set and exit" <<endl;
			cout << "\t 2)	test Integrator and exit" << endl;
			cout << "\t 3)	calculate Acceptance Weights and exit" << endl;
			cout << "\t 4)	calculate Acceptance Weights With Swave and exit" << endl;
			cout << "\t 5)	calculate Per Event Acceptance and exit" << endl;
			cout << "\t 6)	test Projection Plots and exit" << endl;
			cout << "\t 7)	perform a Toy Study and exit" << endl;
			cout << "\t 8)	perform an MC Study and exit" << endl;
			cout << "\t 9)	Actually perform a fit according to the XML" << endl;
			cout << "\t 10)	If requested perform 1D LL Scans" << endl;
			cout << "\t 11)	If requested perform 2D LL Scans" << endl;
			cout << "\t 12)	Perform a FC scan if requested" << endl;
			cout << "\t 13)	Universal Exit" << endl;
			cout << endl;
			cout << endl;
			cout << "For more information on the arguments RapidFit takes: " << endl<<endl;
			cout << "\t\t" << argv[0] << "\t --help" << endl << endl;
			goto exit_RapidFit;
		}
		else if ( currentArgument == "--help" )
		{
			++argumentIndex;

			cout << "QUICK HELP FOR RAPIDFIT" << endl ;

			cout << endl;
			cout << "--about" << endl;
			cout << "	Gives some quick info on the priority of Command Line arguments in RapidFit" << endl;

			cout << endl ;
			cout << " -f <filename>   " << endl ;
			cout << "	Specifies the XML file to drive RapidFit." <<endl ;

			cout << endl ;
			cout << " -repeats n   " << endl ;
			cout << "	Specifies the number of repeats for a Systematic study." <<endl ;

			cout << endl ;

			cout << " --doLLscan  " << endl;
			cout << "	Causes a set of LL scans to be performed around the fit minimum" << endl ;
			cout << "	The parameters which will be scanned must be specified in the xml file as " << endl ;
			cout << "	The range may be defined as hard numerical limits as shown below or as a " <<endl;
			cout << "	Range of n sigma from the fit minima, as shown in Y_Param for LLcontour" << endl;
			cout << "\n	<Output> " << endl ;
			cout << "		<Scan>" << endl;
			cout << "			<Name>parameterName</Name> " << endl ;
			cout << "			<Maximum>parameterMax</Maximum> " << endl ;
			cout << "			<Minimum>parameterMin</Minimum> " << endl ;
			cout << "			<Points>parameterPoints</Points> " << endl ;
			cout << "		</Scan>" << endl;
			cout << "	</Output>" << endl ;

			cout << endl ;

			cout << " --doLLcontour  " << endl;
			cout << "	Causes a set of LL scans to be performed around the fit minimum" << endl ;
			cout << "	The parameters which will be scanned must be specified in the xml file as " << endl ;
			cout << "\n	<Output> " << endl ;
			cout << "		<2DScan>" << endl;
			cout << "			<X_Param>" << endl;
			cout << "				<Name>parameterName</Name> " << endl ;
			cout << "				<Maximum>Maxima of x axis</Maximum> " << endl ;
			cout << "				<Minimum>Minima of x axis</Minimum> " << endl ;
			cout << "				<Points>number of points in x axis</Points> " << endl ;
			cout << "			</Y_Param>" <<endl;
			cout << "				<Name>parameterName</Name> " << endl ;
			cout << "				<Sigma>number of sigma from fit minima</Sigma> " << endl ;
			cout << "				<Points>number of points in y axis</Points> " << endl ;
			cout << "			</Y_Param>" <<endl;
			cout << "		</2DScan>" << endl;
			cout << "	</Output> " << endl ;

			cout << endl ;
			cout << " --saveOneDataSet <filename>   " << endl ;
			cout << "	Causes one Toy dataset to be written out to a file " <<endl ;
			cout << "	filename = file to write the dataset to" << endl ;

			cout << endl ;
			cout << " --testIntegrator   " << endl ;
			cout << "	Useful feature which only tests the numerical<=>analytic integrator for each PDF then exits " <<endl ;

			cout << endl;
			cout << " --SetSeed 12345" << endl;
			cout << "	Set the Random seed to 12345 if you wish to make the output reproducable. Useful on Batch Systems" << endl;

			cout << endl;
			cout << "--PhysParam Phi_s,-3.1,3.1,5,Free" << endl;
			cout << "	Define Phi_s at runtime to be 5 with a range -3.1 < Phi_s < 3.1 and to be Free" <<endl;
			cout << "	This takes priority over what is in the XML file" << endl;

			cout << endl;
			cout << "--defineContour param1,min1,max1,res1 param2,min2,max2,res2" << endl;
			cout << "	This defines a contour in the above parameters with the above ranges and resolutions" << endl;
			cout << "	Using this forces RapidFit to ignore everything else in the XML relating to contours" <<endl;
			cout << "	NB: When scripting be sure to pass strings to these parameters especially coming from python!!!" <<endl;

			cout << endl;
			cout << "--MCStudy" << endl;
			cout << "	Perform an MC style Study which takes an Ntuple and sequentially processes it in sequential steps" << endl;

			cout << endl;
			cout << "--MCStepSize <step_size>" << endl;
			cout << "	The Step size that should be taken by the MCStudy when running over a large NTuple" << endl;
			cout << "	This can either be a single step size (useful for 1 large dataset)," << endl;
			cout << "		 and so <step_size is just a number." << endl;
			cout << "	Or, this can be a unique step size for each ntuple mentioned in a file (useful for multiple datasets)" << endl;
			cout << "\n	i.e.	--MCStepSize	5000		to step through a file every 5000 entries" << endl;
			cout << "		--MCStepSize	2000,3000	step through the first file every 2000 entries and the second every 3000" << endl;

			cout << endl;
			cout << "--MCStartEntry <start_entry>" << endl;
			cout << "	This Start Entry decides where in a file an MC Style Study should start (useful for breaking up the workload)" << endl;
			cout << "	This has a similar structure to the MCStepSize and the style is the same" << endl;
			cout << "\n	i.e.	--MCStartEntry	3000		to start an MC style Study at the 3000th entry in the file(s)" << endl;
			cout << "		--MCStartEntry	30000,50000,222	to start an MC style Study at the 30000th,50000th & 222nd entries in the files" <<endl;

			cout << endl;
			cout << "--ForceScan" << endl;
			cout << "	This forces a scan to be performed on a less than perfectly fit dataset" << endl;

			cout << endl;
			goto exit_RapidFit;
		}
		else if ( currentArgument == "-f" )
		{
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				configFileName = argv[argumentIndex];
				configFileNameFlag = true;
			}
			else
			{
				cerr << "Configuration file name not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "-p" )
		{
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				parameterTemplates.push_back( argv[argumentIndex] );
				parameterTemplateFlag = true;
			}
			else
			{
				cerr << "Parameter template not supplied" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "-minimiser" )
		{
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				theMinimiser = new MinimiserConfiguration( argv[argumentIndex] );
				theMinimiserFlag = true;
			}
			else
			{
				cerr << "Minimiser name not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "-function" )
		{
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				theFunction = new FitFunctionConfiguration( argv[argumentIndex] );
				theFunctionFlag = true;
			}
			else
			{
				cerr << "Function name not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--saveOneDataSet" )
		{
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				saveOneDataSetFileName = argv[argumentIndex];
				saveOneDataSetFlag = true;
			}
			else
			{
				cerr << "Data file name not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "-PDF" )
		{
			if ( argumentIndex + 3 < argc )
			{
				++argumentIndex;
				string pdfName = argv[argumentIndex];
				++argumentIndex;
				string dataSource = argv[argumentIndex];
				++argumentIndex;
				string phaseSpace = argv[argumentIndex];

				pdfsAndData.push_back( InputParsing::MakePDFWithData( pdfName, dataSource, phaseSpace ) );
			}
			else
			{
				cerr << "PDF not specified correctly: must be PDF name, data source and phase space prototype" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "-repeats" )
		{
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				numberRepeats = atoi( argv[argumentIndex] );
				numberRepeatsFlag = true;
			}
			else
			{
				cerr << "Number of repeats not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--OverrideXML" )
		{
			if ( argumentIndex + 2 < argc )
			{
				++argumentIndex;
				string path = argv[argumentIndex];
				++argumentIndex;
				string value = argv[argumentIndex];
				pair<string,string> temp_pair;
				temp_pair.first = path;
				temp_pair.second = value;
				cout << endl;
				cout << "WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!WARNING" << endl;
				cout << "Overriding an XMLTag with:" << endl;
				cout << "\t" << path << "\t" << value << endl;
				cout << "NB:\t this is an expert option for complex situations which haven't been forseen in the design yet" << endl;
				cout << "WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!WARNING" << endl << endl;
				XMLOverrideList->push_back( temp_pair );
			}
			else
			{
				cerr << "XML Override not correctly defined!" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--PhysParam" )
		{
			if ( argumentIndex +1 < argc )
			{
				++argumentIndex;
				string PhysParam = argv[argumentIndex];
				CommandLineParam.push_back( PhysParam );
			}
			else
			{
				cerr << "PhysParam not correctly defined at Runtime" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--testPlot" )
		{
			testPlotFlag = true;
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				plotFileName = argv[argumentIndex];
			}
			else
			{
				cerr << "Path to store plot file not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--doPlotting" )
		{
			doPlottingFlag = true;
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				plotFileName = argv[argumentIndex];
			}
			else
			{
				cerr << "Path to store plot file not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--doPulls" )
		{
			doPullsFlag = true;
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				pullFileName = argv[argumentIndex];
			}
			else
			{
				cerr << "Path to store pull plot file not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--defineScan" )
		{
			if ( argumentIndex + 1 < argc )
			{
				defineScanFlag = true;
				Scan_X.push_back( argv[argumentIndex] );
			} else {
				cerr << "Scan Not Correctly Formatted" <<endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--defineContour" )
		{
			if ( argumentIndex + 2 < argc )
			{
				defineContourFlag = true;
				++argumentIndex;
				Contour_X.push_back( argv[argumentIndex] );
				++argumentIndex;
				Contour_Y.push_back( argv[argumentIndex] );
			} else {
				cerr << "Contour Not Correctly Formatted" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--NuisenceModel" )
		{
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				Nuisencemodel = atoi( argv[argumentIndex] );
			} else {
				cerr << "Not correctly Formatted Nuisence Input" <<endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--SetSeed" )
		{
			if ( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				RuntimeSeed.push_back( atoi( argv[argumentIndex] ) );
			} else {
				cerr << "Seed Not Correctly Defined at Runtime" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--MCStepSize" )
		{
			if( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				MCStepSize = ( argv[argumentIndex] );
			} else {
				cerr << "Badly Defined Step size" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--MCStartEntry" )
		{
			if( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				MCStartEntry = ( argv[argumentIndex] );
			} else {
				cerr << "Badly Defined Start Entry" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--OutputLevel" )
		{
			if( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				OutputLevel = atoi( argv[argumentIndex] );
			} else {
				cerr << "Badly Defined Output Level" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--OutputLevelScans" )
		{
			if( argumentIndex + 1 < argc )
			{
				++argumentIndex;
				OutputLevel2 = atoi( argv[argumentIndex] );
			} else {
				cerr << "Badly Defined Output Level" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if ( currentArgument == "--FC_LL_PART" )
		{
			FC_LL_PART_Flag = true;
			if( argumentIndex + 5 < argc )
			{
				++argumentIndex;
				FC_Input_File = ( argv[argumentIndex] );
				++argumentIndex;
				FC_INPUT_MIN = atoi( argv[argumentIndex] );
				++argumentIndex;
				FC_INPUT_MAX = atoi( argv[argumentIndex] );
				++argumentIndex;
				FC_INPUT_PARAM1 = argv[argumentIndex];
				++argumentIndex;
				FC_INPUT_PARAM2 = argv[argumentIndex];
			} else {
				cerr << "Please Provide a result file with the result of a 2DLL Scan." << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
                else if ( currentArgument == "--GOF" )
                {
                        GOF_Flag = true;
                        if ( argumentIndex + 2 < argc )
                        {
                                ++argumentIndex;
                                jobNum = atoi( argv[argumentIndex] );
                                ++argumentIndex;
                                nData  = atoi( argv[argumentIndex] );
                        }
                        else
                        {
                                cerr << "Job number for GOF not specified" << endl;
                                return BAD_COMMAND_LINE_ARG;
                        }
                }

		//	The Parameters beyond here are for setting boolean flags
		else if ( currentArgument == "--testIntegrator" )			{	testIntegratorFlag = true;			}
		else if ( currentArgument == "--testRapidIntegrator" )			{	testRapidIntegratorFlag = true;			}
		else if ( currentArgument == "--calculateAcceptanceWeights" )		{	calculateAcceptanceWeights = true;		}
		else if ( currentArgument == "--calculateAcceptanceWeightsWithSwave" )	{	calculateAcceptanceWeightsWithSwave = true;	}
		else if ( currentArgument == "--calculatePerEventAcceptance" )		{	calculatePerEventAcceptance = true;		}
		else if ( currentArgument == "--doLLscan" )				{	doLLscanFlag = true;				}
		else if ( currentArgument == "--donewLLscan" )                          {       donewLLscanFlag = true;                         }
		else if ( currentArgument == "--doLLcontour" )				{	doLLcontourFlag = true;				}
		else if ( currentArgument == "--donewLLcontour" )                       {       donewLLcontourFlag = true;                      }
		else if ( currentArgument == "--doFCscan" )				{	doFC_Flag = true;				}
		else if ( currentArgument == "--useUUID" )				{	UUID_Flag = true;				}
		else if ( currentArgument == "--debugFC" )				{	FC_Debug_Flag = true;				}
		else if ( currentArgument == "--BurnToROOT" )				{	BurnToROOTFlag = true;				}
		else if ( currentArgument == "--MCStudy" )				{	MCStudyFlag = true;				}
		else if ( currentArgument == "--ForceScan" )				{	Force_Scan_Flag = true;				}
		else if ( currentArgument == "--DontStartAtCenter" )			{	StartAtCenterFlag = false;			}

		//	We didn't understand the argument to end up here
		else
		{
			cerr << "Unrecognised argument: " << currentArgument << endl;
			exit(BAD_COMMAND_LINE_ARG);
		}
	}


	//Load a config file
	if ( configFileNameFlag )
	{
		xmlFile = new XMLConfigReader(configFileName, XMLOverrideList);
		if ( xmlFile->IsLoaded() )
		{
			cout << "XML config file " << configFileName << " loaded" << endl;
		}
		else
		{
			cerr << "XML config file " << configFileName << " not found" << endl;
			return BAD_XML;
		}
	}


	//	CONFIGURE RAPIDFIT PARAMETERS


	//	I am unaware of anything reasonable you can do with RapidFit without providing a driving XML
	if ( !configFileNameFlag )
	{
		cerr << "No Input XML defined" << endl;
		return NO_XML;
	}

	//	This Section of Code deals with the seed value that is to be used to start whatever your studying
	if( ( !RuntimeSeed.empty() && xmlFile->IsLoaded() ) && !UUID_Flag )
	{
		cout << "Setting Seed At Runtime to be: " << RuntimeSeed[0] << endl;
		xmlFile->SetSeed( unsigned(RuntimeSeed[0]) );
	} else if ( !RuntimeSeed.empty() && UUID_Flag && xmlFile->IsLoaded() )	{
		//  I have used the TRandom3 code as inspiration for 'Salting the Seed' of the UUID at runtime
		cout << "Using a Random seed of UUID x RuntimeSeed to be 'unique' for all machines"<<endl;
		TUUID uid;
		UChar_t uuid[16];
		uid.GetUUID(uuid);
		RuntimeSeed[0] = RuntimeSeed[0] * ( uuid[ 2*(RuntimeSeed[0]%8) ]*256 +uuid[ 2*(RuntimeSeed[0]%8) ] );
		xmlFile->SetSeed( unsigned(RuntimeSeed[0]) );
	}

	//	Command line arguments are passed and interpreted within the parser to override what is read from the XMLFile
	if ( parameterTemplateFlag )	{	argumentParameterSet.push_back( InputParsing::MakeParameterSet( parameterTemplates ) );		}

	//Actually configure a fit: first configure the output
	makeOutput = xmlFile->GetOutputConfiguration();

	//	Command line arguments override the config file

	//	If the Minimiser wasn't defined at runtime consult the XML
	if ( !theMinimiserFlag )
	{
		theMinimiser = xmlFile->GetMinimiserConfiguration();
	}
	//	If the range for a Contour Scan was provided at runtime
	if ( defineContourFlag )
	{
		for( unsigned short int i=0; i < Contour_X.size(); ++i)
		{
			makeOutput->Clear2DScanList();
			makeOutput->AddContour( Contour_X[i], Contour_Y[i] );
		}
	}
	//	If the range for a Scan was provided at runtime
	if ( defineScanFlag )
	{
		for( unsigned short int i=0; i < Scan_X.size(); ++i )
		{
			makeOutput->AddScan( Scan_X[i] );
		}
	}
	//	The Minimisation Function hasn't been defined yet, request one from the XML
	if ( !theFunctionFlag )
	{
		theFunction = xmlFile->GetFitFunctionConfiguration();
		// If weights were specified then we need to let the output plotting know
		if( theFunction->GetWeightsWereUsed() ) makeOutput->SetWeightsWereUsed( theFunction->GetWeightName() ) ;
	}
	//	No Parameter Template had been provided, look in the XML
	if ( !parameterTemplateFlag )
	{
		argumentParameterSet = xmlFile->GetFitParameters( CommandLineParam );
	}
	//	The number of repeats hasn't been defined at runtime, look in the XML
	if ( !numberRepeatsFlag )
	{
		numberRepeats = xmlFile->GetNumberRepeats();
	}
	//	PDFs used in fits and the dataset to be fit to hasn't been defined at the command line
	if ( pdfsAndData.size() == 0 )
	{
		//	Read in from XML
		pdfsAndData = xmlFile->GetPDFsAndData();
		//	If we are performing a scan we want to check for Data Generation instances and generate/store the data in a cache for future use
		if( doLLscanFlag || ( doLLcontourFlag || doFC_Flag ) )
		{
			//	Loop over all ToFits containing data
			for( unsigned int pdf_num=0; pdf_num< pdfsAndData.size(); ++pdf_num )
			{
				//	Get the DataSet in this PDFWithData
				vector<DataSetConfiguration*> DataConfigs = pdfsAndData[pdf_num]->GetAllDataSetConfigs();
				//	Loop over all DataSetConfiguration objects existing within this PDFWithData
				for( unsigned int config_num=0; config_num < DataConfigs.size(); ++config_num )
				{
					//	If it's a file the data is never changed/destroyed
					if( DataConfigs[config_num]->GetSource() != "File" )
					{
						cout << "SCAN REQUESTED, GENERATING AND CACHING DATA" << endl;
						pdfsAndData[pdf_num]->SetDelete( false );			//	Don't remove cache after one use
						pdfsAndData[pdf_num]->SetPhysicsParameters( xmlFile->GetFitParameters() );
						vector<IDataSet*> gen_data;
						gen_data.push_back( pdfsAndData[pdf_num]->GetDataSet() );	//	Generate a Foam DataSet
						pdfsAndData[pdf_num]->AddCachedData( gen_data );		//	Cache it for the lifetime of a scan
					}
				}
			}
		}
	}

	//	Do Pulls?
	if ( doPullsFlag )
	{
		makeOutput->SetPullFileName(pullFileName);
	} else {
		doPullsFlag = makeOutput->DoPullPlots();
	}
	//	Do Plotting
	if ( doPlottingFlag )
	{
		makeOutput->MakeAllPlots(plotFileName);
	}

	//	Seperate MC portion and fitting portion of a FC study
	if ( FC_LL_PART_Flag )
	{
		if( xmlFile == NULL )
		{
			cerr << "I'm not that good, please provide me with an XML" <<endl;
			exit(-3546);
		}
		FC_Input = DataSetConfiguration::LoadFitResult( FC_Input_File, argumentParameterSet.back() );
		GlobalResult = FC_Input->GetFitResult( 0 );
		GlobalFitResult = new FitResultVector( FC_Input->GetAllNames() );
		GlobalFitResult->AddFitResult( GlobalResult, false );
		GlobalFitResult->AddCPUTime( FC_Input->GetCPUTime(0) );
		GlobalFitResult->AddRealTime( FC_Input->GetRealTime(0) );

		vector<FitResultVector*> FC_Inputz;
		for( int i= FC_INPUT_MIN; i<= FC_INPUT_MAX; ++i )
		{
			FitResultVector* temp_vec = new FitResultVector( FC_Input->GetAllNames() );
			temp_vec->AddFitResult( FC_Input->GetFitResult( i ), false );
			temp_vec->AddCPUTime( FC_Input->GetCPUTime( i ) );
			temp_vec->AddRealTime( FC_Input->GetRealTime( i ) );
			FC_Inputz.push_back( temp_vec );
		}
		_2DResultForFC = new FitResultVector( FC_Inputz );
		pdfsAndData = xmlFile->GetPDFsAndData();
		for( unsigned int i=0; i< pdfsAndData.size(); ++i )
		{
			pdfsAndData[i]->SetPhysicsParameters( xmlFile->GetFitParameters() );
		}

		//	All Information about where to perform scans is stored in the _2DResultForFC object
		//
		//	Apart from the names of the controls which the user has to provide
		makeOutput->Clear2DScanList();
		FC_INPUT_PARAM1.Append(",0,0,0");
		FC_INPUT_PARAM2.Append(",0,0,0");
		makeOutput->AddContour( FC_INPUT_PARAM1.Data(), FC_INPUT_PARAM2.Data() );
	}

	//	FINISHED CONFIGURING RAPIDFIT PARAMETERS



	//	Work out what we want to do, now that everything is configured


	//	1)	save One Data Set and exit
	//	2)	test Integrator and exit
	//	3)	calculate Acceptance Weights and exit
	//	4)	calculate Acceptance Weights With Swave and exit
	//	5)	calculate Per Event Acceptance and exit
	//	6)	test Projection Plots and exit
	//	7)	perform a Toy Study and exit
	//	8)	perform an MC Study and exit

	//	9)	Actually perform a fit according to the XML

	//	10)	If requested perform 1D LL Scans
	//	12)	If requested perform 2D LL Scans
	//	13)	Perform a FC scan if requested

	//	13)	Universal Exit


	//	1)
	//	If we want to save one dataset and have an input XML
	if ( saveOneDataSetFlag )
	{
		//Make a file containing toy data from the PDF
		vector<PDFWithData*> quickDataGen = xmlFile->GetPDFsAndData();
		for( unsigned int i=0; i< quickDataGen.size(); ++i )
		{
			vector<IDataSet*> quickDataSet;
			quickDataGen[i]->SetPhysicsParameters( xmlFile->GetFitParameters( CommandLineParam ) );
			IDataSet* temp_dataSet = quickDataGen[i]->GetDataSet();
			quickDataSet.push_back( temp_dataSet );
			string ext_dot=".";
			vector<string> temp_strings = StringProcessing::SplitString( saveOneDataSetFileName, *(ext_dot.c_str()) );
			TString FileName_Pre_Suffix = StringProcessing::CondenseStrings( temp_strings, 0, int(temp_strings.size() -1) );
			TString number;

			if( quickDataGen.size() > 1){ number.Append("_"); number+=i; }

			TString real_saveOneDataSetFileName = TString( FileName_Pre_Suffix + number + ".root" );
			ResultFormatter::MakeRootDataFile( string(real_saveOneDataSetFileName.Data()), quickDataSet );
			delete temp_dataSet;
		}
		goto exit_RapidFit;

	}

	//	2)
	if (testIntegratorFlag && configFileNameFlag)
	{
		//Compare numerical and analytical integration
		PDFWithData * quickData = xmlFile->GetPDFsAndData()[0];
		quickData->SetPhysicsParameters( xmlFile->GetFitParameters() );
		IDataSet * quickDataSet = quickData->GetDataSet();
		RapidFitIntegrator * testIntegrator = new RapidFitIntegrator( quickData->GetPDF() );
		testIntegrator->Integral( quickDataSet->GetDataPoint(0), quickDataSet->GetBoundary() );
		delete testIntegrator;
		delete quickData;
		goto exit_RapidFit;
	}

	//	3)
	if ( calculateAcceptanceWeights && configFileNameFlag )
	{
		// Calculate the acceptance weights from MC
		PDFWithData * pdfAndData = xmlFile->GetPDFsAndData()[0];
		pdfAndData->SetPhysicsParameters( xmlFile->GetFitParameters() );
		IDataSet * dataSet = pdfAndData->GetDataSet();
		IPDF * pdf = pdfAndData->GetPDF();
		Mathematics::calculateAcceptanceWeights(dataSet, pdf);
		delete pdfAndData;
		goto exit_RapidFit;
	}

	//	4)
	if ( calculateAcceptanceWeightsWithSwave && configFileNameFlag )
	{
		// Calculate the acceptance weights from MC
		PDFWithData * pdfAndData = xmlFile->GetPDFsAndData()[0];
		pdfAndData->SetPhysicsParameters( xmlFile->GetFitParameters() );
		IDataSet * dataSet = pdfAndData->GetDataSet();
		IPDF * pdf = pdfAndData->GetPDF();
		Mathematics::calculateAcceptanceWeightsWithSwave(dataSet, pdf);
		delete pdfAndData;
		goto exit_RapidFit;
	}

	//	5)
	if (calculatePerEventAcceptance)
	{
		PerEventAngularAcceptance* a = new PerEventAngularAcceptance("jpsikmc09_loose.root","Bu2JpsiKTuple/DecayTree", "out2.root");
		for (int iter = 1; iter <= 3; ++iter)
		{
			a->fillEffHistos( iter );
			a->loopOnReconstructedBs();
		}
		a->writeHistos();
		delete a;
		goto exit_RapidFit;
	}

	//	6)
	if ( testPlotFlag && configFileNameFlag )
	{
		//Project the PDF onto the data
		PDFWithData * quickData = xmlFile->GetPDFsAndData()[0];
		quickData->SetPhysicsParameters( xmlFile->GetFitParameters() );
		IDataSet * quickDataSet = quickData->GetDataSet();
		Plotter * testPlotter = new Plotter( quickData->GetPDF(), quickDataSet );
		testPlotter->PlotAllObservables(plotFileName);
		delete testPlotter;
		delete quickData;
		goto exit_RapidFit;
	}

	//	7)	Toy Study
	//Pick a toy study if there are repeats, or if pull plots are wanted
	if ( ( ( numberRepeats > 1 ) || doPullsFlag ) && !( ( ( doLLcontourFlag || doFC_Flag ) || doLLscanFlag ) || MCStudyFlag ) )
	{
		vector< ConstraintFunction* > XMLConstraints = xmlFile->GetConstraints();
		//Do the toy study
		ToyStudy* newStudy = new ToyStudy( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, XMLConstraints, numberRepeats );
		newStudy->DoWholeStudy();
		FitResultVector* fitResults = newStudy->GetStudyResult();

		//Output results
		makeOutput->OutputToyResult( fitResults );
		makeOutput->OutputFitResult( fitResults->GetFitResult(0) );

		while( !XMLConstraints.empty() )
		{
			delete XMLConstraints.back();
			XMLConstraints.pop_back();
		}
		delete newStudy;

		goto exit_RapidFit;
	}

	//	8)	MC Study
	if( MCStudyFlag )
	{
		//	Process user input;
		vector<string> true_MCStepSize = StringProcessing::SplitString( MCStepSize, ',' );
		vector<string> true_MCStartEntry = StringProcessing::SplitString( MCStartEntry, ',' );

		//	Convert the components into ints
		vector<int> MCStep_int;
		vector<int> MCStart_int;
		for( vector<string>::iterator step_i = true_MCStepSize.begin(); step_i != true_MCStepSize.end(); ++step_i )
		{
			MCStep_int.push_back( atoi( step_i->c_str() ) );
		}
		for( vector<string>::iterator start_i = true_MCStartEntry.begin(); start_i != true_MCStartEntry.end(); ++start_i )
		{
			MCStart_int.push_back( atoi( start_i->c_str() ) );
		}

		//	Create Toy Study
		newMCStudy = new MCStudy( xmlFile );

		//	Setup Toy Study
		newMCStudy->SetNumRepeats( numberRepeats );
		if( MCStep_int.size() != 0 ) { newMCStudy->SetStartingEntry( MCStart_int ); }
		if( MCStart_int.size() != 0 ) { newMCStudy->SetNumEvents( MCStep_int ); }


		//	Perform Toy Study
		newMCStudy->DoWholeStudy();

		ResultFormatter::WriteFlatNtuple( string( "MC_Study.root" ), newMCStudy->GetStudyResult() );

		delete newMCStudy;
		goto exit_RapidFit;
	}


	//	9)
	if( !FC_LL_PART_Flag )
	{	//	This is for code collapsing and to clearly outline the 'Fit' step of this file
		//		This is re-used for FC scans and forms FC Step 1

		cout << "\n\n\t\tStarting Fit to Find Global Minima!\n"<<endl;

		//Do the fit to find GLOBAL MINIMA
		GlobalFitResult = new FitResultVector( argumentParameterSet.back()->GetAllNames() );
		GlobalFitResult->StartStopwatch();

		XMLConstraints = xmlFile->GetConstraints();

		GlobalResult = FitAssembler::DoSafeFit( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, XMLConstraints, OutputLevel );

		GlobalFitResult->AddFitResult( GlobalResult );

		//ResultFormatter::FlatNTuplePullPlots( string("Global_Fit.root"), GlobalFitResult );

		cout << "\n\n\t\tFit Output:" <<endl;

		if( OutputLevel >= 0 )
		{
			//Output results
			makeOutput->SetInputResults( GlobalResult->GetResultParameterSet() );
			makeOutput->OutputFitResult( GlobalFitResult->GetFitResult(0) );
		}
		ResultFormatter::ReviewOutput( GlobalResult );

		//	If requested write the central value to a single file
		if( BurnToROOTFlag )
		{
			cout << "If you get any errors here and/or an empty Global_Fit_Result.root file, check your xml ONLY contains parameters in the fit." << endl;
			cout << "(There is a known bug which requires a bit of work to fix)" << endl;
			ResultFormatter::WriteFlatNtuple( string( "Global_Fit_Result.root" ), GlobalFitResult );
		}

		if ( GOF_Flag ) {
			cout << "Starting GOF" << endl;
			PDFWithData * pdfAndData = xmlFile->GetPDFsAndData()[0];
			pdfAndData->SetPhysicsParameters( xmlFile->GetFitParameters( CommandLineParam ) );
			IPDF * pdf = pdfAndData->GetPDF();
			vector<IPDF*> vectorPDFs;
			vectorPDFs.push_back(pdf);
			IDataSet * data = pdfAndData->GetDataSet();
			PhaseSpaceBoundary * phase = xmlFile->GetPhaseSpaceBoundaries()[0];
			EdStyle * greigFormat = new EdStyle();
			greigFormat->SetStyle();

			// Set up to be able to generate some MC data
			DataSetConfiguration * dataConfig = pdfAndData->GetDataSetConfig();
			dataConfig->SetSource( "Foam" );
			TH1D * pvalueHist = new TH1D("pvalues", "pvalues", 10, 0, 1);
			double pvalue = 0.;
			bool model = true;
			for ( int i = jobNum; i < jobNum + 100; i++ ) {

				cout << "Ensemble " << i << endl;

				// First, generate some data
				pdfAndData->SetPhysicsParameters( xmlFile->GetFitParameters( CommandLineParam ) );
				pdf->SetRandomFunction( i );
				pdf->SetMCCacheStatus( false );
				MemoryDataSet * subset = (MemoryDataSet*)dataConfig->MakeDataSet( phase, pdf, nData );
				vector<IDataSet*> vectorData;
				vectorData.push_back(subset);

				if ( model ) {
					// Use the same PDF for the large MC dataset (the Model approach)
					pdfAndData->SetPhysicsParameters( xmlFile->GetFitParameters( CommandLineParam ) );
				}
				else {
					// Fit the generated data (the Fit I approach)
					FitResult * gofResult = FitAssembler::DoFit( theMinimiser, theFunction, argumentParameterSet, vectorPDFs, vectorData, xmlFile->GetConstraints() );
					vector< ParameterSet * > parSet;
					parSet.push_back(gofResult->GetResultParameterSet()->GetDummyParameterSet());
					pdfAndData->SetPhysicsParameters( parSet );
				}
				// Generate large sample of MC
				pdf->SetRandomFunction( i*i + 1 );
				pdf->SetMCCacheStatus( false );
				MemoryDataSet * mcData = (MemoryDataSet*)dataConfig->MakeDataSet( phase, pdf, 10*nData );

				// Finally calculate the p-value relative to the large MC sample
				pvalue = GoodnessOfFit::pValueFromPoint2PointDissimilarity( subset, mcData );
				//pvalue = GoodnessOfFit::pValueFromPoint2PointDissimilarity( data, mcData );
				pvalueHist->Fill( pvalue );
				cout << "p-value " << pvalue << endl;
				delete subset;
				delete mcData;
			}
			TFile * outputFile = new TFile("pvalues.root", "RECREATE");
			pvalueHist->Write();
			outputFile->Close();
			delete outputFile;
		}

		if( GlobalResult->GetFitStatus() != 3 )
		{
			cerr << "--------------------------------------------------------------" << endl;
			cerr << "---------------------FIT RESULT IS NOT 3----------------------" << endl;
			cerr << "--------------------------------------------------------------" << endl;
			cerr << "--------------------------------------------------------------" << endl;
			cerr << "---------If this is a Foam study, change seed and re-run------" << endl;
			cerr << "--------------------------------------------------------------" << endl;
			cerr << "--------------------------------------------------------------" << endl;
			cerr << "-If your sure you want to continue employ the following flag--" << endl;
			cerr << "--------------------------------------------------------------" << endl;
			cerr << "--------------      \'--ForceScan\'            -----------------" << endl;
			cerr << "--------------------------------------------------------------" << endl;
			if( !Force_Scan_Flag || ( GlobalResult->GetMinimumValue() < 0 )  )
			{
				goto exit_RapidFit;
			}
		}
	}
	//	10)
	//	Do LL scan
	if( doLLscanFlag )
	{
		vector<FitResultVector*> scanSoloResults;

		//  Store
		vector<string> LLscanList = makeOutput->GetScanList();

		for(unsigned int scan_num=0; scan_num < LLscanList.size() ; ++scan_num)
		{
			if( StartAtCenterFlag )
			{
				vector<ParameterSet*>  param_set;
				param_set.push_back( GlobalResult->GetResultParameterSet()->GetDummyParameterSet() );
				for(unsigned int i=0; i< argumentParameterSet.back()->GetAllNames().size(); ++i )
				{
					argumentParameterSet.back()->GetPhysicsParameter( argumentParameterSet.back()->GetAllNames()[i] )->SetBlindedValue( param_set.back()->GetPhysicsParameter( argumentParameterSet.back()->GetAllNames()[i] )->GetValue() );
				}
			}
			FitResultVector* scan_result = FitAssembler::SingleScan( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints(), makeOutput, LLscanList[scan_num], OutputLevel2 );
			scanSoloResults.push_back( scan_result );
		}

		for(unsigned int scan_num=0; scan_num < LLscanList.size(); ++scan_num )
		{
			TString output_scan_dat( "LLScanData" );
			output_scan_dat.Append( LLscanList[scan_num] );
			TString new_output_scan_dat( output_scan_dat );
			new_output_scan_dat.Append( "_newformat" );
			output_scan_dat.Append( ".root" );
			new_output_scan_dat.Append(".root");
			vector<FitResultVector*> ammended_format;
			for( int i=0; i< scanSoloResults[scan_num]->NumberResults(); ++i )
			{
				ammended_format.push_back( GlobalFitResult );
				FitResultVector* temp_vec = new FitResultVector( scanSoloResults[scan_num]->GetAllNames() );
				temp_vec->AddFitResult( scanSoloResults[scan_num]->GetFitResult( i ), false );
				temp_vec->AddRealTime( scanSoloResults[scan_num]->GetRealTime(i) );
				temp_vec->AddCPUTime( scanSoloResults[scan_num]->GetCPUTime(i) );
				ammended_format.push_back( temp_vec );
			}
			FitResultVector* corrected_format = new FitResultVector( ammended_format );
			cout << output_scan_dat << "\t\t" << corrected_format->NumberResults() << endl;
			ResultFormatter::WriteFlatNtuple( string( output_scan_dat ), scanSoloResults[scan_num] );
			ResultFormatter::WriteFlatNtuple( string( new_output_scan_dat ), corrected_format );
		}
	}

	//	New LL scan technology
	if( donewLLscanFlag )
	{
		vector<FitResultVector*> scanSoloResults;

		//  Store
		vector<string> LLscanList = makeOutput->GetScanList();

		for(unsigned int scan_num=0; scan_num < LLscanList.size() ; ++scan_num)
		{

			//	Start at the global values
			vector<ParameterSet*>  param_set;
			param_set.push_back( GlobalResult->GetResultParameterSet()->GetDummyParameterSet() );
			for(unsigned int i=0; i< argumentParameterSet.back()->GetAllNames().size(); ++i )
			{	argumentParameterSet.back()->GetPhysicsParameter( argumentParameterSet.back()->GetAllNames()[i] )->SetBlindedValue( param_set.back()->GetPhysicsParameter( argumentParameterSet.back()->GetAllNames()[i] )->GetValue() );	}

			ScanParam* local_param = makeOutput->GetScanParam( LLscanList[scan_num] );

			FitResultVector* scan_result = VectorScan::Scan1D( theMinimiser, makeOutput->GetScanParam( LLscanList[scan_num] ) );

			scanSoloResults.push_back( scan_result );
		}

		for(unsigned int scan_num=0; scan_num < LLscanList.size(); ++scan_num )
		{
			TString output_scan_dat( "LLScanData" );
			output_scan_dat.Append( LLscanList[scan_num] );
			TString new_output_scan_dat( output_scan_dat );
			new_output_scan_dat.Append( "_newformat" );
			output_scan_dat.Append( ".root" );
			new_output_scan_dat.Append(".root");
			vector<FitResultVector*> ammended_format;
			for( int i=0; i< scanSoloResults[scan_num]->NumberResults(); ++i )
			{
				ammended_format.push_back( GlobalFitResult );
				FitResultVector* temp_vec = new FitResultVector( scanSoloResults[scan_num]->GetAllNames() );
				temp_vec->AddFitResult( scanSoloResults[scan_num]->GetFitResult( i ), false );
				temp_vec->AddRealTime( scanSoloResults[scan_num]->GetRealTime(i) );
				temp_vec->AddCPUTime( scanSoloResults[scan_num]->GetCPUTime(i) );
				ammended_format.push_back( temp_vec );
			}
			FitResultVector* corrected_format = new FitResultVector( ammended_format );
			cout << output_scan_dat << "\t\t" << corrected_format->NumberResults() << endl;
			ResultFormatter::WriteFlatNtuple( string( output_scan_dat ), scanSoloResults[scan_num] );
			ResultFormatter::WriteFlatNtuple( string( new_output_scan_dat ), corrected_format );
		}
	}

	if( donewLLcontourFlag )
	{
		vector<pair<string, string> > _2DLLscanList = makeOutput->Get2DScanList();
		for(unsigned int i=0; i < _2DLLscanList.size() ; ++i )
		{
			string name1 = _2DLLscanList[i].first;
			string name2 = _2DLLscanList[i].second;

			pair< ScanParam*, ScanParam* > Param_Set = makeOutput->Get2DScanParams( name1, name2 );

			FitResultVector* scan_result = VectorScan::Scan2D( theMinimiser, Param_Set );

			string output_scan_dat("LLcontourScanData.root");
			if( _2DLLscanList.size() > 0 )
			{
				TString _scan_dat("LLcontourScanData_");
				_scan_dat+= _2DLLscanList.size()+1;
				_scan_dat.Append(".root");
				output_scan_dat=_scan_dat;
			}

			//      Add the Global Results and 'Linearize' the output
			vector<FitResultVector*> TempContourResults;
			TempContourResults.push_back( GlobalFitResult );
			TempContourResults.push_back( scan_result );
			FitResultVector* TempContourResults2 = new FitResultVector( TempContourResults );
			ResultFormatter::WriteFlatNtuple( output_scan_dat , TempContourResults2 );

		}
	}

	//		This is re-used for FC scans and forms FC Step 2

	//	11)
	//Do 2D LL scan for deltaGamma and phis
	if( ( doLLcontourFlag || doFC_Flag ) && ( !FC_LL_PART_Flag ) )
	{
		vector<pair<string, string> > _2DLLscanList = makeOutput->Get2DScanList();

		unsigned int initial_scan=0;
		if( ( ( _2DLLscanList.size() > 1 ) && doFC_Flag ) || defineContourFlag )
		{
			//cerr << "\n\nPERFORMING ONLY ONE 2D SCAN, CHECK THIS IS EXPECTED!" <<endl;
			initial_scan = int( _2DLLscanList.size()-1 );
		}

		if( ( _2DLLscanList.size() == 0 ) && ( doFC_Flag) )
		{
			cerr << "\n\n\t\tNO 2D SCAN DATA, I'M NOT GOING TO DO FC, GO AWAY!\n\n"<<endl;
			exit(-54);
		}

		vector<TString> LLcontourFileNamez;

		for(unsigned int ii=initial_scan; ii < _2DLLscanList.size() ; ++ii )
		{
			string name1 = _2DLLscanList[ii].first;
			string name2 = _2DLLscanList[ii].second;
			//theMinimiser = xmlFile->GetMinimiserConfiguration();
			if( StartAtCenterFlag )
			{
				vector<ParameterSet*> param_set;
				param_set.push_back( GlobalResult->GetResultParameterSet()->GetDummyParameterSet() ); 
				for(unsigned int i=0; i< argumentParameterSet.back()->GetAllNames().size(); ++i )
				{
					argumentParameterSet.back()->GetPhysicsParameter( argumentParameterSet.back()->GetAllNames()[i] )->SetBlindedValue( param_set.back()->GetPhysicsParameter( argumentParameterSet.back()->GetAllNames()[i] )->GetValue() );
				}
			}

			vector<FitResultVector*> Temp_Results = FitAssembler::ContourScan( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints(), makeOutput, name1, name2, OutputLevel2 );

			GlobalResult->GetResultParameterSet()->GetResultParameter( name1 )->ForcePullValue( -9999 );
			GlobalResult->GetResultParameterSet()->GetResultParameter( name1 )->ForceOriginalValue( -9999 );
			GlobalResult->GetResultParameterSet()->GetResultParameter( name2 )->ForcePullValue( -9999 );
			GlobalResult->GetResultParameterSet()->GetResultParameter( name2 )->ForceOriginalValue( -9999 );

			TString ext_string( "-" );
			ext_string.Append( name1 );
			ext_string.Append( "_" );
			ext_string.Append( name2 );
			ext_string.Append( ".root" );
			TString TempName = LLcontourFileName;
			TempName.Append( ext_string );
			LLcontourFileNamez.push_back( TempName );

			//	Linearize the data and then Set the generate values to a default
			SoloContourResults.push_back( new FitResultVector( Temp_Results ) );
			//	This should probably be in the main as I don't want to pollute the FitResult with
			//	'hard coded defaults that seem sensible now'
			for( unsigned short int point_num=0; point_num < SoloContourResults.back()->NumberResults(); ++point_num )
			{
				SoloContourResults.back()->GetFitResult(point_num)->GetResultParameterSet()->GetResultParameter( name1 )->ForcePullValue( -9999 );
				SoloContourResults.back()->GetFitResult(point_num)->GetResultParameterSet()->GetResultParameter( name1 )->ForceOriginalValue( -9999 );
				SoloContourResults.back()->GetFitResult(point_num)->GetResultParameterSet()->GetResultParameter( name2 )->ForcePullValue( -9999 );
				SoloContourResults.back()->GetFitResult(point_num)->GetResultParameterSet()->GetResultParameter( name2 )->ForceOriginalValue( -9999 );
			}
		}

		//	12)
		//  Don't output the scan files when doing a FC scan
		if( doFC_Flag )  {
			_2DResultForFC = SoloContourResults.back();
		} else  {
			for(unsigned int scanNum=0; scanNum < _2DLLscanList.size(); ++scanNum )
			{

				//	This output is A LOT safer as it intended to be fed to a sperate plotting program.
				//	Formatting SHOULD be seperate to generation!
				string output_scan_dat("LLcontourScanData.root");
				if( scanNum > 0 )
				{
					TString _scan_dat("LLcontourScanData_");
					_scan_dat+= scanNum+1;
					_scan_dat.Append(".root");
					output_scan_dat=_scan_dat;
				}

				//	Add the Global Results and 'Linearize' the output
				vector<FitResultVector*> TempContourResults;
				TempContourResults.push_back( GlobalFitResult );
				TempContourResults.push_back( SoloContourResults[scanNum] );
				FitResultVector* TempContourResults2 = new FitResultVector( TempContourResults );
				ResultFormatter::WriteFlatNtuple( output_scan_dat , TempContourResults2 );
			}
		}
	}

	//	12b)
	//	Do the main work of the FC scan
	if( doFC_Flag )
	{
		vector<unsigned int> numberRepeatsVec;
		if( numberRepeatsFlag ) numberRepeatsVec.push_back( unsigned(numberRepeats) );

		//	Do FC scan
		FitResultVector* AllFCResults = FitAssembler::FeldmanCousins( GlobalFitResult, _2DResultForFC, numberRepeatsVec, unsigned(int(Nuisencemodel)), FC_Debug_Flag, makeOutput, theMinimiser, theFunction,  xmlFile, pdfsAndData );

		//		STORE THE OUTPUT OF THE TOY STUDIES
		ResultFormatter::WriteFlatNtuple( "FCOutput.root", AllFCResults );
	}

	//	Should only happen under the condition that no CV fit was performed or anything else
	if( GlobalFitResult == NULL )
	{
		//	Default action - presumably a fit or a toy study
		cerr << "No action performed" << endl;
		cerr << "Not sure how I got here, please email a maintainer!..." <<endl;
	}

	//	13)	Exit
	//	This is executed once everything else has finished
exit_RapidFit:


	while( !XMLConstraints.empty() )
	{
		delete XMLConstraints.back();
		XMLConstraints.pop_back();
	}
	delete GlobalResult;
	delete GlobalFitResult;

	//	These have to be cleaned here such that they are still within scope during FC scans
	while( !SoloContourResults.empty() )
	{
		delete SoloContourResults.back();
		SoloContourResults.pop_back();
	}

	//	Clean UP!
	while ( !pdfsAndData.empty() )
	{
		delete pdfsAndData.back();
		pdfsAndData.pop_back();
	}
	delete xmlFile;


	//Exit blurb
	time(&timeNow);
	cout << endl << "RapidFit" << endl;
	cout << "Ending time: " << ctime( &timeNow ) << endl;

	cout << "Goodbye :)" << endl;
	return 0;
}
