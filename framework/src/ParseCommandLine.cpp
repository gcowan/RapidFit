
#include "RapidFitConfiguration.h"
#include "ParseCommandLine.h"
#include "InputParsing.h"
#include "StringProcessing.h"
#include "ComponentPlotter.h"

#include <vector>
#include <string>
#include <iostream>

using namespace::std;

void ParseCommandLine::RapidFitHelp()
{
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
	cout << "		<TwoDScan>" << endl;
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
	cout << "		</TwoDScan>" << endl;
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
	cout << "--calculateAcceptanceWeights" << endl;
	cout << "       Calculate the Angular Acceptance Weights and Histograms for both Helicity and Transverse Basis and store the results in a .root file" << endl;

	cout << endl;
	cout << "--calculateAcceptanceCoefficients" << endl;
	cout << "       Calculate the Angular Acceptance Coefficients for Helicity Basis" << endl;

	cout << endl;
	cout << "--WeightDataSet" << endl;
	cout << "       If a Preprocessor section is present in the XML the input data is weighted and the results are stored in the relevant file" << endl;

	cout << endl;
	cout << "--saveFitXML" << endl;
	cout << "       Generates an XML roughly equivalent to the one you just input with final fit parameters and errors stored as fixed in the ParameterSet, good for projections" << endl;

	cout << endl;
	cout << "--generateToyXML" << endl;
	cout << "       Generates an XML roughly equivalent to the one you just input with final fit parameters stored in a format compatible with Toy Studies" << endl;

	cout << endl;
	cout << "--BuildConstraints" << endl;
	cout << "       Generate an XML which stores the parameters in a free XML result with constraints stored as an external Constraint XML segment" << endl;

	cout << endl;
	cout << "--saveAllToys" << endl;
	cout << "       When used in conjunction with a toy/MC/FC study this will cause all toys generated to be saved." << endl;
	cout << "       This creates a LOT of extra .root files from the toy study, but is a useful way to generate multiple toy datasets from one XML." << endl;

	cout << endl;
	cout << "--testIntegrator" << endl;
	cout << "       This allows you to test the Numerical vs Analytical Integrals from an XML" << endl;

	cout << endl;
	cout << "--helpProjections" << endl;
	cout << "       This will print a lot of options available for the Projections or ComponentProjections of a fit to data" << endl;

	cout << endl;
	cout << "--SendOutput <folder name>" << endl;
	cout << "	Write the output to a folder with the given name." << endl;

	cout << endl;

}

void ParseCommandLine::RapidFitAbout( string name )
{
	cout << endl <<"RapidFit is a fitter which is configured by a control XML to be provided at runtime" << endl << endl;

	cout << "RapidFit has many arguments, some which have higher priority than others." << endl;
	cout << endl <<"The order RapidFit will execute different arguments is:" << endl <<endl;

	cout << "\t 1)	save One Data Set or generate template XML and exit" <<endl;
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
	cout << "\t\t" << name << "\t --help" << endl << endl;
}

int ParseCommandLine::ParseThisCommandLine( RapidFitConfiguration& config, vector<string> argv )
{
	if( argv.size() == 1 )
	{
		cout << endl << "This is the Edinburgh, RapidFit Fitter!" << endl;
		cout << endl <<"for the priority of command line arguments run:" << endl <<endl<<"\t"<< argv[0] << " --about" << endl;
		cout << endl <<"run: " << argv[0] << " --help\t for more information" << endl << endl;
		return exit_RapidFit;
	}

	cout << "Current Runtime Arguments: " << endl;
	for( unsigned int i=0; i< argv.size(); ++i )
	{
		cout << argv[i] << " ";
	}
	cout << endl;

	//Parse command line arguments

	for( unsigned int argumentIndex = 1; argumentIndex < argv.size(); ++argumentIndex )
	{
		string currentArgument = string( argv[argumentIndex] );

		if( currentArgument == "--about" )
		{
			RapidFitAbout( argv[0] );
			return exit_RapidFit;
		}
		else if( currentArgument == "--help" )
		{
			RapidFitHelp();
			return exit_RapidFit;
		}
		else if( currentArgument == "--helpProjections" )
		{
			cout << endl << endl;
			cout << "Options available for Plotting Projections are:" << endl << endl;
			cout << ComponentPlotter::XML() << endl;
			return exit_RapidFit;
		}
		else if( currentArgument == "-f" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.configFileName = argv[argumentIndex];
				config.configFileNameFlag = true;
			}
			else
			{
				cerr << "Configuration file name not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "-p" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.parameterTemplates.push_back( argv[argumentIndex] );
				config.parameterTemplateFlag = true;
			}
			else
			{
				cerr << "Parameter template not supplied" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "-minimiser" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.theMinimiser = new MinimiserConfiguration( argv[argumentIndex] );
				config.theMinimiserFlag = true;
			}
			else
			{
				cerr << "Minimiser name not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "-function" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.theFunction = new FitFunctionConfiguration( argv[argumentIndex] );
				config.theFunctionFlag = true;
			}
			else
			{
				cerr << "Function name not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--saveOneDataSet" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.saveOneDataSetFileName = argv[argumentIndex];
				config.saveOneDataSetFlag = true;
			}
			else
			{
				cerr << "Data file name not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--saveOneFoamDataSet" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.saveOneDataSetFileName = argv[argumentIndex];
				config.saveOneFoamDataSetFlag = true;
			}
			else
			{
				cerr << "Data file name not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "-PDF" )
		{
			if( argumentIndex + 3 < argv.size() )
			{
				++argumentIndex;
				string pdfName = argv[argumentIndex];
				++argumentIndex;
				string dataSource = argv[argumentIndex];
				++argumentIndex;
				string phaseSpace = argv[argumentIndex];

				config.pdfsAndData.push_back( InputParsing::MakePDFWithData( pdfName, dataSource, phaseSpace ) );
			}
			else
			{
				cerr << "PDF not specified correctly: must be PDF name, data source and phase space prototype" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "-repeats" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.numberRepeats = atoi( argv[argumentIndex].c_str() );
				config.numberRepeatsFlag = true;
			}
			else
			{
				cerr << "Number of repeats not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--OverrideXML" )
		{
			if( argumentIndex + 2 < argv.size() )
			{
				++argumentIndex;
				string path = string( argv[argumentIndex] );
				++argumentIndex;
				string value = string( argv[argumentIndex] );
				pair<string,string> temp_pair;
				temp_pair.first = path;
				temp_pair.second = value;
				cout << endl;
				cout << "ATTENTION!!ATTENTION!!ATTENTION!!ATTENTION!!ATTENTION!!ATTENTION!!ATTENTION" << endl;
				cout << "Overriding an XMLTag with:" << endl;
				cout << "\t" << path << "\t" << value << endl;
				cout << "NB:\t this is an expert option for complex situations which haven't been forseen in the design yet" << endl;
				cout << "ATTENTION!!ATTENTION!!ATTENTION!!ATTENTION!!ATTENTION!!ATTENTION!!ATTENTION" << endl << endl;
				config.XMLOverrideList->push_back( temp_pair );
			}
			else
			{
				cerr << "XML Override not correctly defined!" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--PhysParam" )
		{
			if( argumentIndex +1 < argv.size() )
			{
				++argumentIndex;
				string PhysParam = argv[argumentIndex];
				config.CommandLineParamvector.push_back( PhysParam );
			}
			else
			{
				cerr << "PhysParam not correctly defined at Runtime" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--testComponentPlot" )
		{
			config.testComponentPlotFlag = true;
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.plotFileName = argv[argumentIndex];
			}
			else
			{
				cerr << "Path to store plot file not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--observableName" )
		{
			config.observableNameFlag = true;
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.observableName = argv[argumentIndex];
			}
			else
			{
				cerr << "Observable Name not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--doPlotting" )
		{
			config.doPlottingFlag = true;
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.plotFileName = argv[argumentIndex];
			}
			else
			{
				cerr << "Path to store plot file not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--doPulls" )
		{
			config.doPullsFlag = true;
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.pullFileName = argv[argumentIndex];
			}
			else
			{
				cerr << "Path to store pull plot file not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--defineScan" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.defineScanFlag = true;
				config.Scan_X.push_back( argv[argumentIndex] );
			}
			else
			{
				cerr << "Scan Not Correctly Formatted" <<endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--defineContour" )
		{
			if( argumentIndex + 2 < argv.size() )
			{
				config.defineContourFlag = true;
				++argumentIndex;
				config.Contour_X.push_back( argv[argumentIndex] );
				++argumentIndex;
				config.Contour_Y.push_back( argv[argumentIndex] );
			}
			else
			{
				cerr << "Contour Not Correctly Formatted" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--NuisenceModel" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.Nuisencemodel = (unsigned)atoi( argv[argumentIndex].c_str() );
			}
			else
			{
				cerr << "Not correctly Formatted Nuisence Input" <<endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--SetSeed" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.RuntimeSeed.push_back( atoi( argv[argumentIndex].c_str() ) );
			}
			else
			{
				cerr << "Seed Not Correctly Defined at Runtime" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--MCStepSize" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.MCStepSize = ( argv[argumentIndex] );
			}
			else
			{
				cerr << "Badly Defined Step size" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--MCStartEntry" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.MCStartEntry = ( argv[argumentIndex] );
			}
			else
			{
				cerr << "Badly Defined Start Entry" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--OutputLevel" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.OutputLevel = atoi( argv[argumentIndex].c_str() );
				config.OutputLevelSet=true;
			}
			else
			{
				cerr << "Badly Defined Output Level" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--OutputLevelScans" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				config.OutputLevel2 = atoi( argv[argumentIndex].c_str() );
				config.OutputLevelSet=true;
			}
			else
			{
				cerr << "Badly Defined Output Level" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--GOF" )
		{
			config.GOF_Flag = true;
			if( argumentIndex + 2 < argv.size() )
			{
				++argumentIndex;
				config.jobNum = atoi( argv[argumentIndex].c_str() );
				++argumentIndex;
				config.nData  = atoi( argv[argumentIndex].c_str() );
			}
			else
			{
				cerr << "Job number for GOF not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--JackKnife" )
		{
			config.JackKnife_Flag = true;
			if( argumentIndex + 2 < argv.size() )
			{
				++argumentIndex;
				config.jackStartNum = atoi( argv[argumentIndex].c_str() );
				++argumentIndex;
				config.jackStopNum  = atoi( argv[argumentIndex].c_str() );
			}
			else
			{
				cerr << "Start and stop event numbers for JackKnife not specified" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--MakeTemplate" )
		{
			config.makeTemplateXML = true;
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				for( ; argumentIndex < argv.size(); ++argumentIndex )
				{
					config.templatePDFs.push_back( argv[argumentIndex] );
				}
			}
			else
			{
				cerr << "Required to give at least 1 PDF as input fot an XML template" << endl;
				return BAD_COMMAND_LINE_ARG;
			}
		}
		else if( currentArgument == "--DebugClass" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				vector<string> class_names = StringProcessing::SplitString( argv[argumentIndex], ':' );
				DebugClass* thisDebug = new DebugClass( false );
				thisDebug->SetClassNames( class_names );

				if( config.debug != NULL ) delete config.debug;
				config.debug = thisDebug;
			}
		}
		else if( currentArgument == "--SendOutput" )
		{
			if( argumentIndex + 1 < argv.size() )
			{
				++argumentIndex;
				ResultFormatter::SetOutputFolder(argv[argumentIndex]);	
			}
			else{
				cerr << "Required to give output folder name" << endl;
				return BAD_COMMAND_LINE_ARG;
			}	
		}

		//	The Parameters beyond here are for setting boolean flags
		else if( currentArgument == "--testIntegrator" )			{	config.testIntegratorFlag = true;			}
		else if( currentArgument == "--testRapidIntegrator" )			{	config.testRapidIntegratorFlag = true;			}
		else if( currentArgument == "--calculateFitFractions" )			{	config.calculateFitFractionsFlag = true;		}
		else if( currentArgument == "--calculateAcceptanceWeights" )		{	config.calculateAcceptanceWeights = true;		}
		else if( currentArgument == "--calculateAcceptanceCoefficients" )		{	config.calculateAcceptanceCoefficients = true;		}
		else if( currentArgument == "--calculateAcceptanceWeightsWithSwave" )	{	config.calculateAcceptanceWeightsWithSwave = true;	}
		else if( currentArgument == "--calculatePerEventAcceptance" )		{	config.calculatePerEventAcceptance = true;		}
		else if( currentArgument == "--doLLscan" )				{	config.doLLscanFlag = true;				}
		else if( currentArgument == "--doLLcontour" )				{	config.doLLcontourFlag = true;				}
		else if( currentArgument == "--doFCscan" )				{	config.doFC_Flag = true;				}
		else if( currentArgument == "--useUUID" )				{	config.UUID_Flag = true;				}
		else if( currentArgument == "--BurnToROOT" )				{	config.BurnToROOTFlag = true;				}
		else if( currentArgument == "--MCStudy" )				{	config.MCStudyFlag = true;				}
		else if( currentArgument == "--ForceContinue" )				{	config.Force_Continue_Flag = true;			}
		else if( currentArgument == "--DontStartAtCenter" )			{	config.StartAtCenterFlag = false;			}
		else if( currentArgument == "--WeightDataSet" ) 			{       config.WeightDataSet=true;				}
		else if( currentArgument == "--saveFitXML" )				{	config.saveFitXML = true;				}
		else if( currentArgument == "--generateToyXML" )			{	config.generateToyXML = true;				}
		else if( currentArgument == "--fixedTotalToys" )			{	config.fixedTotalToys = true;				}
		else if( currentArgument == "--saveAllToys" )				{	config.saveAllToys = true;				}
		else if( currentArgument == "--Debug" )	{  if( config.debug != NULL ){ delete config.debug;}	config.debug = new DebugClass(true);	}
		else if( currentArgument == "--BuildConstraints" )			{	config.BuildConstraints = true;				}
		else if( currentArgument == "--disableLatexOutput" )			{	config.disableLatexOutput = true;			}

		//	We didn't understand the argument to end up here
		else
		{
			cerr << "Unrecognised argument: " << currentArgument << endl;
			exit(-9834);
			return BAD_COMMAND_LINE_ARG;
		}
	}

	return 0;
}
