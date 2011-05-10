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
#include "ToyStudyResult.h"
#include "LLscanResult.h"
#include "LLscanResult2D.h"
#include "main.h"
#include "DataSetConfiguration.h"
#include "IDataSet.h"
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

#ifndef __CINT__
int RapidFit( int argc, char * argv[] );

int main( int argc, char * argv[] )
{
	return RapidFit( argc, argv );
}
#endif

int RapidFit( int argc, char * argv[] )
{
	//ProfilerStart("profile.out");
	//TProof * proof = TProof::Open( "" ); // idea of using PROOF for paralellisation
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

	if ( argc == 1 )
	{
		//Default behaviour if no command line arguments
		ToyStudy * testStudy = new ToyStudy("../config/MC09Toy_noBkg.xml");
		ToyStudyResult * fitResults = testStudy->DoWholeStudy( true );
		ResultFormatter::MakePullPlots( "SeparateParameter", "pullPlots.root", fitResults );
		ResultFormatter::LatexOutputFitResult( fitResults->GetFitResult(0) );
		delete testStudy;
	}
	else
	{
		//Variables to store command line arguments
		int numberRepeats = 0;
		int Nuisencemodel=2;
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
		vector<pair<string, string> >* XMLOverideList = new vector<pair<string,string> >;

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

		//Parse command line arguments
		string currentArgument;
		for ( int argumentIndex = 1; argumentIndex < argc; ++argumentIndex )
		{
			currentArgument = argv[argumentIndex];

			if ( currentArgument == "--help" )
			{
					++argumentIndex;

				cout << "QUICK HELP FOR RAPIDFIT" << endl ;

				cout << endl ;
				cout << " -f <filename>   " << endl ;
				cout << "      Specifies the XML file to drive RapidFit " <<endl ;

				cout << endl ;
				cout << " -repeats n   " << endl ;
				cout << "      Specifies the number of repeats for a Toy study " <<endl ;

				cout << endl ;
				cout << " --doLLscan  " << endl;
				cout << "      Causes a set of LL scans to be perfomed around the fit minimum" << endl ;
				cout << "      The parameters which will be scanned must be specified in the xml file as " << endl ;
				cout << "      The range may be defined as hard numerical limits as shown below or as a " <<endl;
				cout << "      Range of n sigma from the fit minima, as shown in Y_Param for LLcontour" << endl;
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
				cout << "      Causes a set of LL scans to be perfomed around the fit minimum" << endl ;
				cout << "      The parameters which will be scanned must be specified in the xml file as " << endl ;
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
				cout << "      Causes one Toy dataset to be written out to a file " <<endl ;
				cout << "      filename = file to write the dataset to" << endl ;

				cout << endl ;
				cout << " --testIntegrator   " << endl ;
				cout << "      Usefl feature which only tests the numerical<=>analytic integrator for each PDF then exits " <<endl ;

				cout << endl;
				cout << " --SetSeed 12345" << endl;
				cout << "      Set the Random seed to 12345 if you wish to make the output reproducable. Useful on Batch Systems" << endl;

				return 1;

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
					return 1;
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
					return 1;
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
					return 1;
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
					return 1;
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
					return 1;
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
					return 1;
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
					return 1;
				}
			}
                        else if ( currentArgument == "--OverideXML" )
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
					cout << "Overiding an XMLTag with:" << endl;
					cout << "\t" << path << "\t" << value << endl;
					cout << "WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!WARNING" << endl << endl;
					XMLOverideList->push_back( temp_pair );
				}
				else
				{
					cerr << "XML Overide not correctly defined!" << endl;
					return 1;
				}
			}
			else if ( currentArgument == "--testIntegrator" )
			{
				testIntegratorFlag = true;
			}
			else if ( currentArgument == "--testRapidIntegrator" )
			{
				testRapidIntegratorFlag = true;
			}
			else if ( currentArgument == "--calculateAcceptanceWeights" )
			{
				calculateAcceptanceWeights = true;
			}
			else if ( currentArgument == "--calculateAcceptanceWeightsWithSwave" )
                        {
                                calculateAcceptanceWeightsWithSwave = true;
                        }
			else if ( currentArgument == "--calculatePerEventAcceptance" )
			{
				calculatePerEventAcceptance = true;
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
					return 1;
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
					return 1;
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
					return 1;
				}
			}
			else if ( currentArgument == "--doLLscan" )
			{
				doLLscanFlag = true;

			}
			else if ( currentArgument == "--doLLcontour" )
			{
				doLLcontourFlag = true;
			}
			else if ( currentArgument == "--doFCscan" )
			{
				doFC_Flag = true;
			}
			else if ( currentArgument == "--useUUID" )
			{
				UUID_Flag = true;
			}
			else if ( currentArgument == "--debugFC" )
			{
				FC_Debug_Flag = true;
			}
			else if ( currentArgument == "--defineScan" )
			{
				if ( argumentIndex + 1 < argc )
				{
					defineScanFlag = true;
					Scan_X.push_back( argv[argumentIndex] );
				} else {
					cerr << "Scan Not Correctly Formatted" <<endl;
					return 1;
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
					return 1;
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
					return 1;
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
					return 1;
				}
			}
			else if ( currentArgument == "--BurnToROOT" )
			{
				BurnToROOTFlag = true;
			}
			else
			{
				cerr << "Unrecognised argument: " << currentArgument << endl;
				exit(1);
			}
		}

		//Load a config file
		XMLConfigReader* xmlFile=NULL;
		if (configFileNameFlag)
		{
			xmlFile = new XMLConfigReader(configFileName, XMLOverideList);
			if ( xmlFile->IsLoaded() )
			{
				cout << "XML config file " << configFileName << " loaded" << endl;
			}
			else
			{
				cerr << "XML config file " << configFileName << " not found" << endl;
				return 1;
			}
		}
		if( ( !RuntimeSeed.empty() && xmlFile->IsLoaded() ) && !UUID_Flag )
		{
			cout << "Setting Seed At Runtime to be: " << RuntimeSeed[0] << endl;
			xmlFile->SetSeed( unsigned(RuntimeSeed[0]) );
		} else if ( !RuntimeSeed.empty() && UUID_Flag && xmlFile->IsLoaded() )
		{
			//  I have used the TRandom3 code as inspiration for 'Salting the Seed' of the UUID at runtime
			cout << "Using a Random seed of UUID x RuntimeSeed to be 'unique' for all machines"<<endl;
			TUUID uid;
			UChar_t uuid[16];
			uid.GetUUID(uuid);
			RuntimeSeed[0] = RuntimeSeed[0] * ( uuid[ 2*(RuntimeSeed[0]%8) ]*256 +uuid[ 2*(RuntimeSeed[0]%8) ] );
			xmlFile->SetSeed( unsigned(RuntimeSeed[0]) );
		}

		//Create a parameter set
		vector<ParameterSet*> argumentParameterSet;
		if (parameterTemplateFlag)
		{
			argumentParameterSet.push_back( InputParsing::MakeParameterSet(parameterTemplates) );
		}

		//Choose what action to take, now that you've configured everything
		if (saveOneDataSetFlag)
		{
			if (configFileNameFlag)
			{
				//Make a file containing toy data from the PDF
				vector<PDFWithData*> quickDataGen = xmlFile->GetPDFsAndData();
				vector<IDataSet*> quickDataSet;
				for( unsigned int i=0; i< quickDataGen.size(); ++i )
				{
					quickDataGen[i]->SetPhysicsParameters( xmlFile->GetFitParameters() );
					IDataSet* temp_dataSet = quickDataGen[i]->GetDataSet();
					quickDataSet.push_back( temp_dataSet );
				}

				if( quickDataGen.size() != quickDataSet.size() )
				{	
					cerr << "Unexpected, Unexplained, Error... Goodbye!" << endl;
					exit(-53);
				}

				ResultFormatter::MakeRootDataFile( saveOneDataSetFileName, quickDataSet );

				while( !quickDataGen.empty() )
				{
					delete quickDataGen.back();
					delete quickDataSet.back();
					quickDataGen.pop_back();
					quickDataSet.pop_back();
				}
			}
			else
			{
				cerr << "No data set specified" << endl;
				return 1;
			}
		}
		else if (testIntegratorFlag)
		{
			if (configFileNameFlag)
			{
				//Compare numerical and analytical integration
				PDFWithData * quickData = xmlFile->GetPDFsAndData()[0];
				quickData->SetPhysicsParameters( xmlFile->GetFitParameters() );
				IDataSet * quickDataSet = quickData->GetDataSet();
				RapidFitIntegrator * testIntegrator = new RapidFitIntegrator( quickData->GetPDF() );
				testIntegrator->Integral( quickDataSet->GetDataPoint(0), quickDataSet->GetBoundary() );
				delete testIntegrator;
				delete quickData;
			}
			else
			{
				cerr << "No data set specified" << endl;
				return 1;
			}
		}
    else if (calculateAcceptanceWeights)
    {
      if (configFileNameFlag)
      {
        // Calculate the acceptance weights from MC
        PDFWithData * pdfAndData = xmlFile->GetPDFsAndData()[0];
        pdfAndData->SetPhysicsParameters( xmlFile->GetFitParameters() );
        IDataSet * dataSet = pdfAndData->GetDataSet();
        IPDF * pdf = pdfAndData->GetPDF();
        Mathematics::calculateAcceptanceWeights(dataSet, pdf);
	delete pdfAndData;
      }
      else
      {
        cerr << "No data set specified" << endl;
        return 1;
      }
    }
	else if (calculateAcceptanceWeightsWithSwave)
    {
      if (configFileNameFlag)
      {
        // Calculate the acceptance weights from MC
        PDFWithData * pdfAndData = xmlFile->GetPDFsAndData()[0];
        pdfAndData->SetPhysicsParameters( xmlFile->GetFitParameters() );
        IDataSet * dataSet = pdfAndData->GetDataSet();
        IPDF * pdf = pdfAndData->GetPDF();
        Mathematics::calculateAcceptanceWeightsWithSwave(dataSet, pdf);
	delete pdfAndData;
      }
      else
      {
        cerr << "No data set specified" << endl;
        return 1;
      }
    }


	  else if (calculatePerEventAcceptance)
    {
      PerEventAngularAcceptance* a = new PerEventAngularAcceptance("jpsikmc09_loose.root","Bu2JpsiKTuple/DecayTree", "out2.root");
      for (int iter = 1; iter <= 3; ++iter)
      {
        a->fillEffHistos( iter );
        a->loopOnReconstructedBs();
      }
      a->writeHistos();
      delete a;
    }
		else if (testPlotFlag)
		{
			if (configFileNameFlag)
			{
				//Project the PDF onto the data
				PDFWithData * quickData = xmlFile->GetPDFsAndData()[0];
				quickData->SetPhysicsParameters( xmlFile->GetFitParameters() );
				IDataSet * quickDataSet = quickData->GetDataSet();
				Plotter * testPlotter = new Plotter( quickData->GetPDF(), quickDataSet );
				testPlotter->PlotAllObservables(plotFileName);
				delete testPlotter;
				delete quickData;
			}
			else
			{
				cerr << "No data set specified" << endl;
				return 1;
			}
		}
		else if (configFileNameFlag)
		{
			//Actually configure a fit: first configure the output
			OutputConfiguration * makeOutput = xmlFile->GetOutputConfiguration();

			//Command line argments override the config file
			if (!theMinimiserFlag)
			{
				theMinimiser = xmlFile->GetMinimiserConfiguration();
			}
			if ( defineContourFlag )
			{
				for( unsigned short int i=0; i < Contour_X.size(); ++i)
				{
					makeOutput->Clear2DScanList();
					makeOutput->AddContour( Contour_X[i], Contour_Y[i] );
				}
			}
			if ( defineScanFlag )
			{
				for( unsigned short int i=0; i < Scan_X.size(); ++i )
				{
					makeOutput->AddScan( Scan_X[i] );
				}
			}
			if (!theFunctionFlag)
			{
				theFunction = xmlFile->GetFitFunctionConfiguration();
				// If weights were specified then we need to let the output plotting know
				if( theFunction->GetWeightsWereUsed() ) makeOutput->SetWeightsWereUsed( theFunction->GetWeightName() ) ;
			}
			if (!parameterTemplateFlag)
			{
				argumentParameterSet = xmlFile->GetFitParameters();
			}
			if (!numberRepeatsFlag)
			{
				numberRepeats = xmlFile->GetNumberRepeats();
			}
			if ( pdfsAndData.size() == 0 )
			{
				pdfsAndData = xmlFile->GetPDFsAndData();

				//	If we are performing a scan we want to check for Data Generation instances and generate/store the data in a cache for future use
				if( doLLscanFlag || ( doLLcontourFlag || doFC_Flag ) )
				{
					for( unsigned int pdf_num=0; pdf_num< pdfsAndData.size(); ++pdf_num )
					{
						vector<DataSetConfiguration*> DataConfigs = pdfsAndData[pdf_num]->GetAllDataSetConfigs();
						for( unsigned int config_num=0; config_num < DataConfigs.size(); ++config_num )
						{
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

			//Pulls
			if ( doPullsFlag )
			{
				makeOutput->SetPullFileName(pullFileName);
			}
			else
			{
				doPullsFlag = makeOutput->DoPullPlots();
			}

			//Plotting
			if ( doPlottingFlag )
			{
				makeOutput->MakeAllPlots(plotFileName);
			}

			//Pick a toy study if there are repeats, or if pull plots are wanted
			if ( ( ( numberRepeats > 1 ) && ((!doFC_Flag) && (!doLLcontourFlag) && (!doLLscanFlag) ) ) || doPullsFlag )
			{
				vector< ConstraintFunction* > XMLConstraints = xmlFile->GetConstraints();
				//Do the toy study
				ToyStudy newStudy( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, XMLConstraints, numberRepeats );
				ToyStudyResult * fitResults = newStudy.DoWholeStudy( true );

				//Output results
				makeOutput->OutputToyResult(fitResults);
				makeOutput->OutputFitResult( fitResults->GetFitResult(0) );

				while( !XMLConstraints.empty() )
				{
					delete XMLConstraints.back();
					XMLConstraints.pop_back();
				}
			}
			else
			{
				//		This is re-used for FC scans and forms FC Step 1
				cout << "\n\n\t\tStarting Fit to Find Global Minima!\n"<<endl;
				//Do the fit to find GLOBAL MINIMA
				ToyStudyResult* GlobalFitResult = new ToyStudyResult( argumentParameterSet.back()->GetAllNames() );
				GlobalFitResult->StartStopwatch();

				vector< ConstraintFunction* > XMLConstraints = xmlFile->GetConstraints();

				FitResult * GlobalResult = FitAssembler::DoSafeFit( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, XMLConstraints );

				GlobalFitResult->AddFitResult( GlobalResult );

				//ResultFormatter::FlatNTuplePullPlots( string("Global_Fit.root"), GlobalFitResult );

				cout << "\n\n\t\tFit Output:" <<endl;
				//Output results
				makeOutput->SetInputResults( GlobalResult->GetResultParameterSet() );
				makeOutput->OutputFitResult( GlobalResult );
				ResultFormatter::ReviewOutput( GlobalResult );
				if( BurnToROOTFlag )
				{
					ResultFormatter::WriteFlatNtuple( string( "Global_Fit_Result.root" ), GlobalFitResult );
				}

				//Do LL scan
				if( doLLscanFlag ) {
					//  Temporary Holder for the scan outputs to plot
					LLscanResult * llResult;
					vector<ToyStudyResult*> scanSoloResult;
					//  Store
					vector<LLscanResult*> scanResults;
					vector<string> LLscanList = makeOutput->GetScanList();

					for(unsigned int ii=0; ii < LLscanList.size() ; ++ii)
					{
						ToyStudyResult* scan_result = FitAssembler::SingleScan( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints(), makeOutput, LLscanList[ii] );
						llResult = ResultFormatter::LLScan( scan_result, LLscanList[ii] );
						scanResults.push_back( llResult );
						scanSoloResult.push_back( scan_result );
					}
					makeOutput->SetLLscanFileName( LLscanFileName );
					makeOutput->OutputLLscanResult( scanResults ) ;
					for(unsigned int ii=0; ii < LLscanList.size(); ++ii )
					{
						TString output_scan_dat("LLScanData");
						output_scan_dat.Append(LLscanList[ii]);
						output_scan_dat.Append(".root" );
						//ResultFormatter::FlatNTuplePullPlots( string( output_scan_dat ), scanSoloResult[ii] );
						ResultFormatter::WriteFlatNtuple( string( output_scan_dat ), scanSoloResult[ii] );
					}
				}

				//		This is re-used for FC scans and forms FC Step 2
				ToyStudyResult* _2DResultForFC=NULL;
				//Do 2D LL scan for deltaGamma and phis
				if( doLLcontourFlag || doFC_Flag ) {

					//  Soon to be gladly deprecated!
					LLscanResult2D * llContourResult ;
					vector<LLscanResult2D*> contourResults ;

					//  Array of individual Results
					vector<ToyStudyResult*> SoloContourResults;

					vector<pair<string, string> > _2DLLscanList = makeOutput->Get2DScanList();

					unsigned int initial_scan=0;
					if( ( ( _2DLLscanList.size() > 1 ) && doFC_Flag ) || defineContourFlag )
					{
						cerr << "\n\nPERFORMING ONLY ONE 2D SCAN, CHECK THIS IS EXPECTED!" <<endl;
						initial_scan = int(_2DLLscanList.size()-1);
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
						vector<ToyStudyResult*> Temp_Results = FitAssembler::ContourScan( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints(), makeOutput, name1, name2 );

						if( !doFC_Flag )
						{
							llContourResult = ResultFormatter::LLScan2D( Temp_Results, name1, name2 );
							contourResults.push_back( llContourResult );
						}

						GlobalResult->GetResultParameterSet()->GetResultParameter( name1 )->ForcePullValue( -9999 );
						GlobalResult->GetResultParameterSet()->GetResultParameter( name1 )->ForceOriginalValue( -9999 );
						GlobalResult->GetResultParameterSet()->GetResultParameter( name2 )->ForcePullValue( -9999 );
						GlobalResult->GetResultParameterSet()->GetResultParameter( name2 )->ForceOriginalValue( -9999 );

						TString ext_string("-");
						ext_string.Append( name1 );
						ext_string.Append("_");
						ext_string.Append( name2 );
						ext_string.Append( ".root" );
						TString TempName = LLcontourFileName;
						TempName.Append( ext_string );
						LLcontourFileNamez.push_back( TempName );

						//	Linearize the data and then Set the generate values to a default
						SoloContourResults.push_back( new ToyStudyResult( Temp_Results ) );
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

					//  Don't output the scan files when doing a FC scan (This is only noise on Condor btw!)
					if( doFC_Flag )  {
						_2DResultForFC = SoloContourResults.back();
					} else  {
						for(unsigned int ii=0; ii < _2DLLscanList.size(); ++ii )
						{
							//	STOP RELYING ON THIS OUTPUT, IT'S BADLY WRITTEN, AND IT WILL BE DISABLED SOONER RATHER THAN LATER
							makeOutput->SetLLcontourFileName( LLcontourFileNamez[ii].Data() );
							makeOutput->OutputLLcontourResult( contourResults ) ;

							//	This output is A LOT safer as it intended to be fed to a sperate plotting program.
							//	Formatting SHOULD be seperate to generation!
							string output_scan_dat("LLcontourScanData.root");
							if( ii > 0 )
							{
								TString _scan_dat("LLcontourScanData_");
								_scan_dat+= ii+1;
								_scan_dat.Append(".root");
								output_scan_dat=_scan_dat;
							}
							//	Add the Global Results and 'Linearize' the output
							vector<ToyStudyResult*> TempContourResults;
							TempContourResults.push_back( GlobalFitResult );
							TempContourResults.push_back( SoloContourResults[ii] );
							ToyStudyResult* TempContourResults2 = new ToyStudyResult( TempContourResults );
							ResultFormatter::WriteFlatNtuple( output_scan_dat , TempContourResults2 );
						}
					}

					while( !SoloContourResults.empty() )
					{
						delete SoloContourResults.back();
						SoloContourResults.pop_back();
					}
				}
				//	Do the main work of the FC scan
				if( doFC_Flag )
				{
					vector<unsigned int> numberRepeatsVec;
					if( numberRepeatsFlag ) numberRepeatsVec.push_back( unsigned(numberRepeats) );
					
					//	Do FC scan
					ToyStudyResult* AllFCResults = FitAssembler::FeldmanCousins( GlobalFitResult, _2DResultForFC, numberRepeatsVec, unsigned(int(Nuisencemodel)), FC_Debug_Flag, makeOutput, theMinimiser, theFunction,  xmlFile, pdfsAndData );

					//		STORE THE OUTPUT OF THE TOY STUDIES
					ResultFormatter::WriteFlatNtuple( "FCOutput.root", AllFCResults );
				}
			
				while( !XMLConstraints.empty() )
				{
					delete XMLConstraints.back();
					XMLConstraints.pop_back();
				}
				delete GlobalResult;
				delete GlobalFitResult;
			}
		}
		else
		{
			//Default action - presumably a fit or a toy study
			cerr << "No action performed" << endl;
			return 1;
		}

	//	Clean UP!
	while ( !pdfsAndData.empty() )
	{
		delete pdfsAndData.back();
		pdfsAndData.pop_back();
	}
	delete xmlFile;

	}

	//Exit blurb
	time(&timeNow);
	cout << endl << "RapidFit" << endl;
	cout << "Ending time: " << ctime( &timeNow ) << endl;
	//ProfilerStop();

	return 0;
}
