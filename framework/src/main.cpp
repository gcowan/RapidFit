/**
  @file main.cpp

  Entry point for RapidFit

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

#include <string>
#include <vector>
#include <iostream>
#include <ctime>
#include "Mathematics.h"
#include "FitAssembler.h"
#include "ToyStudy.h"
#include "XMLConfigReader.h"
#include "ResultFormatter.h"
#include "InputParsing.h"
#include "RapidFitIntegrator.h"
#include "Plotter.h"
#include <stdio.h>
#include <stdlib.h>
#include "MakeFoam.h"
#include "PerEventAngularAcceptance.h"
#include "OutputConfiguration.h"
//#include "ScanParam.h"
#include "TString.h"
#include "ToyStudyResult.h"
#include "LLscanResult.h"
#include "LLscanResult2D.h"

using namespace std;

int main( int argc, char * argv[] )
{

	//TProof * proof = TProof::Open( "" ); // idea of using PROOF for paralellisation
	//Welcome blurb
	time_t timeNow;
	time(&timeNow);
	cout << endl << "RapidFit" << endl;
	cout << "Starting time: " << ctime( &timeNow ) << endl;

	if ( argc == 1 )
	{
		//Default behaviour if no command line arguments
		ToyStudy * testStudy = new ToyStudy("../config/MC09Toy_noBkg.xml");
		ToyStudyResult * fitResults = testStudy->DoWholeStudy();
		ResultFormatter::MakePullPlots( "SeparateParameter", "pullPlots.root", fitResults );
		ResultFormatter::LatexOutputFitResult( fitResults->GetFitResult(0) );
	}
	else
	{
		//Variables to store command line arguments
		int numberRepeats = 0;
		string configFileName = "";
		vector<string> parameterTemplates;
		MinimiserConfiguration * theMinimiser;
		FitFunctionConfiguration * theFunction;
		string saveOneDataSetFileName = "";
		string plotFileName = "FitPlots.root";
		string pullFileName = "PullPlots.root";
		string LLscanFileName = "LLscanPlots.root";
		string LLcontourFileName = "LLcontourPlots";
		vector< PDFWithData* > pdfsAndData;
		int numberLLscanPoints = 10 ;
		double LLscanRange = 2 ;
		vector<int> RuntimeSeed;

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

		//Parse command line arguments
		string currentArgument;
		for ( int argumentIndex = 1; argumentIndex < argc; argumentIndex++ )
		{
			currentArgument = argv[argumentIndex];

			if ( currentArgument == "--help" )
			{
					argumentIndex++;

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
				cout << "				<Points>number of points in y axis</Points> " << endl ;
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
					argumentIndex++;
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
					argumentIndex++;
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
					argumentIndex++;
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
					argumentIndex++;
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
					argumentIndex++;
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
					argumentIndex++;
					string pdfName = argv[argumentIndex];
					argumentIndex++;
					string dataSource = argv[argumentIndex];
					argumentIndex++;
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
					argumentIndex++;
					numberRepeats = atoi( argv[argumentIndex] );
					numberRepeatsFlag = true;
				}
				else
				{
					cerr << "Number of repeats not specified" << endl;
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
					argumentIndex++;
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
					argumentIndex++;
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
					argumentIndex++;
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

//			else if ( currentArgument == "-numberLLscanPoints" )
//			{
//				if ( argumentIndex + 1 < argc )
//				{
//					argumentIndex++;
//					numberLLscanPoints = atoi( argv[argumentIndex] );
//					if( numberLLscanPoints <= 0 ) { cerr << "Number of ll scan points not >=0" << endl; return 1; }
//				}
//				else
//				{
//					cerr << "Number of ll scan points not specified" << endl;
//					return 1;
//				}
//			}
			else if ( currentArgument == "--doLLcontour" )
			{
				doLLcontourFlag = true;
			}

			else if ( currentArgument == "--SetSeed" )
			{
				if ( argumentIndex + 1 < argc )
				{
					argumentIndex++;
					RuntimeSeed.push_back( atoi( argv[argumentIndex] ) );
				} else {
					cerr << "Seed Not Correctly Defined at Runtime" << endl;
					return 1;
				}
			}
			else
			{
				cerr << "Unrecognised argument: " << currentArgument << endl;
				exit(1);
			}
		}

		//Load a config file
		XMLConfigReader * xmlFile;
		if (configFileNameFlag)
		{
			xmlFile = new XMLConfigReader(configFileName);
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
		if( !RuntimeSeed.empty() && xmlFile->IsLoaded() )
		{
			cout << "Setting Seed At Runtime to be: " << RuntimeSeed[0] << endl;
			xmlFile->SetSeed( RuntimeSeed[0] );
		}

		//Create a parameter set
		ParameterSet * argumentParameterSet;
		if (parameterTemplateFlag)
		{
			argumentParameterSet = InputParsing::MakeParameterSet(parameterTemplates);
		}

		//Choose what action to take, now that you've configured everything
		if (saveOneDataSetFlag)
		{
			if (configFileNameFlag)
			{
				//Make a file containing toy data from the PDF
				PDFWithData * quickData = xmlFile->GetPDFsAndData()[0];
				quickData->SetPhysicsParameters( xmlFile->GetFitParameters() );
				ResultFormatter::MakeRootDataFile( saveOneDataSetFileName, quickData->GetDataSet() );
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
      }
      else
      {
        cerr << "No data set specified" << endl;
        return 1;
      }
    }


	  else if (calculatePerEventAcceptance)
    {
      PerEventAngularAcceptance a = PerEventAngularAcceptance("jpsikmc09_loose.root","Bu2JpsiKTuple/DecayTree", "out2.root");
      for (int iter = 1; iter <= 3; iter++)
      {
        a.fillEffHistos( iter );
        a.loopOnReconstructedBs();
      }
      a.writeHistos();
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
			if ( numberRepeats > 1 || doPullsFlag )
			{
				//Do the toy study
				ToyStudy newStudy( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints(), numberRepeats );
				ToyStudyResult * fitResults = newStudy.DoWholeStudy();

				//Output results
				makeOutput->OutputToyResult(fitResults);
				makeOutput->OutputFitResult( fitResults->GetFitResult(0) );

			}
			else
			{
				//Do the fit
				FitResult * oneResult = FitAssembler::DoFit( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints() );

				//Output results
				makeOutput->SetInputResults( oneResult->GetResultParameterSet() );
				makeOutput->OutputFitResult( oneResult );


				//Do LL scan
				if( doLLscanFlag ) {
					//  Temporary Holder for the scan outputs to plot
					LLscanResult * llResult;
					vector<ToyStudyResult*> scanSoloResult;
					//  Store
					vector<LLscanResult*> scanResults;
					vector<string> LLscanList = makeOutput->GetScanList();

					for(int ii=0; ii < LLscanList.size() ; ii++)
					{
						ToyStudyResult* scan_result = FitAssembler::SingleScan( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints(), makeOutput, LLscanList[ii] );
						llResult = ResultFormatter::LLScan( scan_result, LLscanList[ii] );
						scanResults.push_back( llResult );
						scanSoloResult.push_back( scan_result );
					}
					makeOutput->SetLLscanFileName( LLscanFileName );
					makeOutput->OutputLLscanResult( scanResults ) ;
					for( int ii=0; ii < LLscanList.size(); ii++ )
					{
						TString output_scan_dat("LLScanData");
						output_scan_dat.Append(LLscanList[ii]);
						output_scan_dat.Append(".root" );
						ResultFormatter::FlatNTuplePullPlots( string( output_scan_dat ), scanSoloResult[ii] );
					}
				}

				//Do 2D LL scan for deltaGamma and phis
				if( doLLcontourFlag ) {
					LLscanResult2D * llContourResult ;
					vector<LLscanResult2D*> contourResults ;

					vector<ToyStudyResult*> SoloContourResults;
					vector<ToyStudyResult*> AllContourFitResults;
					ToyStudyResult* ContourLinearResults;

					vector<pair<string, string> > _2DLLscanList = makeOutput->Get2DScanList();

					for( int ii=0; ii < _2DLLscanList.size() ; ii++ )
					{
						string name1 = _2DLLscanList[ii].first;
						string name2 = _2DLLscanList[ii].second;
						SoloContourResults = FitAssembler::ContourScan( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints(), makeOutput, name1, name2 );

						llContourResult = ResultFormatter::LLScan2D( SoloContourResults, name1, name2 );
						contourResults.push_back( llContourResult );

						TString ext_string("-");
						ext_string.Append( name1 );
						ext_string.Append("_");
						ext_string.Append( name2 );
						ext_string.Append( ".root" );
						LLcontourFileName.append( ext_string );

						ContourLinearResults = new ToyStudyResult( SoloContourResults );

						AllContourFitResults.push_back( ContourLinearResults );

					}
					makeOutput->SetLLcontourFileName( LLcontourFileName );
					makeOutput->OutputLLcontourResult( contourResults ) ;
					for( int ii=0; ii < _2DLLscanList.size(); ii++ )
					{
						TString output_scan_dat( "LLcontourScanData");
						TString ext("_");
						ext.Append(_2DLLscanList[ii].first);
						ext.Append("_");
						ext.Append(_2DLLscanList[ii].second);
						ext.Append(".root" );
						output_scan_dat.Append(ext);
						ResultFormatter::FlatNTuplePullPlots( string( output_scan_dat ), AllContourFitResults[ii] );
					}
					
					//ResultFormatter::FlatNTuplePullPlots( string(output_scan_dat), AllContourFitResults.back() );
				}

			}
		}
		else
		{
			//Default action - presumably a fit or a toy study
			cerr << "No action performed" << endl;
			return 1;
		}
	}

	//Exit blurb
	time(&timeNow);
	cout << endl << "RapidFit" << endl;
	cout << "Ending time: " << ctime( &timeNow ) << endl;
	return 0;
}
