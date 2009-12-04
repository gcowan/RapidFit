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
		vector< PDFWithData* > pdfsAndData;

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
		bool testRapidIntegratorFlag = false;

		//Parse command line arguments
		string currentArgument;
		for ( int argumentIndex = 1; argumentIndex < argc; argumentIndex++ )
		{
			currentArgument = argv[argumentIndex];

			if ( currentArgument == "-f" )
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
			else
			{
				cerr << "Unrecognised argument: " << currentArgument << endl;
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
		else if (testRapidIntegratorFlag)
		{
			if (configFileNameFlag)
			{
				//Compare numerical and analytical integration
				PDFWithData * quickData = xmlFile->GetPDFsAndData()[0];
				quickData->SetPhysicsParameters( xmlFile->GetFitParameters() );
				IDataSet * quickDataSet = quickData->GetDataSet();
				MakeFoam * testFoam = new MakeFoam( quickData->GetPDF(), quickDataSet->GetBoundary(), quickDataSet->GetDataPoint(0) );
				testFoam->Debug();
			}
			else
			{
				cerr << "No data set specified" << endl;
				return 1;
			}
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
				ToyStudy newStudy( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, numberRepeats );
				ToyStudyResult * fitResults = newStudy.DoWholeStudy();

				//Output results
				makeOutput->OutputToyResult(fitResults);
				makeOutput->OutputFitResult( fitResults->GetFitResult(0) );
			}
			else
			{
				//Do the fit
				FitResult * oneResult = FitAssembler::DoFit( theMinimiser, theFunction, argumentParameterSet, pdfsAndData );

				//Output results
				makeOutput->OutputFitResult(oneResult);
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
