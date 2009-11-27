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
#include "ClassLookUp.h"
#include "RapidFitIntegrator.h"
#include "Plotter.h"
#include <stdio.h>
#include <stdlib.h>
#include "MakeFoam.h"

using namespace std;

void OutputOneFit( FitResult*, bool, string );

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
		ResultFormatter::MakePullPlots( "pullPlots.root", fitResults );
		ResultFormatter::LatexOutputFitResult( fitResults->GetFitResult(0) );
	}
	else
	{
		//Variables to store command line arguments
		int numberRepeats = 0;
		string configFileName = "";
		vector<string> parameterTemplates;
		string minimiserName = "";
		FitFunction * theFunction;
		string saveOneDataSetFileName = "";
		string plotFileName = "FitPlots.root";
		string pullFileName = "PullPlots.root";
		vector< PDFWithData* > pdfsAndData;

		//Flags for which arguments have been received
		bool numberRepeatsFlag = false;
		bool configFileNameFlag = false;
		bool parameterTemplateFlag = false;
		bool minimiserNameFlag = false;
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
					minimiserName = argv[argumentIndex];
					minimiserNameFlag = true;
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
					theFunction = ClassLookUp::LookUpFitFunctionName( argv[argumentIndex], "" );
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
				XMLConfigReader quickConfig(configFileName);
				PDFWithData * quickData = quickConfig.GetPDFsAndData()[0];
				quickData->SetPhysicsParameters( quickConfig.GetFitParameters() );
				ResultFormatter::MakeRootDataFile( saveOneDataSetFileName, quickData->GetDataSet() );
			}
			else
			{
				cerr << "No data set specified" << endl;
			}
		}
		else if (testIntegratorFlag)
		{
			if (configFileNameFlag)
			{
				//Compare numerical and analytical integration
				XMLConfigReader quickConfig(configFileName);
				PDFWithData * quickData = quickConfig.GetPDFsAndData()[0];
				quickData->SetPhysicsParameters( quickConfig.GetFitParameters() );
				IDataSet * quickDataSet = quickData->GetDataSet();
				RapidFitIntegrator * testIntegrator = new RapidFitIntegrator( quickData->GetPDF() );
				testIntegrator->Integral( quickDataSet->GetDataPoint(0), quickDataSet->GetBoundary() );
			}
			else
			{
				cerr << "No data set specified" << endl;
			}
		}
		else if (testRapidIntegratorFlag)
		{
			if (configFileNameFlag)
			{
				//Compare numerical and analytical integration
				XMLConfigReader quickConfig(configFileName);
				PDFWithData * quickData = quickConfig.GetPDFsAndData()[0];
				quickData->SetPhysicsParameters( quickConfig.GetFitParameters() );
				IDataSet * quickDataSet = quickData->GetDataSet();
				MakeFoam * testFoam = new MakeFoam( quickData->GetPDF(), quickDataSet->GetBoundary(), quickDataSet->GetDataPoint(0) );
				testFoam->Debug();
			}
			else
			{
				cerr << "No data set specified" << endl;
			}
		}
		else if (testPlotFlag)
		{
			if (configFileNameFlag)
			{
				//Project the PDF onto the data
				XMLConfigReader quickConfig(configFileName);
				PDFWithData * quickData = quickConfig.GetPDFsAndData()[0];
				quickData->SetPhysicsParameters( quickConfig.GetFitParameters() );
				IDataSet * quickDataSet = quickData->GetDataSet();
				Plotter * testPlotter = new Plotter( quickData->GetPDF(), quickDataSet );
				testPlotter->PlotAllObservables(plotFileName);
			}
			else
			{
				cerr << "No data set specified" << endl;
			}
		}
		else if (configFileNameFlag)
		{
			//Actually configure a fit

			//Command line argments override the config file
			if (!minimiserNameFlag)
			{
				minimiserName = xmlFile->GetMinimiserName();
			}
			if (!theFunctionFlag)
			{
				theFunction = xmlFile->GetFitFunction();
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

			//Pick a toy study if there are repeats, or if pull plots are wanted
			if ( numberRepeats > 1 || doPullsFlag )
			{
				ToyStudy newStudy( minimiserName, theFunction, argumentParameterSet, pdfsAndData, numberRepeats );
				ToyStudyResult * fitResults = newStudy.DoWholeStudy();
				ResultFormatter::MakePullPlots( pullFileName, fitResults );

				OutputOneFit( fitResults->GetFitResult(0), doPlottingFlag, plotFileName );
			}
			else
			{
				FitResult * oneResult = FitAssembler::DoFit( minimiserName, theFunction, argumentParameterSet, pdfsAndData );
				OutputOneFit( oneResult, doPlottingFlag, plotFileName );
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

void OutputOneFit( FitResult * OneResult, bool DoPlotting, string PlotFileName )
{
	ResultFormatter::LatexOutputFitResult( OneResult );
	ResultFormatter::LatexOutputCovarianceMatrix( OneResult );
	ResultFormatter::PlotFitContours( OneResult, "contours" );

	if (DoPlotting)
	{
		PhysicsBottle * resultBottle = OneResult->GetPhysicsBottle();

		//Loop over all PDFs, and plot
		for ( int resultIndex = 0; resultIndex < resultBottle->NumberResults(); resultIndex++ )
		{
			time_t timeNow1;
			time(&timeNow1);
			cout << "Plotting parameter: " << ctime( &timeNow1 ) << endl;

			Plotter * testPlotter = new Plotter( resultBottle->GetResultPDF(resultIndex), resultBottle->GetResultDataSet(resultIndex) );
			char fileNumber[100];
			sprintf( fileNumber, "%d.", resultIndex );
			testPlotter->PlotAllObservables( fileNumber + PlotFileName );
		}
	}
}
