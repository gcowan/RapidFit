/**
  @file main.cpp

  Entry point for RapidFit

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

//  Root Headers
#include <TString.h>

//  System Headers
#include <string>
#include <vector>
#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

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
					argumentIndex++;
					Contour_X.push_back( argv[argumentIndex] );
					argumentIndex++;
					Contour_Y.push_back( argv[argumentIndex] );
				} else {
					cerr << "Contour Not Correctly Formatted" << endl;
					return 1;
				}
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
		XMLConfigReader * xmlFile=NULL;
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
		if( ( !RuntimeSeed.empty() && xmlFile->IsLoaded() ) && !UUID_Flag )
		{
			cout << "Setting Seed At Runtime to be: " << RuntimeSeed[0] << endl;
			xmlFile->SetSeed( RuntimeSeed[0] );
		} else if ( !RuntimeSeed.empty() && UUID_Flag && xmlFile->IsLoaded() )
		{
			cout << "Using a Random seed of UUID x RuntimeSeed to be 'unique' for all machines"<<endl;
			TUUID uid;
			UChar_t uuid[16];
			uid.GetUUID(uuid);
			RuntimeSeed[0] = RuntimeSeed[0] * ( uuid[ 2*(RuntimeSeed[0]%8) ]*256 +uuid[ 2*(RuntimeSeed[0]%8) ] );
		}

		//Create a parameter set
		ParameterSet * argumentParameterSet=NULL;
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
			if ( defineContourFlag )
			{
				for( unsigned short int i=0; i < Contour_X.size(); i++)
				{
					makeOutput->AddContour( Contour_X[i], Contour_Y[i] );
				}
			}
			if ( defineScanFlag )
			{
				for( unsigned short int i=0; i < Scan_X.size(); i++ )
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
			if ( ( ( numberRepeats > 1 ) && ( !doFC_Flag ) ) || doPullsFlag )
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
				//		This is re-used for FC scans and forms FC Step 1
				cout << "\n\n\t\tStarting Fit to Find Global Minima!\n"<<endl;
				//Do the fit to find GLOBAL MINIMA
				ToyStudyResult* GlobalFitResult = new ToyStudyResult( argumentParameterSet->GetAllNames() );
				GlobalFitResult->StartStopwatch();
				FitResult * GlobalResult = FitAssembler::DoSafeFit( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints() );

				GlobalFitResult->AddFitResult( GlobalResult );

				cout << "\n\n\t\tFit Output:" <<endl;
				//Output results
				makeOutput->SetInputResults( GlobalResult->GetResultParameterSet() );
				makeOutput->OutputFitResult( GlobalResult );


				//Do LL scan
				if( doLLscanFlag ) {
					//  Temporary Holder for the scan outputs to plot
					LLscanResult * llResult;
					vector<ToyStudyResult*> scanSoloResult;
					//  Store
					vector<LLscanResult*> scanResults;
					vector<string> LLscanList = makeOutput->GetScanList();

					for(unsigned int ii=0; ii < LLscanList.size() ; ii++)
					{
						ToyStudyResult* scan_result = FitAssembler::SingleScan( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints(), makeOutput, LLscanList[ii] );
						llResult = ResultFormatter::LLScan( scan_result, LLscanList[ii] );
						scanResults.push_back( llResult );
						scanSoloResult.push_back( scan_result );
					}
					makeOutput->SetLLscanFileName( LLscanFileName );
					makeOutput->OutputLLscanResult( scanResults ) ;
					for(unsigned int ii=0; ii < LLscanList.size(); ii++ )
					{
						TString output_scan_dat("LLScanData");
						output_scan_dat.Append(LLscanList[ii]);
						output_scan_dat.Append(".root" );
						ResultFormatter::FlatNTuplePullPlots( string( output_scan_dat ), scanSoloResult[ii] );
					}
				}

				//		This is re-used for FC scans and forms FC Step 2
				ToyStudyResult* _2DResultForFC=NULL;
				//Do 2D LL scan for deltaGamma and phis
				if( doLLcontourFlag || doFC_Flag ) {
					LLscanResult2D * llContourResult ;
					vector<LLscanResult2D*> contourResults ;

					vector<ToyStudyResult*> SoloContourResults;
//					vector<ToyStudyResult*> AllContourFitResults;
					vector<ToyStudyResult*> ContourLinearResults;

					vector<pair<string, string> > _2DLLscanList = makeOutput->Get2DScanList();

					if( _2DLLscanList.size() > 1 )
					{
						cerr << "\n\n\t\tI WILL NOT DO ONE 2D SCAN PER XML FILE, GO AWAY AND FIX YOURSELF\n" << endl;
					}
//					for(unsigned int ii=0; ii < _2DLLscanList.size() ; ii++ )
//					{
						string name1 = _2DLLscanList.back().first;
						string name2 = _2DLLscanList.back().second;
						SoloContourResults = FitAssembler::ContourScan( theMinimiser, theFunction, argumentParameterSet, pdfsAndData, xmlFile->GetConstraints(), makeOutput, name1, name2 );


						llContourResult = ResultFormatter::LLScan2D( SoloContourResults, name1, name2 );
						contourResults.push_back( llContourResult );

						TString ext_string("-");
						ext_string.Append( name1 );
						ext_string.Append("_");
						ext_string.Append( name2 );
						ext_string.Append( ".root" );
						LLcontourFileName.append( ext_string );

						ContourLinearResults.push_back( SoloContourResults[0] );

//						AllContourFitResults.push_back( ContourLinearResults );
//					}

					if( !doFC_Flag )
					{
					makeOutput->SetLLcontourFileName( LLcontourFileName );
					makeOutput->OutputLLcontourResult( contourResults ) ;
//					for(unsigned int ii=0; ii < _2DLLscanList.size(); ii++ )
//					{
//						TString output_scan_dat( "LLcontourScanData");
//						TString ext("_");
//						ext.Append(_2DLLscanList[ii].first);
//						ext.Append("_");
//						ext.Append(_2DLLscanList[ii].second);
//						ext.Append(".root" );
//						output_scan_dat.Append(ext);
						string output_scan_dat("LLcontourScanData.root");
						_2DResultForFC = new ToyStudyResult( ContourLinearResults );
						ResultFormatter::WriteFlatNtuple( output_scan_dat , ContourLinearResults[0] );
//					}
					}
				}


				if( doFC_Flag )
				{
					vector<ToyStudyResult*> AllResults;

					//		Want to loop over all points on a 2D Plot that have been allocated to this instance
					for( int iFC=0; iFC < _2DResultForFC->NumberResults(); iFC++ )
					{

						//		GET INPUT Data from fit Results
						//
						//  Get a ParameterSet that contains the input parameters from the output of a fit
					  	vector<pair<string, string> > _2DLLscanList = makeOutput->Get2DScanList();
						string name1 = _2DLLscanList[0].first;
						string name2 = _2DLLscanList[0].second;

						//		Use the inputs from Canonical Step 1
						ParameterSet* InputParamSet = GlobalResult->GetResultParameterSet()->GetDummyParameterSet();
						double lim1 = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetResultParameter( name1 )->GetValue();
						double lim2 = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetResultParameter( name2 )->GetValue();
						InputParamSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
						InputParamSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
						InputParamSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
						InputParamSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
						ParameterSet* ControlParamSet = GlobalResult->GetResultParameterSet()->GetDummyParameterSet();
						ControlParamSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
						ControlParamSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
						ParameterSet* InputFreeSet = GlobalResult->GetResultParameterSet()->GetDummyParameterSet();
						InputFreeSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
						InputFreeSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
						InputFreeSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
						InputFreeSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );

						//		Use the inputs from Canonical Step 2
						//ParameterSet* InputParamSet = _2DResultForFC[iFC]->GetResultParameterSet()->GetDummyParameterSet();
						//ParameterSet* InputFreeSet = _2DResultForFC[iFC]->GetResultParameterSet()->GetDummyParameterSet();

						//		Just for clarity
						//		also when using values from Step1 these are not set correctly for the study
						InputParamSet->GetPhysicsParameter( name1 )->SetType( "Fixed" );
						InputParamSet->GetPhysicsParameter( name2 )->SetType( "Fixed" );
						InputFreeSet->GetPhysicsParameter( name1 )->SetType( "Free" );
						InputFreeSet->GetPhysicsParameter( name2 )->SetType( "Free" );



						//	Collect all of the relevent Data from the XML
						//	Note: most of these had to be written for FCscans
						MinimiserConfiguration * ToyStudyMinimiser = theMinimiser;
						FitFunctionConfiguration* ToyStudyFunction = theFunction;
						vector<vector<IPrecalculator*> > ToyPrecalculators = xmlFile->GetPrecalculators();
						vector<PhaseSpaceBoundary*> PhaseSpaceForToys = xmlFile->GetPhaseSpaceBoundaries();
						vector<ConstraintFunction*> ConstraintsForToys = xmlFile->GetConstraints();

						//	Read the number of events from the existing DataSets
						//	I think there needs be an equivalent function drafted to use sWeights
						vector<unsigned int> EventsPerPDF;
						for( unsigned int pdf_num=0; pdf_num < pdfsAndData.size(); pdf_num++ )
						{
							EventsPerPDF.push_back( pdfsAndData[pdf_num]->GetDataSet()->GetDataNumber() );
						}



						//		Number of datasets wanted i.e. Number of Toys
						int wanted_number_of_toys = 100;
						if( numberRepeatsFlag ) wanted_number_of_toys = numberRepeats;



						//		GENERATE DATA
						//	We want to now generate a new DataSet for a Toy Study
						//	I choose a Foam DataSet with no Cuts or extra arguments as these are not required
						//
						//	Generate Once, Fit twice
						//		This was a bit of a pain to code up.
						//		Although I like the idea of coding it up into a GenerateToyData object
						vector<PDFWithData*> PDFsWithDataForToys;
						vector<vector<IDataSet*> > Memory_Data;
						Memory_Data.resize( pdfsAndData.size() );

						//	We may have multiple PDFs in the XML
						for( unsigned short int pdf_num=0; pdf_num < pdfsAndData.size(); pdf_num++ )
						{
							//	Generate the Data for this Study using the given the XML PhaseSpace and PDF
							IPDF* PDF_from_XML = pdfsAndData[pdf_num]->GetPDF();

							//	Create a DataSetConfiguration to make the datasets
							vector<string> empty_args;
							vector<string> empty_arg_names;
							DataSetConfiguration* Toy_Foam_DataSet= new DataSetConfiguration( "Foam", EventsPerPDF[pdf_num], "", empty_args, empty_arg_names, PDF_from_XML );
							//	Set the Input Physics Parameters to be as defined above ControlSet
							//	This is to seperate changing bottles for Physics from that used at Generation
							//	as we want the data to have a wider scope than 1 fit with the data
							Toy_Foam_DataSet->SetPhysicsParameters( ControlParamSet );


							cout << "\n\n\t\tCACHING DATA FOR TOYS, THIS WILL LIKELY TAKE A LONG TIME!\n\n" <<endl;
							cout << "Caching Data for: " << wanted_number_of_toys << " toys."<<endl;
							for( short int dataset_num=0; dataset_num < wanted_number_of_toys; dataset_num++ )
							{
								//	Make the data
								cout << "Making Data...\t";
								//	...I can't seem to workout why these outputs are not linear
								//	even without compiler optimisations this code is a single for loop?...
								IDataSet* new_dataset = Toy_Foam_DataSet->MakeDataSet( PhaseSpaceForToys[pdf_num], PDF_from_XML );
								//	Store the Data Object in Memory
								Memory_Data[pdf_num].push_back( new_dataset );
								cout << "Cached PDF: " << pdf_num << " DataSet: " << dataset_num <<"\t";
							}
							cout << endl;

							//	Store the information for each PDFWithData
							vector<DataSetConfiguration*> DataSetConfigForToys;
							DataSetConfigForToys.push_back( Toy_Foam_DataSet );
							PDFWithData* ToyPDFWithData = new PDFWithData( PDF_from_XML, PhaseSpaceForToys[pdf_num], DataSetConfigForToys, ToyPrecalculators[pdf_num] );

							//	Store this PDFWithData
							ToyPDFWithData->SetPhysicsParameters( InputParamSet );
							PDFsWithDataForToys.push_back( ToyPDFWithData );
						}

						//		Result Vectors for the Fit for clarity and Output Formatting
						ToyStudyResult* study1Results = new ToyStudyResult( GlobalResult->GetResultParameterSet()->GetAllNames() );
						ToyStudyResult* study2Results = new ToyStudyResult( GlobalResult->GetResultParameterSet()->GetAllNames() );


						//	Try to make minuit as quiet as possible here... does someone have a way to gag it's output?!?
						ToyStudyMinimiser->GetMinimiser( int(InputParamSet->GetAllNames().size()) )->SetOutputLevel(-1);


						cout << "\n\n\t\tPerforming Fits to Toys in FC\n"<<endl;
						//		Perform Fits
						//	This will record ONLY Data from working fits i.e. FitStatus==3
						for( unsigned short int dataset_num=0; dataset_num < Memory_Data[0].size(); dataset_num++ )
						{
							
							FitResult* fit1Result = NULL;
							FitResult* fit2Result = NULL;
							bool toy_failed = false;


							//	Assuming Nuisence Parameters set to Global Minima
							//	We need to set this for EVERY FIT in order to have correct generation/pull values
							ParameterSet* LocalInputFreeSet = GlobalResult->GetResultParameterSet()->GetDummyParameterSet();
							LocalInputFreeSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
							LocalInputFreeSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
							LocalInputFreeSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
							LocalInputFreeSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
							LocalInputFreeSet->GetPhysicsParameter( name1 )->SetType( "Free" );
							LocalInputFreeSet->GetPhysicsParameter( name2 )->SetType( "Free" );
							ParameterSet* LocalInputFixedSet = GlobalResult->GetResultParameterSet()->GetDummyParameterSet();
							LocalInputFixedSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
							LocalInputFixedSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
							LocalInputFixedSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
							LocalInputFixedSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
							LocalInputFixedSet->GetPhysicsParameter( name1 )->SetType( "Fixed" );
							LocalInputFixedSet->GetPhysicsParameter( name2 )->SetType( "Fixed" );

							//	Assuming Nuisence Parameters set to Local Minima
							//
							//ParameterSet* LocalInputFreeSet = _2DResultForFC[iFC]->GetResultParameterSet()->GetDummyParameterSet();
							//LocalInputFreeSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
							//LocalInputFreeSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
							//LocalInputFreeSet->GetPhysicsParameter( name1 )->SetType( "Free" );
							//LocalInputFreeSet->GetPhysicsParameter( name2 )->SetType( "Free" );
							//ParameterSet* LocalInputFixedSet = _2DResultForFC[iFC]->GetResultParameterSet()->GetDummyParameterSet();
							//LocalInputFixedSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
							//LocalInputFixedSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
							//LocalInputFixedSet->GetPhysicsParameter( name1 )->SetType( "Fixed" );
							//LocalInputFixedSet->GetPhysicsParameter( name2 )->SetType( "Fixed" );


							//	We need to set some factors before we perform the fit
							for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); pdf_num++ )
							{
								vector<IDataSet*> wanted_set;
								wanted_set.push_back( Memory_Data[pdf_num][dataset_num] );
								PDFsWithDataForToys[pdf_num]->AddCachedData( wanted_set );
								PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( LocalInputFreeSet );
							}

							cout << "\n\n\t\tPerforming Fit To Toy "<< dataset_num <<" of " << Memory_Data[0].size() << endl;
							//	Fit once with control parameters Free
							fit1Result = FitAssembler::DoSafeFit( ToyStudyMinimiser, ToyStudyFunction, LocalInputFreeSet, PDFsWithDataForToys, ConstraintsForToys );

							//	Only Fit again to this dataset if it fits well with +2 dof
							//	This has the obvious savings in CPU resources
							if( ( fit1Result->GetFitStatus() == 3 ) || FC_Debug_Flag )
							{
								//  Fit second with control parameters Fixed
								for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); pdf_num++ )
								{
									vector<IDataSet*> wanted_set;
									wanted_set.push_back( Memory_Data[pdf_num][dataset_num] );
									PDFsWithDataForToys[pdf_num]->AddCachedData( wanted_set );
									PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( LocalInputFixedSet );
								}

								cout << "\n\n\t\tFirst Fit Successful, Performing a Second Fit " << dataset_num << " of " << Memory_Data[0].size() <<endl;
								fit2Result = FitAssembler::DoSafeFit( ToyStudyMinimiser, ToyStudyFunction, LocalInputFixedSet, PDFsWithDataForToys, ConstraintsForToys );

								//	If either Fit Failed we want to 'dump the results' and run an extra Fit.
								if( (fit2Result->GetFitStatus() != 3) || (fit2Result->GetFitStatus() != 3) ) toy_failed = true;
							} else {
								toy_failed = true;
							}

							//	Do we want to store the Data OR run another toy to get a better Fit
							if( toy_failed )
							{
								cerr << "\n\n\t\tA Single Toy Study Failed... Requesting More Data for another pass.\n" << endl;
								for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); pdf_num++ )
								{
									cerr << "Adding More Data for pdf: " << pdf_num << "\t";
									PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( ControlParamSet );
									IDataSet* new_dataset = PDFsWithDataForToys[pdf_num]->GetDataSetConfig()->MakeDataSet( PhaseSpaceForToys[pdf_num], PDFsWithDataForToys[pdf_num]->GetPDF() );
									Memory_Data[pdf_num].push_back( new_dataset );
								}
							}
							if( FC_Debug_Flag || !toy_failed ){
								if( toy_failed )
								{
									fit1Result->ForceFitStatus(-2);
									fit2Result->ForceFitStatus(-2);
								}
								study1Results->AddFitResult( fit1Result );
								study2Results->AddFitResult( fit2Result );
							}
						}


						//			STEP 6
						//

						//		Standard Output Format which makes the results the same running either
						//		the whole scan on one machine or running the whole set on a batch system
						ToyStudyResult* ThisStudy = new ToyStudyResult( GlobalResult->GetResultParameterSet()->GetAllNames() );
						ThisStudy->AddFitResult( _2DResultForFC->GetFitResult( iFC ), false );
						ThisStudy->AddCPUTimes( _2DResultForFC->GetAllCPUTimes() );
						ThisStudy->AddRealTimes( _2DResultForFC->GetAllRealTimes() );
						//	The Generated Value for the Global and Local fit are best defined as -9999 as a sensible default
						for( unsigned short int num=0; num < GlobalResult->GetResultParameterSet()->GetAllNames().size(); num++ )
						{
							string name = GlobalResult->GetResultParameterSet()->GetAllNames()[num];
							GlobalResult->GetResultParameterSet()->GetResultParameter( name )->ForcePullValue( -9999 );
							GlobalResult->GetResultParameterSet()->GetResultParameter( name )->ForceOriginalValue( -9999 );
							ThisStudy->GetFitResult(0)->GetResultParameterSet()->GetResultParameter( name )->ForcePullValue( -9999 );
							ThisStudy->GetFitResult(0)->GetResultParameterSet()->GetResultParameter( name )->ForceOriginalValue( -9999 );
						}
						AllResults.push_back( GlobalFitResult );
						AllResults.push_back( ThisStudy );
						AllResults.push_back( study1Results );
						AllResults.push_back( study2Results );

					}

					//		STORE THE OUTPUT OF THE TOY STUDIES
					ToyStudyResult* AllFlatResult = new ToyStudyResult( AllResults );
					ResultFormatter::WriteFlatNtuple( "FCOutput.root", AllFlatResult );

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
