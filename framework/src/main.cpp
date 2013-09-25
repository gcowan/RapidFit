/**
 * @file main.cpp
 *
 * Entry point for RapidFit
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

///  Root Headers
#include "TString.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TSystem.h"
#include "TROOT.h"
///  RapidFit Headers
#include "Mathematics.h"
#include "FitAssembler.h"
#include "ToyStudy.h"
#include "XMLConfigReader.h"
#include "ResultFormatter.h"
#include "InputParsing.h"
#include "RapidFitIntegrator.h"
#include "MakeFoam.h"
#include "PerEventAngularAcceptance.h"
#include "OutputConfiguration.h"
#include "FitResultVector.h"
#include "main.h"
#include "DataSetConfiguration.h"
#include "IPDF.h"
#include "IDataSet.h"
#include "MemoryDataSet.h"
#include "StringProcessing.h"
#include "MCStudy.h"
#include "GoodnessOfFit.h"
#include "ComponentPlotter.h"
#include "ParameterSet.h"
#include "IConstraint.h"
#include "ClassLookUp.h"
#include "ScanStudies.h"
#include "VectoredFeldmanCousins.h"
#include "PDFConfigurator.h"
#include "RapidFitConfiguration.h"
#include "ParseCommandLine.h"
#include "JackKnife.h"
#include "RapidRun.h"
#include "ResultFormatter.h"
///  System Headers
#include <string>
#include <vector>
#include <iostream>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

using namespace::std;

//	Standard opening into RapidFit when not in CINT/pyROOT
#ifndef __CINT__
int main( int argc, char * argv[] )
{
	vector<string> input;
	for( int i=0; i< argc; ++i )
	{
		input.push_back( argv[i] );
	}
	return RapidFit( input );
}
#endif

void RapidFitWelcome()
{
	//Welcome blurb
	time_t timeNow;
	time(&timeNow);

	cout << endl << "RapidFit" << endl;
	cout << "Framework SVN Rev:\t" << STR(SVN_REV) << endl;
	cout << "PDF SVN REV:\t\t" << STR(SVN_PDF_REV) << endl;

#ifdef __RAPIDFIT_USE_GSL
	cout << "Compiled With Optional GSL Components" << endl;
#endif

	string mod="M";
	size_t found = string( STR(SVN_REV) ).find( mod );
	if( found != string::npos )
	{
		cout << endl << "\t\t!!!YOU HAVE LOCAL CODE MODIFICATIONS!!!" << endl << endl;
	}
	cout << "Build Date:\t\t" << STR(BUILD_DATE) << endl;
	cout << "Built Against ROOT:\t" << STR(ROOT_RELEASE) << endl;
	cout << "Starting time:\t\t" << ctime( &timeNow ) << endl << endl;
}

void RapidFitExit()
{
	time_t timeNow;
	//Exit blurb
	time(&timeNow);
	cout << endl << "RapidFit" << endl;
	cout << "Ending time: " << ctime( &timeNow ) << endl;

	cout << "Goodbye :)" << endl;
}

int BuildTemplateXML( RapidFitConfiguration* config );

int ConfigureRapidFit( RapidFitConfiguration* config );

int saveOneDataSet( RapidFitConfiguration* config );

int testIntegrator( RapidFitConfiguration* config );

int testComponentPlot( RapidFitConfiguration* config );

int calculateFitFractions( RapidFitConfiguration* config );

int calculateAcceptanceWeights( RapidFitConfiguration* config );

int calculateAcceptanceCoefficients( RapidFitConfiguration* config );

int calculateAcceptanceWeightsWithSwave( RapidFitConfiguration* config );

int calculatePerEventAcceptance( RapidFitConfiguration* config );

int PerformToyStudy( RapidFitConfiguration* config );

int PerformMCStudy( RapidFitConfiguration* config );

int PerformFCStudy( RapidFitConfiguration* config );

int PerformMainFit( RapidFitConfiguration* config );

int PerformLLScan( RapidFitConfiguration* config );

int Perform2DLLScan( RapidFitConfiguration* config );

int PerformJackKnife( RapidFitConfiguration* config );

void SaveXML( RapidFitConfiguration* config );

//	The 'meat' of the Code
int RapidFit( vector<string> input )
{
	//  reduce root verbosity:
	//  gErrorIgnoreLevel = kInfo; // The default is kError
	//  gErrorIgnoreLevel = kFatal; // Explicitly remove all messages
	//  gErrorIgnoreLevel = kError; // Only error message (default)
	gErrorIgnoreLevel = kWarning; // error and warning message
	//  gErrorIgnoreLevel = kNote; // error, warning and note
	//  gErrorIgnoreLevel = kInfo; // Display all information (same as -v)

	RapidFitWelcome();

	RapidFitConfiguration* thisConfig = new RapidFitConfiguration();

	thisConfig->runtimeArgs = input;

	int command_check = ParseCommandLine::ParseThisCommandLine( *thisConfig, input );

	if( thisConfig->debug != NULL )
	{
		if( thisConfig->debug->DebugThisClass( "main" ) )
		{
			cout << endl;
			cout << "Finished Procesing Command Line" << endl;
			cout << "main: Debugging Enabled" << endl;
		}
	}


	if( command_check != 0 )
	{
		return command_check;
	}

	if( thisConfig->makeTemplateXML )
	{
		BuildTemplateXML( thisConfig );
		exit(0);
	}

	if( thisConfig->debug != NULL )
	{
		if( thisConfig->debug->DebugThisClass( "main" ) )
		{
			cout << endl;
			cout << "About to Configure RapidFit with the given XML and command line options" << endl;
		}
	}

	ConfigureRapidFit( thisConfig );

	if( thisConfig->debug != NULL )
	{
		if( thisConfig->debug->DebugThisClass( "main" ) )
		{
			cout << endl;
			cout << "RapidFit has been Configured and initial objects in XML/CommandLine have been constructed" << endl;
		}
	}

	int main_fitResult=0;

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
	//	9b)	Weight the dataset if requested to do so and save the output(s)

	//	10)	If requested perform 1D LL Scans
	//	12)	If requested perform 2D LL Scans
	//	13)	Perform a FC scan if requested

	//	13)	Universal Exit


	//	1)
	//	If we want to save one dataset and have an input XML
	if( thisConfig->saveOneDataSetFlag ) saveOneDataSet( thisConfig );

	//	2)
	else if( thisConfig->testIntegratorFlag && thisConfig->configFileNameFlag) testIntegrator( thisConfig );

	//	3)
	else if( thisConfig->calculateAcceptanceWeights && thisConfig->configFileNameFlag ) calculateAcceptanceWeights( thisConfig );
	else if( thisConfig->calculateAcceptanceCoefficients && thisConfig->configFileNameFlag ) calculateAcceptanceCoefficients( thisConfig );

	//	4)
	else if( thisConfig->calculateAcceptanceWeightsWithSwave && thisConfig->configFileNameFlag ) calculateAcceptanceWeightsWithSwave( thisConfig );

	//	5)
	else if( thisConfig->calculatePerEventAcceptance) calculatePerEventAcceptance( thisConfig );

	//	6)
	else if( thisConfig->testComponentPlotFlag && thisConfig->configFileNameFlag && thisConfig->observableNameFlag ) testComponentPlot( thisConfig );

	else if( thisConfig->JackKnife_Flag ) PerformJackKnife( thisConfig );

	//	7)	Toy Study
	//Pick a toy study if there are repeats, or if pull plots are wanted
	else if( (thisConfig->numberRepeats > 1 || thisConfig->doPullsFlag ) && !( thisConfig->doFC_Flag || thisConfig->doLLscanFlag || thisConfig->doLLcontourFlag ) )
	{
		PerformToyStudy( thisConfig );
	}

	//	8)	MC Study
	//	9)
	else if( !thisConfig->FC_LL_PART_Flag ) main_fitResult = PerformMainFit( thisConfig );

	//	10)
	//	Do LL scan
	if( main_fitResult>-1 && ( thisConfig->doLLscanFlag || ( thisConfig->doFC_Flag && thisConfig->makeOutput->Get2DScanList().empty() ) ) ) PerformLLScan( thisConfig );

	//		This is re-used for FC scans and forms FC Step 2
	//	11)
	//Do 2D LL scan
	if( main_fitResult>-1 )
		if( thisConfig->doLLcontourFlag || ( thisConfig->doFC_Flag && !thisConfig->makeOutput->Get2DScanList().empty() ) )
			if( !thisConfig->FC_LL_PART_Flag && !thisConfig->doLLscanFlag ) Perform2DLLScan( thisConfig );

	//	12b)
	//	Do the main work of the FC scan
	if( main_fitResult>-1 && (thisConfig->doFC_Flag && !thisConfig->doLLscanFlag )) PerformFCStudy( thisConfig );

	if( main_fitResult > -1 && thisConfig->calculateFitFractionsFlag ) calculateFitFractions( thisConfig );

	//	Should only happen under the condition that no CV fit was performed or anything else
	//if( thisConfig->GlobalFitResult == NULL )
	//{
	//	//	Default action - presumably a fit or a toy study
	//	cerr << "No action performed" << endl;
	//	cerr << "Not sure how I got here, please email a maintainer!..." <<endl;
	//}


	cout << "Any Fit Results and Projection Outputs are Stored in: " << ResultFormatter::GetOutputFolder().c_str() << endl;

	//	thisConfig performs a cleanup on the ResultFormatter outputFolder object.
	//	This is safe as long as it's not needed after here!
	delete thisConfig;

	RapidFitExit();

	return 0;
}

int PerformJackKnife( RapidFitConfiguration* config )
{
	JackKnife::jackknife( config->xmlFile, config->theMinimiser, config->theFunction, config->argumentParameterSet, config->CommandLineParamvector, config->jackStartNum, config->jackStopNum );
	return 0;
}

int Perform2DLLScan( RapidFitConfiguration* config )
{
	vector<pair<string, string> > _2DLLscanList = config->makeOutput->Get2DScanList();

	unsigned int initial_scan=0;
	if( ( ( _2DLLscanList.size() > 1 ) && config->doFC_Flag ) || config->defineContourFlag )
	{
		initial_scan = unsigned( _2DLLscanList.size()-1 );
	}

	if( ( ( _2DLLscanList.size() == 0 ) && ( config->doFC_Flag) ) && !config->makeOutput->GetScanList().empty() )
	{
		cerr << "\n\n\t\tNO 2D SCAN DATA, I'M NOT GOING TO DO FC, GO AWAY!\n\n"<<endl;
		exit(-54);
	}

	for(unsigned int ii=initial_scan; ii < _2DLLscanList.size() ; ++ii )
	{
		string name1 = _2DLLscanList[ii].first;
		string name2 = _2DLLscanList[ii].second;
		if( config->StartAtCenterFlag )
		{
			ParameterSet* param_set = config->GlobalResult->GetResultParameterSet()->GetDummyParameterSet();
			for(unsigned int i=0; i< config->argumentParameterSet->GetAllNames().size(); ++i )
			{
				config->argumentParameterSet->GetPhysicsParameter( config->argumentParameterSet->GetAllNames()[i] )
					->SetBlindedValue( param_set->GetPhysicsParameter( config->argumentParameterSet->GetAllNames()[i] )->GetValue() );
			}
		}

		vector<FitResultVector*> Temp_Results = ScanStudies::ContourScan( config->theMinimiser, config->theFunction, config->argumentParameterSet,

				config->pdfsAndData, config->xmlFile->GetConstraints(), config->makeOutput, name1, name2, config->OutputLevel2, config->debug, config->Force_Continue_Flag );

		vector<FitResultVector*> Ordered_Results;

		for( vector<FitResultVector*>::iterator result_i = Temp_Results.begin(); result_i != Temp_Results.end(); ++result_i )
		{
			Ordered_Results.push_back( config->GlobalFitResult );
			Ordered_Results.push_back( *result_i );
		}

		config->SoloContourResults.push_back( new FitResultVector( Ordered_Results ) );
	}

	//	12a)
	//  Don't output the scan files when doing a FC scan
	if( config->doFC_Flag )
	{
		config->_2DResultForFC = config->SoloContourResults.back();
	}
	else
	{
		for(unsigned int scanNum=0; scanNum < _2DLLscanList.size(); ++scanNum )
		{
			TString output_scan_dat("LLScanData");
			TString time_stamped_name( output_scan_dat ); time_stamped_name.Append( "_" ); time_stamped_name.Append( StringProcessing::TimeString() );
			output_scan_dat.Append(".root"); time_stamped_name.Append(".root");
			//	Add the Global Results and 'Linearize' the output
			vector<FitResultVector*> TempContourResults;
			config->GlobalResult->GetResultParameterSet()->GetResultParameter( _2DLLscanList[scanNum].first )->SetScanStatus( true );
			config->GlobalResult->GetResultParameterSet()->GetResultParameter( _2DLLscanList[scanNum].second )->SetScanStatus( true );
			ResultFormatter::WriteFlatNtuple( string( output_scan_dat ), config->SoloContourResults[scanNum], config->xmlFile->GetXML(), config->runtimeArgs );
			ResultFormatter::WriteFlatNtuple( string( time_stamped_name ), config->SoloContourResults[scanNum], config->xmlFile->GetXML(), config->runtimeArgs );
			config->GlobalResult->GetResultParameterSet()->GetResultParameter( _2DLLscanList[scanNum].first )->SetScanStatus( false );
			config->GlobalResult->GetResultParameterSet()->GetResultParameter( _2DLLscanList[scanNum].second )->SetScanStatus( false );
		}
	}
	return 0;
}

int PerformLLScan( RapidFitConfiguration* config )
{
	vector<FitResultVector*> scanSoloResults;

	//  Store
	vector<string> LLscanList = config->makeOutput->GetScanList();

	for( unsigned int scan_num=0; scan_num < LLscanList.size() ; ++scan_num)
	{
		cout << "Scanning: " << LLscanList[scan_num] << endl;
	}

	for(unsigned int scan_num=0; scan_num < LLscanList.size() ; ++scan_num)
	{
		if( config->StartAtCenterFlag )
		{
			ParameterSet*  param_set = config->GlobalResult->GetResultParameterSet()->GetDummyParameterSet();
			for(unsigned int i=0; i< config->argumentParameterSet->GetAllNames().size(); ++i )
			{
				config->argumentParameterSet->GetPhysicsParameter( config->argumentParameterSet->GetAllNames()[i] )->SetBlindedValue(
						param_set->GetPhysicsParameter( config->argumentParameterSet->GetAllNames()[i] )->GetValue() );
			}
		}

		FitResultVector* scan_result = ScanStudies::SingleScan( config->theMinimiser, config->theFunction, config->argumentParameterSet, config->pdfsAndData,

				config->xmlFile->GetConstraints(), config->makeOutput, LLscanList[scan_num], config->OutputLevel2, config->debug, config->Force_Continue_Flag );

		scanSoloResults.push_back( scan_result );

		cout << "Scan Finished" << endl;

		if( config->doFC_Flag ) cout << "Wanting FC Scan" << endl;
		else cout << "NOT Wanting FC Scan" << endl;

		if( config->doFC_Flag )
		{
			cout << "Performing 1DFC Scan" << endl;
			FitResultVector* new_1D = new FitResultVector( scanSoloResults );

			//config->GlobalResult->GetResultParameterSet()->GetResultParameter( LLscanList[scan_num] )->SetScanStatus( true );

			if( config->debug != NULL )
			{
				if( config->debug->DebugThisClass( "main" ) )
				{
					cout << "main: Constructing FeldmanCousins" << endl;
				}
			}

			for( unsigned int i=0; i< config->pdfsAndData.size(); ++i )	config->pdfsAndData[i]->ClearCache();

			VectoredFeldmanCousins* new_study = new VectoredFeldmanCousins( config->GlobalFitResult, new_1D, config->Nuisencemodel, config->makeOutput, config->theMinimiser,
					config->theFunction, config->xmlFile, config->pdfsAndData );
			new_study->SetNumRepeats( config->numberRepeats );

			if( config->debug != NULL )
			{
				if( config->debug->DebugThisClass( "main" ) )
				{
					cout << "main: Starting FeldmanCousins" << endl;
				}
			}

			new_study->DoWholeStudy( config->OutputLevel2 );
			//	Doesn't hurt to be sure we obay the file format standard
			vector<FitResultVector*> file_output;
			file_output.push_back( config->GlobalFitResult );
			FitResultVector* study_output = new_study->GetStudyResult();
			file_output.push_back( study_output );
			FitResultVector* for_file = new FitResultVector( file_output );
			//	Making the assumption the user isn't running more than one of these at a time and isn't an idiot
			ResultFormatter::WriteFlatNtuple( "1DLL_FCScan.root", for_file, config->xmlFile->GetXML(), config->runtimeArgs );
			ResultFormatter::WriteFlatNtuple( "FCScan.root", for_file, config->xmlFile->GetXML(), config->runtimeArgs );
			config->GlobalResult->GetResultParameterSet()->GetResultParameter( LLscanList[scan_num] )->SetScanStatus( false );
		}

		for(unsigned int this_scan_num=0; this_scan_num < LLscanList.size(); ++this_scan_num )
		{
			TString new_output_scan_dat( "LLScanData" );
			TString time_stamped_name( new_output_scan_dat ); time_stamped_name.Append( "_" ); time_stamped_name.Append( StringProcessing::TimeString() );
			new_output_scan_dat.Append(".root"); time_stamped_name.Append(".root");
			vector<FitResultVector*> ammended_format;
			config->GlobalResult->GetResultParameterSet()->GetResultParameter( string(LLscanList[this_scan_num]) )->SetScanStatus( true );

			//	The output file format is [0] = Global_CV, [1] = Scan_CV_1, [2] = Global_CV, [3] = Scan_CV_2 ...
			for( int i=0; i< scanSoloResults[this_scan_num]->NumberResults(); ++i )
			{
				ammended_format.push_back( config->GlobalFitResult );
				FitResultVector* temp_vec = new FitResultVector( scanSoloResults[this_scan_num]->GetAllNames() );
				temp_vec->AddFitResult( scanSoloResults[this_scan_num]->GetFitResult( i ), false );
				temp_vec->AddRealTime( scanSoloResults[this_scan_num]->GetRealTime(i) );
				temp_vec->AddCPUTime( scanSoloResults[this_scan_num]->GetCPUTime(i) );
				temp_vec->AddGLTime( scanSoloResults[this_scan_num]->GetGLTime(i) );
				ammended_format.push_back( temp_vec );
			}
			FitResultVector* corrected_format = new FitResultVector( ammended_format );
			ResultFormatter::WriteFlatNtuple( string( new_output_scan_dat ), corrected_format, config->xmlFile->GetXML(), config->runtimeArgs );
			ResultFormatter::WriteFlatNtuple( string( time_stamped_name ), corrected_format, config->xmlFile->GetXML(), config->runtimeArgs );
			config->GlobalResult->GetResultParameterSet()->GetResultParameter( LLscanList[this_scan_num] )->SetScanStatus( false );
		}
	}
	return 0;
}

void SaveXML( RapidFitConfiguration* config )
{
	stringstream full_xml;

	full_xml << "<RapidFit>" << endl;
	full_xml << endl;
	if( config->generateToyXML == true )
	{
		full_xml << config->GlobalResult->GetResultParameterSet()->ToyXML();
		for( vector<PDFWithData*>::iterator toFit_i = config->pdfsAndData.begin(); toFit_i != config->pdfsAndData.end(); ++toFit_i )
		{
			vector<DataSetConfiguration*> gen_list = (*toFit_i)->GetAllDataSetConfigs();
			for( vector<DataSetConfiguration*>::iterator gen_i = gen_list.begin(); gen_i != gen_list.end(); ++gen_i )
			{
				(*gen_i)->SetSource("Foam");
			}
		}
	}
	else
	{
		full_xml << config->GlobalResult->GetResultParameterSet()->FitXML();
	}
	full_xml << endl;
	full_xml << config->theMinimiser->XML() << endl;
	full_xml << endl;
	full_xml << config->theFunction->XML() << endl;
	full_xml << endl;
	for( vector<ConstraintFunction*>::iterator const_i = config->XMLConstraints.begin(); const_i != config->XMLConstraints.end(); ++const_i )
	{
		full_xml << (*const_i)->XML() << endl;
		full_xml << endl;
	}
	full_xml << endl;
	for( vector<PDFWithData*>::iterator toFit_i = config->pdfsAndData.begin(); toFit_i != config->pdfsAndData.end(); ++toFit_i )
	{
		full_xml << (*toFit_i)->XML() << endl;
		full_xml << endl;
	}
	full_xml << "</RapidFit>" << endl;

	string fileName = ResultFormatter::GetOutputFolder();
	fileName.append("/");
	string xml_filename = "outputXMLFile";
	xml_filename.append( StringProcessing::TimeString() );
	xml_filename.append( ".xml" );
	fileName.append(xml_filename);
	ofstream output_xmlFile;
	output_xmlFile.open( fileName.c_str() );

	output_xmlFile << full_xml.str() ;

	output_xmlFile.close();

	cout << endl << "Output XML Stored in:\t" << fileName << endl << endl;
}

int PerformMainFit( RapidFitConfiguration* config )
{
	//	This is for code collapsing and to clearly outline the 'Fit' step of this file
	//		This is re-used for FC scans and forms FC Step 1

	cout << "\n\n\t\tStarting Fit to Find Global Minima!\n"<<endl;

	if( config->debug != NULL )
	{
		if( config->debug->DebugThisClass( "main" ) )
		{
			cout << "main: Starting Timer" << endl;
		}
	}

	//Do the fit to find GLOBAL MINIMA
	config->GlobalFitResult = new FitResultVector( config->argumentParameterSet->GetAllNames() );
	config->GlobalFitResult->StartStopwatch();

	if( config->debug != NULL )
	{
		if( config->debug->DebugThisClass( "main" ) )
		{
			cout << "main: Requesting External Constraints from XML" << endl;
		}
	}

	config->XMLConstraints = config->xmlFile->GetConstraints();

	config->GlobalResult = FitAssembler::DoSafeFit( config->theMinimiser, config->theFunction, config->argumentParameterSet,

			config->pdfsAndData, config->XMLConstraints, config->Force_Continue_Flag, config->OutputLevel, config->debug );

	cout << "Finished Fitting :D" << endl;

	if( config->debug != NULL )
	{
		if( config->debug->DebugThisClass( "main" ) )
		{
			cout << "main: Stopping Timer" << endl;
		}
	}

	config->GlobalFitResult->AddFitResult( config->GlobalResult );

	if( config->saveFitXML == true || config->generateToyXML == true )
	{
		SaveXML( config );
	}

	if( config->BuildConstraints )
	{
		stringstream full_xml;

		ResultParameterSet* resultParameters = config->GlobalResult->GetResultParameterSet();

		vector<string> tobeConstrained = resultParameters->GetAllFloatNames();

		full_xml << endl;

		vector<ResultParameter*> ConstrainedParams;

		for( unsigned int i=0; i< tobeConstrained.size(); ++i )
		{
			ResultParameter* thisParam = resultParameters->GetResultParameter( tobeConstrained[i] );

			full_xml << thisParam->ToyXML();
			full_xml << endl;

			ConstrainedParams.push_back( thisParam );
		}

		full_xml << endl;

		for( unsigned int i=0; i< ConstrainedParams.size(); ++i )
		{
			full_xml << ConstrainedParams[i]->ConstraintXML();
			full_xml << endl;
		}

		full_xml << endl;


		string fileName = ResultFormatter::GetOutputFolder();
		fileName.append("/");
		string xml_filename = "constraintXMLFile";
		xml_filename.append( StringProcessing::TimeString() );
		xml_filename.append( ".xml" );
		fileName.append(xml_filename);
		ofstream output_xmlFile;
		output_xmlFile.open( fileName.c_str() );

		output_xmlFile << full_xml.str() ;

		output_xmlFile.close();

		cout << endl << "Output XML Stored in:\t" << fileName << endl << endl;
	}

	if( config->GlobalResult->GetFitStatus() < 3 )
	{
		cout << "--------------------------------------------------------------" << endl;
		cout << "---------------------FIT RESULT IS NOT 3----------------------" << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << "---------If this is a Foam study, change seed and re-run------" << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << "-If your sure you want to continue employ the following flag--" << endl;
		cout << "--------------------------------------------------------------" << endl;
		cout << "--------------      \'--ForceContinue\'      -----------------" << endl;
		cout << "--------------------------------------------------------------" << endl;
		if( !config->Force_Continue_Flag  )
		{
			return -1;
		}
	}

	if( config->GlobalResult->GetFitStatus() == 2 )
	{
		cout << "ALERT ALERT ALERT " << endl << endl;
		cout << "Fit Status was 2, this means your correlation matrix was forced +ve def" << endl;
		cout << "Carrying on Regardless, but, BUYER BEWARE!" << endl << endl;
	}

	if( config->saveOneFoamDataSetFlag )
	{
		//For each dataset in the fit
		for( unsigned int i=0; i< config->pdfsAndData.size(); ++i )
		{
			//Get the phasespace from the fit
			PhaseSpaceBoundary* temp_boundary = config->pdfsAndData[i]->GetDataSet()->GetBoundary();

			//Get the PDF from the fit with the CV from the fit still passed
			IPDF* temp_pdf = config->pdfsAndData[i]->GetPDF();

			/*
			//Tell the PDF you want to copy the form of all continuous objects on the don't integrate list
			for( unsigned int j=0; j< temp_pdf->GetDoNotIntegrateList().size(); ++j )
			{
				string observableName = temp_pdf->GetDoNotIntegrateList()[j];
				IConstraint* temp_constr = temp_boundary->GetConstraint( observableName );
				if( temp_constr != NULL )
				{
					if( temp_constr->IsDiscrete() == false )
					{
						PDFConfigurator* pdf_config = new PDFConfigurator();
						pdf_config->addObservableToModel( observableName, config->pdfsAndData[i]->GetDataSet() );
						//pdf_config->setFitFunc( "pol9" );
						IPDF* observable_modeller = ClassLookUp::LookUpPDFName( "Observable_1D_distribution", pdf_config );
						temp_pdf->SetObservableDistribution( observableName, observable_modeller );
					}
				}
			}
			*/

			//Construct the configuration to give a Foam generator
			DataSetConfiguration* temp_config = new DataSetConfiguration( "Foam", config->pdfsAndData[i]->GetDataSet()->GetDataNumber(), "", vector<string>(), vector<string>(), temp_pdf );
			vector<DataSetConfiguration*> data_config; data_config.push_back( temp_config );
			PDFWithData* pdf_data_to_fit = new PDFWithData( temp_pdf, temp_boundary, data_config );

			ParameterSet* result_set = config->GlobalResult->GetResultParameterSet()->GetDummyParameterSet();
			pdf_data_to_fit->SetPhysicsParameters( result_set );

			vector<IDataSet*> quickDataSet;
			IDataSet* temp_dataSet = pdf_data_to_fit->GetDataSet();
			quickDataSet.push_back( temp_dataSet );
			string ext_dot=".";
			vector<string> temp_strings = StringProcessing::SplitString( config->saveOneDataSetFileName, *(ext_dot.c_str()) );
			TString FileName_Pre_Suffix = StringProcessing::CondenseStrings( temp_strings, 0, int(temp_strings.size() -1) );
			TString number; number+=i;
			TString real_saveOneDataSetFileName = TString( FileName_Pre_Suffix + "-" + number + ".root" );

			ResultFormatter::MakeRootDataFile( string(real_saveOneDataSetFileName.Data()), quickDataSet );
			delete temp_config;
			delete pdf_data_to_fit;
		}
	}

	if( config->WeightDataSet == true )
	{
		cout << "Weighting DataSet(s) with: " << endl;
		if( config->calcConfig == NULL )
		{
			cerr << "No PreCalculator Defined in the XML" << endl;
			cerr << "You need to add a section to your XML to tell this tool what PreCalculator to calculate event weights" << endl;
			exit(7239);
		}
		else
		{
			//		Generate Weighted Datasets
			vector<IDataSet*> WeightedDataSets;
			for( unsigned int i=0; i< config->pdfsAndData.size(); ++i )
			{
				cout << endl << "Weighting DataSet " << i+1 << endl;
				config->calculator = config->calcConfig->GetPreCalculator( config->GlobalResult );
				IDataSet* weightedDataSet = config->calculator->ProcessDataSet( config->pdfsAndData[i]->GetDataSet(), config->pdfsAndData[i]->GetPDF() );
				vector<IDataSet*> weightedDataSet_v; weightedDataSet_v.push_back( weightedDataSet );
				TString filename = config->calcConfig->GetFileName();
				filename.Append("_"); filename+=i; filename.Append(".root");
				cout << "Saving Weighted DataSet " << i << " as: " << filename << endl;
				ResultFormatter::MakeRootDataFile( filename.Data(), weightedDataSet_v );
				//delete weightedDataSet;
				//delete calculator;
				WeightedDataSets.push_back( weightedDataSet );
			}

			//		Do Something With them :D
			double sum = 0.;
			double sum_sq = 0.;
			ObservableRef WeightsName( "sWeight" );
			ObservableRef WeightsSquared( "sWeightSq" );

			Observable* testObservable = WeightedDataSets[0]->GetDataPoint( 0 )->GetObservable( WeightsName );

			if( WeightsName.GetIndex() >= 0 )
			{
				for( unsigned int i=0; i< config->pdfsAndData.size(); ++i )
				{
					for( unsigned int j=0; j< WeightedDataSets[i]->GetDataNumber(); ++j )
					{
						sum += WeightedDataSets[i]->GetDataPoint( j )->GetObservable( WeightsName )->GetValue();
						sum_sq += WeightedDataSets[i]->GetDataPoint( j )->GetObservable( WeightsSquared )->GetValue();
					}
				}
				cout << endl;
				cout << "Total of All sWeights: " << sum << endl;
				cout << "With an Error of:      " << sqrt( sum_sq ) << endl;
				cout << endl;
			}

			//		Delete them as they're local objects
			for( unsigned int i=0; i< config->pdfsAndData.size(); ++i )
			{
				while( !WeightedDataSets.empty() )
				{
					if( WeightedDataSets.back() != NULL ) delete WeightedDataSets.back();
					WeightedDataSets.pop_back();
				}
			}
		}
	}

	//ResultFormatter::FlatNTuplePullPlots( string("Global_Fit.root"), GlobalFitResult );

	cout << "\n\n\t\tFit Output:" <<endl;

	if( ( !RapidRun::isGridified() ) && ( config->Force_Continue_Flag || config->OutputLevel >= 0 ) )
	{
		//Output results
		config->makeOutput->SetInputResults( config->GlobalResult->GetResultParameterSet() );
		if( !config->doLLcontourFlag && ( !config->doFC_Flag && !config->doLLscanFlag  ) )
		{
			if( !config->disableLatexOutput )
			{
				config->makeOutput->OutputFitResult( config->GlobalFitResult->GetFitResult(0) );
			}
		}
		string fileName = ResultFormatter::GetOutputFolder();
		string fileName_ext=string( "Global_Fit_Result_"+StringProcessing::TimeString()+".root" );
		fileName.append("/"); fileName.append( fileName_ext );
		//cout << "Writing Output Tuple to:\t" << fileName << endl;
		//cout << "Output Folder:\t" << ResultFormatter::GetOutputFolder() << endl;
		//cout << "Current Folder:\t" << gSystem->pwd() << endl;
		ResultFormatter::WriteFlatNtuple( fileName, config->GlobalFitResult, config->xmlFile->GetXML(), config->runtimeArgs );
		ResultFormatter::ReviewOutput( config->GlobalResult );
	}

	//	If requested write the central value to a single file
	if( config->BurnToROOTFlag )
	{
		string time = StringProcessing::TimeString();
		string fileName=string( "Global_Fit_Result_"+time+".root" );
		cout << "Fit Output is being saved in: " << fileName << endl;
		cout << endl << "This Contains the fit result in a nTulple output, the runtime and XML used to construct the fit and the final correlation matrix" << endl;
		cout << endl;
		//ResultFormatter::WriteFlatNtuple( string( "Global_Fit_Result"+StringProcessing::TimeString()+".root" ), config->GlobalFitResult, config->xmlFile->GetXML(), config->runtimeArgs );
		ResultFormatter::WriteFlatNtuple( fileName, config->GlobalFitResult, config->xmlFile->GetXML(), config->runtimeArgs );
	}

	if( config->GOF_Flag ) {
		PDFWithData * pdfAndData = config->xmlFile->GetPDFsAndData()[0];
		pdfAndData->SetPhysicsParameters( config->xmlFile->GetFitParameters() );
		IDataSet * data = pdfAndData->GetDataSet();
		IPDF * pdf = pdfAndData->GetPDF();	(void) pdf;
		PhaseSpaceBoundary * phase = data->GetBoundary();	(void) phase;
		//GoodnessOfFit::plotUstatistic( pdf, data, phase, "ustat.pdf" );

		TH1D * pvalueHist = new TH1D("pvalues", "pvalues", 10, 0, 1);
		//double pvalue = GoodnessOfFit::gofLoop( xmlFile, theMinimiser, theFunction, argumentParameterSet, CommandLineParam, nData );
		double pvalue = GoodnessOfFit::fitDataCalculatePvalue( config->xmlFile, config->theMinimiser, config->theFunction, config->argumentParameterSet, config->GlobalResult );
		pvalueHist->Fill( pvalue );
		TFile * outputFile = new TFile("pvalues.root", "RECREATE");
		pvalueHist->Write();
		outputFile->Write();
		outputFile->Close();
		delete pvalueHist;
		delete outputFile;

	}
	return 0;
}

int PerformFCStudy( RapidFitConfiguration* config )
{
	vector<unsigned int> numberRepeatsVec;
	if( config->numberRepeatsFlag ) numberRepeatsVec.push_back( unsigned(config->numberRepeats) );

	if( config->_2DResultForFC == NULL ) return 0;

	//	Do FC scan
	vector<pair<string, string> > _2DLLscanList = config->makeOutput->Get2DScanList();
	VectoredFeldmanCousins* new_study =
		new VectoredFeldmanCousins( config->GlobalFitResult, config->_2DResultForFC, config->Nuisencemodel, config->makeOutput, config->theMinimiser, config->theFunction, config->xmlFile, config->pdfsAndData );
	if( config->numberRepeatsFlag ) new_study->SetNumRepeats( config->numberRepeats );
	new_study->DoWholeStudy( config->OutputLevel2 );
	FitResultVector* study_output = new_study->GetStudyResult();

	config->GlobalResult->GetResultParameterSet()->GetResultParameter( _2DLLscanList.back().first )->SetScanStatus( true );
	config->GlobalResult->GetResultParameterSet()->GetResultParameter( _2DLLscanList.back().second )->SetScanStatus( true );
	//      Making the assumption the user isn't running more than one of these at a time and isn't an idiot
	ResultFormatter::WriteFlatNtuple( "2DLL_FCScan.root", study_output, config->xmlFile->GetXML(), config->runtimeArgs );
	ResultFormatter::WriteFlatNtuple( "FCScan.root", study_output, config->xmlFile->GetXML(), config->runtimeArgs );
	config->GlobalResult->GetResultParameterSet()->GetResultParameter( _2DLLscanList.back().first )->SetScanStatus( false );
	config->GlobalResult->GetResultParameterSet()->GetResultParameter( _2DLLscanList.back().second )->SetScanStatus( false );

	//		STORE THE OUTPUT OF THE TOY STUDIES
	//ResultFormatter::WriteFlatNtuple( "FCOutput.root", AllFCResults );
	delete new_study;
	return 0;
}

int PerformMCStudy( RapidFitConfiguration* config )
{
	//	Process user input;
	vector<string> true_MCStepSize = StringProcessing::SplitString( config->MCStepSize, ',' );
	vector<string> true_MCStartEntry = StringProcessing::SplitString( config->MCStartEntry, ',' );

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
	MCStudy* newMCStudy = new MCStudy( config->xmlFile );

	//	Setup Toy Study
	newMCStudy->SetNumRepeats( config->numberRepeats );
	if( MCStep_int.size() != 0 ) { newMCStudy->SetStartingEntry( MCStart_int ); }
	if( MCStart_int.size() != 0 ) { newMCStudy->SetNumEvents( MCStep_int ); }


	//	Perform Toy Study
	newMCStudy->DoWholeStudy();

	ResultFormatter::WriteFlatNtuple( string( "MC_Study.root" ), newMCStudy->GetStudyResult(), config->xmlFile->GetXML(), config->runtimeArgs );

	delete newMCStudy;
	return 0;
}

int PerformToyStudy( RapidFitConfiguration* config )
{
	vector<ConstraintFunction*> XMLConstraints = config->xmlFile->GetConstraints();

	//Do the toy study
	ToyStudy* newStudy = new ToyStudy( config->theMinimiser, config->theFunction, config->argumentParameterSet, config->pdfsAndData, XMLConstraints, config->numberRepeats );

	if( config->fixedTotalToys ) newStudy->SetFixedNumberToys();
	if( config->saveAllToys ) newStudy->setSaveAllToys();

	if( config->OutputLevelSet == false ) config->OutputLevel = -999;

	newStudy->DoWholeStudy( config->OutputLevel );

	FitResultVector* fitResults = newStudy->GetStudyResult();

	//Output results
	//config->makeOutput->OutputToyResult( fitResults );
	//makeOutput->OutputFitResult( fitResults->GetFitResult(0) );

	ResultFormatter::WriteFlatNtuple( config->makeOutput->GetPullFileName(), fitResults, config->xmlFile->GetXML(), config->runtimeArgs );

	while( !XMLConstraints.empty() )
	{
		delete XMLConstraints.back();
		XMLConstraints.pop_back();
	}
	delete newStudy;

	return 0;
}

int testComponentPlot( RapidFitConfiguration* config )
{
	//Project the PDF onto the data
	PDFWithData * quickData = config->xmlFile->GetPDFsAndData()[0];
	quickData->SetPhysicsParameters( config->xmlFile->GetFitParameters() );
	IDataSet * quickDataSet = quickData->GetDataSet();
	if( config->observableName.empty() ) config->observableName = "time";
	TFile* testFile = new TFile( "testFile.root", "UPDATE" );
	ComponentPlotter * testPlotter = new ComponentPlotter( quickData->GetPDF(), quickDataSet, "testPDF", testFile, config->observableName );
	testPlotter->ProjectObservable();
	delete testPlotter;
	delete quickData;
	return 0;
}

int calculatePerEventAcceptance( RapidFitConfiguration* config )
{
	(void) config;
	PerEventAngularAcceptance* a = new PerEventAngularAcceptance("jpsikmc09_loose.root","Bu2JpsiKTuple/DecayTree", "out2.root");
	for (int iter = 1; iter <= 3; ++iter)
	{
		a->fillEffHistos( iter );
		a->loopOnReconstructedBs();
	}
	a->writeHistos();
	delete a;
	return 0;
}

int calculateAcceptanceWeightsWithSwave( RapidFitConfiguration* config )
{
	// Calculate the acceptance weights from MC
	PDFWithData * pdfAndData = config->xmlFile->GetPDFsAndData()[0];
	pdfAndData->SetPhysicsParameters( config->xmlFile->GetFitParameters() );
	IDataSet * dataSet = pdfAndData->GetDataSet();
	IPDF * pdf = pdfAndData->GetPDF();
	Mathematics::calculateAcceptanceWeightsWithSwave(dataSet, pdf);
	delete pdfAndData;
	return 0;
}

TH1D* ProjectAxis( TH3* input_histo, TString axis, TString name )
{
	TH1D* X_proj = (TH1D*)input_histo->Project3D( axis );
	X_proj->SetName( name );
	X_proj->SetTitle( name );
	X_proj->Write("",TObject::kOverwrite);
	return X_proj;
}

int calculateFitFractions( RapidFitConfiguration* config )
{
	PDFWithData * pdfAndData = config->xmlFile->GetPDFsAndData()[0];
	ParameterSet * parset = config->GlobalResult->GetResultParameterSet()->GetDummyParameterSet();
    parset->SetPhysicsParameter("BkgFraction", 0., 0., 0., 0., "BkgFraction", "");
    pdfAndData->SetPhysicsParameters( parset );//config->GlobalResult->GetResultParameterSet()->GetDummyParameterSet() );

	IDataSet * dataSet = pdfAndData->GetDataSet();
	IPDF * pdf = pdfAndData->GetPDF();

	RapidFitIntegrator * testIntegrator = new RapidFitIntegrator( pdf, true, true );
    double total_integral(0.);
	double integral(0.);
	double fraction(0.);
	double sumOfFractions(0.);

	TFile * f = TFile::Open("fitFractions.root", "RECREATE");
	TTree * tree = new TTree("tree", "tree containing fit fractions");

	vector<double> fitFractions;
	vector<string> doNotIntegrate;
	vector<string> pdfComponents = pdf->PDFComponents();
	pdfComponents = StringProcessing::MoveElementToStart( pdfComponents, "0" );
	for( unsigned int i = 0; i < pdfComponents.size(); ++i )
	{
		ComponentRef * thisRef = new ComponentRef( pdfComponents[i], "dummyObservable" );
		integral = testIntegrator->NumericallyIntegratePhaseSpace( dataSet->GetBoundary(), doNotIntegrate, thisRef );
		if ( pdfComponents[i] == "0" )
        {
            total_integral = integral;
            std::cout << "Total integral: " << total_integral << std::endl;
        }
        delete thisRef;
		fraction = integral/total_integral;
		if ( pdfComponents[i] != "0" && pdfComponents[i] != "1-Z" )
        {
            std::cout << pdfComponents[i] << "\t fraction: " << fraction << std::endl;
		    fitFractions.push_back( fraction );
		    sumOfFractions += fraction;
        }
        if ( pdfComponents[i] == "1-Z" ) std::cout << pdfComponents[i] << "\t fraction: " << 1. - fraction << std::endl;
	}
	std::cout << "Sum of fractions (not necessarily 1!): " << sumOfFractions << std::endl;

	fitFractions.push_back( sumOfFractions );
	tree->Branch("fractions", "std::vector<double>", &fitFractions);
	tree->Fill();
	f->Write();
	f->Close();

	delete testIntegrator;
	delete pdfAndData;
	return 1;
}

int calculateAcceptanceCoefficients( RapidFitConfiguration* config )
{
	PDFWithData * pdfAndData = config->xmlFile->GetPDFsAndData()[0];
	pdfAndData->SetPhysicsParameters( config->xmlFile->GetFitParameters() );
	IDataSet * dataSet = pdfAndData->GetDataSet();
	int nMCEvents = dataSet->GetDataNumber();
	IPDF * pdf = pdfAndData->GetPDF();

	return Mathematics::calculateAcceptanceCoefficients(dataSet, pdf);
}

int calculateAcceptanceWeights( RapidFitConfiguration* config )
{
	// Calculate the acceptance weights from MC
	PDFWithData * pdfAndData = config->xmlFile->GetPDFsAndData()[0];
	pdfAndData->SetPhysicsParameters( config->xmlFile->GetFitParameters() );
	IDataSet * dataSet = pdfAndData->GetDataSet();
	int nMCEvents = dataSet->GetDataNumber();
	IPDF * pdf = pdfAndData->GetPDF();

	string pdf_name = pdf->GetName();
	string Helicity_Switch("UseHelicityBasis:True");
	PDFConfigurator* tr_config = (pdf->GetConfigurator()==NULL)?new PDFConfigurator():pdf->GetConfigurator();
	PDFConfigurator* helicity_config = new PDFConfigurator( *tr_config );
	helicity_config->addConfigurationParameter( Helicity_Switch );
	IPDF* helpdf = ClassLookUp::LookUpPDFName( pdf_name, helicity_config );
	helpdf->UpdatePhysicsParameters( pdf->GetPhysicsParameters() );

	vector<double> weights = Mathematics::calculateAcceptanceWeights(dataSet, pdf);

	TString FileName("acceptance_weights_and_histos_"); FileName.Append( StringProcessing::TimeString() );
	FileName.Append(".root");
	TFile * file = TFile::Open(FileName, "RECREATE");
	TDirectory* output_file = gDirectory;
	TTree * tree = new TTree("tree", "tree containing acceptance weights and histo");
	tree->Branch("weights", "std::vector<double>", &weights);
	tree->Fill();

	// Now calculate the acceptance histograms from the data PDF/xml and MC sample
	DataSetConfiguration * dataConfig = pdfAndData->GetDataSetConfig();
	dataConfig->SetSource( "Foam" );
	PhaseSpaceBoundary* phase = new PhaseSpaceBoundary(*dataSet->GetBoundary());
	PhaseSpaceBoundary* helphase = new PhaseSpaceBoundary(*dataSet->GetBoundary());
	int factor = 10;
	int nToyEvents = nMCEvents;
	MemoryDataSet * toy = NULL;
	file->cd();
	double pi = TMath::Pi();
	output_file->cd();

	int cosThetabin = 30;
	int cosPsibin = 1;
	int phibin = 1;
	int costhetaKbin = 10;
	int costhetaLbin = 24;
	int helphibin = 5;

	TH3D * num = new TH3D("trnum", "trnum", cosPsibin, -1., 1., cosThetabin, -1., 1., phibin, -pi, pi);
	TH2D* num_psith = new TH2D( "num_tr_psitheta", "num_tr_psitheta", cosPsibin, -1., 1., cosThetabin, -1., 1. );
	TH2D* num_psiphi = new TH2D( "num_tr_psiphi", "num_tr_psiphi", cosPsibin, -1., 1., phibin, -pi, pi );
	TH2D* num_phith = new TH2D( "num_tr_phitheta", "num_tr_phi_theta", phibin, -pi, pi, cosThetabin, -1., 1. );

	TH3D * den = new TH3D("trden", "trden", cosPsibin, -1., 1., cosThetabin, -1., 1., phibin, -pi, pi);
	TH2D* den_psith = new TH2D( "den_tr_psitheta", "den_tr_psitheta", cosPsibin, -1., 1., cosThetabin, -1., 1. );
	TH2D* den_psiphi = new TH2D( "den_tr_psiphi", "den_tr_psiphi", cosPsibin, -1., 1., phibin, -pi, pi );
	TH2D* den_phith = new TH2D( "den_tr_phitheta", "den_tr_phitheta", phibin, -pi, pi, cosThetabin, -1., 1. );

	TH3D * acc = new TH3D("tracc", "tracc", cosPsibin, -1., 1., cosThetabin, -1., 1., phibin, -pi, pi);
	TH2D* acc_psith = new TH2D( "acc_tr_psitheta", "acc_tr_psitheta", cosPsibin, -1., 1., cosThetabin, -1, 1. );
	TH2D* acc_psiphi = new TH2D( "acc_tr_psiphi", "acc_tr_psiphi", cosPsibin, -1., 1., phibin, -pi, pi );
	TH2D* acc_phith = new TH2D( "acc_tr_phitheta", "acc_tr_phitheta", phibin, -pi, pi, cosThetabin, -1., 1. );

	TH3D * numh = new TH3D( "helnum", "helnum", costhetaKbin, -1., 1., costhetaLbin, -1., 1., helphibin, -pi, pi);
	TH2D* num_kl = new TH2D( "num_hel_thetaKthetaL", "num_hel_thetaKthetaL", costhetaKbin, -1., 1., costhetaLbin, -1., 1. );
	TH2D* num_kphi = new TH2D( "num_hel_thetaKphi", "num_hel_thetaKphi", costhetaKbin, -1., 1., helphibin, -pi, pi );
	TH2D* num_lphi = new TH2D( "num_hel_thetaLphi", "num_hel_thetaLphi", costhetaLbin, -1., 1., helphibin, -pi, pi );

	TH3D * denh = new TH3D( "helden", "helden", costhetaKbin, -1., 1., costhetaLbin, -1., 1., helphibin, -pi, pi);
	TH2D* den_kl = new TH2D( "den_hel_thetaKthetaL", "den_hel_thetaKthetaL", costhetaKbin, -1., 1., costhetaLbin, -1., 1. );
	TH2D* den_kphi = new TH2D( "den_hel_thetaKphi", "den_hel_thetaKphi", costhetaKbin, -1., 1., helphibin, -pi, pi );
	TH2D* den_lphi = new TH2D( "den_hel_thetaLphi", "den_hel_thetaLphi", costhetaLbin, -1., 1., helphibin, -pi, pi );

	TH3D * acch = new TH3D("helacc", "helacc", costhetaKbin, -1., 1., costhetaLbin, -1., 1., helphibin, -pi, pi);
	TH2D* acc_kl = new TH2D( "acc_hel_thetaKthetaL", "acc_hel_thetaKthetaL", costhetaKbin, -1., 1., costhetaLbin, -1., 1. );
	TH2D* acc_kphi = new TH2D( "acc_hel_thetaKphi", "acc_hel_thetaKphi", costhetaKbin, -1., 1., helphibin, -pi, pi );
	TH2D* acc_lphi = new TH2D( "acc_hel_thetaLphi", "acc_hel_thetaLphi", costhetaLbin, -1., 1., helphibin, -pi, pi );

	num->Sumw2(), numh->Sumw2();
	den->Sumw2(), denh->Sumw2();
	acc->Sumw2(), acch->Sumw2();
	double cosTheta, phi, cosPsi;
	double helcosk, helcosl, helphi;
	for( unsigned int i=0; i< (unsigned)factor; ++i )
	{
		cout << "Generating TOY DataSet: " << i+1 << " of " << factor << endl;
		toy = (MemoryDataSet*)dataConfig->MakeDataSet( phase, pdf, nToyEvents );
		output_file->cd();
		for( int j = 0; j < nToyEvents; ++j ) {
			if (j % 10000 == 0) cout << "Toy event # " << j << "\b\b\b\b\b\b\b\b\r\r\r\r\r\r\r\r\r\r";
			DataPoint * event = toy->GetDataPoint(j);
			cosPsi   = event->GetObservable("cosPsi")->GetValue();
			cosTheta = event->GetObservable("cosTheta")->GetValue();
			phi      = event->GetObservable("phi")->GetValue();
			den->Fill(cosPsi, cosTheta, phi);
			den_psith->Fill( cosPsi, cosTheta );
			den_psiphi->Fill( cosPsi, phi );
			den_phith->Fill( phi, cosTheta );
			//delete event;
		}
		delete toy;	toy = NULL;
	}
	output_file->cd();

	for( unsigned int i=0; i< (unsigned)factor; ++i )
	{
		cout << "Generating TOY DataSet: " << i+1 << " of " << factor << endl;
		toy = (MemoryDataSet*)dataConfig->MakeDataSet( helphase, helpdf, nToyEvents );
		output_file->cd();
		for( int j = 0; j < nToyEvents; ++j ) {
			if( j % 10000 == 0) cout << "Toy event # " << j << "\b\b\b\b\b\b\b\b\r\r\r\r\r\r\r\r\r\r";
			DataPoint * event = toy->GetDataPoint(j);
			helcosk  = event->GetObservable("helcosthetaK")->GetValue();
			helcosl  = event->GetObservable("helcosthetaL")->GetValue();
			helphi   = event->GetObservable("helphi")->GetValue();
			denh->Fill(helcosk, helcosl, helphi);
			den_kl->Fill( helcosk, helcosl );
			den_kphi->Fill( helcosk, helphi );
			den_lphi->Fill( helcosl, helphi );
			//delete event;
		}
		delete toy;	toy = NULL;
	}
	output_file->cd();

	for ( int i = 0; i < nMCEvents; i++ ) {
		if (i % 10000 == 0) cout << "MC event # " << i << "\b\b\b\b\b\b\b\b\r\r\r\r\r\r\r\r\r\r";
		DataPoint * event = dataSet->GetDataPoint(i);
		cosPsi   = event->GetObservable("cosPsi")->GetValue();
		cosTheta = event->GetObservable("cosTheta")->GetValue();
		phi      = event->GetObservable("phi")->GetValue();
		num->Fill(cosPsi, cosTheta, phi);
		num_psith->Fill( cosPsi, cosTheta );
		num_psiphi->Fill( cosPsi, phi );
		num_phith->Fill( phi, cosTheta );
		helcosk  = event->GetObservable("helcosthetaK")->GetValue();
		helcosl  = event->GetObservable("helcosthetaL")->GetValue();
		helphi   = event->GetObservable("helphi")->GetValue();
		numh->Fill(helcosk, helcosl, helphi);
		num_kl->Fill( helcosk, helcosl );
		num_kphi->Fill( helcosk, helphi );
		num_lphi->Fill( helcosl, helphi );
		//delete event;
	}

	acc->Divide( num, den );
	acc_psith->Divide( num_psith, den_psith );
	acc_psiphi->Divide( num_psiphi, den_psiphi );
	acc_phith->Divide( num_phith, den_phith );

	acch->Divide( numh, denh );
	acc_kl->Divide( num_kl, den_kl );
	acc_kphi->Divide( num_kphi, den_kphi );
	acc_lphi->Divide( num_lphi, den_lphi );

	TH1D* num_x = ProjectAxis( num, "x", "num_cosPsi" );
	TH1D* num_y = ProjectAxis( num, "y", "num_cosTheta" );
	TH1D* num_z = ProjectAxis( num, "z", "num_phi" );
	TH1D* den_x = ProjectAxis( den, "x", "den_cosPsi" );
	TH1D* den_y = ProjectAxis( den, "y", "den_cosTheta" );
	TH1D* den_z = ProjectAxis( den, "z", "den_phi" );
	TH1D* acc_x = new TH1D( "acc_cosPsi", "acc_cosPsi", cosPsibin, -1., 1. );
	acc_x->Divide( num_x, den_x, factor, 1. );
	acc_x->Write("",TObject::kOverwrite);
	acc_x->GetYaxis()->SetRangeUser( 0., acc_x->GetYaxis()->GetXmax() );
	TH1D* acc_y = new TH1D( "acc_cosTheta", "acc_cosTheta", cosThetabin, -1., 1. );
	acc_y->Divide( num_y, den_y, factor, 1. );
	acc_y->Write("",TObject::kOverwrite);
	acc_y->GetYaxis()->SetRangeUser( 0., acc_y->GetYaxis()->GetXmax() );
	TH1D* acc_z = new TH1D( "acc_phi", "acc_cosphi", phibin, -pi, pi );
	acc_z->Divide( num_z, den_z, factor, 1. );
	acc_z->Write("",TObject::kOverwrite);
	acc_z->GetYaxis()->SetRangeUser( 0., acc_z->GetYaxis()->GetXmax() );

	//ProjectAxis( acc, "x", "acc_cosPsi" );
	//ProjectAxis( acc, "y", "acc_cosTheta" );
	//ProjectAxis( acc, "z", "acc_phi" );

	TH1D* numh_x = ProjectAxis( numh, "x", "numh_helcosthetaK" );
	TH1D* numh_y = ProjectAxis( numh, "y", "numh_helcosthetaL" );
	TH1D* numh_z = ProjectAxis( numh, "z", "numh_helphi" );
	TH1D* denh_x = ProjectAxis( denh, "x", "denh_helcosthetaK" );
	TH1D* denh_y = ProjectAxis( denh, "y", "denh_helcosthetaL" );
	TH1D* denh_z = ProjectAxis( denh, "z", "denh_helphi" );
	TH1D* acch_x = new TH1D( "acch_helcosthetaK", "acch_helcosthetaK", costhetaKbin, -1., 1. );
	acch_x->Divide( numh_x, denh_x, factor, 1. );
	acch_x->Write("",TObject::kOverwrite);
	acch_x->GetYaxis()->SetRangeUser( 0., acch_x->GetYaxis()->GetXmax() );
	TH1D* acch_y = new TH1D( "acch_helcosthetaL", "acch_helcosthetaL", costhetaLbin, -1., 1. );
	acch_y->Divide( numh_y, denh_y, factor, 1. );
	acch_y->Write("",TObject::kOverwrite);
	acch_y->GetYaxis()->SetRangeUser( 0., acch_y->GetYaxis()->GetXmax() );
	TH1D* acch_z = new TH1D( "acc_helphi", "acc_helphi", helphibin, -pi, pi );
	acch_z->Divide( numh_z, denh_z, factor, 1. );
	acch_z->Write("",TObject::kOverwrite);
	acch_z->GetYaxis()->SetRangeUser( 0., acch_z->GetYaxis()->GetXmax() );

	num->Write("",TObject::kOverwrite);
	num_psith->Write("",TObject::kOverwrite);
	num_psiphi->Write("",TObject::kOverwrite);
	num_phith->Write("",TObject::kOverwrite);

	den->Write("",TObject::kOverwrite);
	den_psith->Write("",TObject::kOverwrite);
	den_psiphi->Write("",TObject::kOverwrite);
	den_phith->Write("",TObject::kOverwrite);

	acc->Write("",TObject::kOverwrite);
	acc_psith->Write("",TObject::kOverwrite);
	acc_psiphi->Write("",TObject::kOverwrite);
	acc_phith->Write("",TObject::kOverwrite);

	numh->Write("",TObject::kOverwrite);
	num_kl->Write("",TObject::kOverwrite);
	num_kphi->Write("",TObject::kOverwrite);
	num_lphi->Write("",TObject::kOverwrite);

	denh->Write("",TObject::kOverwrite);
	den_kl->Write("",TObject::kOverwrite);
	den_kphi->Write("",TObject::kOverwrite);
	den_lphi->Write("",TObject::kOverwrite);

	acch->Write("",TObject::kOverwrite);
	acc_kl->Write("",TObject::kOverwrite);
	acc_kphi->Write("",TObject::kOverwrite);
	acc_lphi->Write("",TObject::kOverwrite);

	tree->Write("",TObject::kOverwrite);

	//file->Write("",TObject::kOverwrite);
	file->Close();
	//delete tree;
	//delete file;
	delete pdfAndData;
	delete helpdf;
	return 0;
}

int testIntegrator( RapidFitConfiguration* config )
{
	vector<PDFWithData*> PDFinXML = config->xmlFile->GetPDFsAndData();
	for( unsigned int i=0; i< PDFinXML.size(); ++i )
	{
		//Compare numerical and analytical integration
		PDFWithData * quickData = PDFinXML[i];
		quickData->SetPhysicsParameters( config->xmlFile->GetFitParameters() );
		IDataSet * quickDataSet = quickData->GetDataSet();
		RapidFitIntegrator * testIntegrator = quickData->GetPDF()->GetPDFIntegrator();
		testIntegrator->ForceTestStatus( false );
		testIntegrator->NumericallyIntegratePhaseSpace( quickDataSet->GetBoundary(), vector<string>(), NULL );//, vector<string>(), NULL, NULL );
	}
	while( !PDFinXML.empty() )
	{
		if( PDFinXML.back() != NULL ) delete PDFinXML.back();
		PDFinXML.pop_back();
	}
	return 0;
}

int saveOneDataSet( RapidFitConfiguration* config )
{
	//Make a file containing toy data from the PDF
	cout << "Saving One DataSet: Constructing PDFs&DataSets " << endl;
	vector<PDFWithData*> quickDataGen = config->xmlFile->GetPDFsAndData();
	for( unsigned int i=0; i< quickDataGen.size(); ++i )
	{
		vector<IDataSet*> quickDataSet;
		quickDataGen[i]->SetPhysicsParameters( config->xmlFile->GetFitParameters( config->CommandLineParamvector ) );
		IDataSet* temp_dataSet = quickDataGen[i]->GetDataSet();
		quickDataSet.push_back( temp_dataSet );
		string ext_dot=".";
		vector<string> temp_strings = StringProcessing::SplitString( config->saveOneDataSetFileName, *(ext_dot.c_str()) );
		TString FileName_Pre_Suffix = StringProcessing::CondenseStrings( temp_strings, 0, int(temp_strings.size() -1) );
		TString number;

		if( quickDataGen.size() > 1){ number.Append("_"); number+=i; }

		TString real_saveOneDataSetFileName = TString( FileName_Pre_Suffix + number + ".root" );
		//quickDataSet.back()->Print();
		ResultFormatter::MakeRootDataFile( string(real_saveOneDataSetFileName.Data()), quickDataSet );
		//delete temp_dataSet;
	}
	while( !quickDataGen.empty() )
	{
		delete quickDataGen.back();
		quickDataGen.pop_back();
	}
	return 0;
}
int ConfigureRapidFit( RapidFitConfiguration* config )
{
	//Load a config file
	if( config->configFileNameFlag )
	{
		config->xmlFile = new XMLConfigReader( config->configFileName, config->XMLOverrideList, config->debug );
		//config->xmlFile->SetDebug( config->debug );
		if( config->xmlFile->IsValid() )
		{
			cout << endl << "XML config file: " << config->configFileName << " loaded Successfully!" << endl << endl;
		}
		else
		{
			cout << endl << "XML config file: " << config->configFileName << " Failed to load!" << endl << endl;
			exit(-8365);
		}
	}

	if( config->WeightDataSet == true )	{	config->calcConfig = config->xmlFile->GetPrecalculatorConfig();		}

	//	CONFIGURE RAPIDFIT PARAMETERS


	//	I am unaware of anything reasonable you can do with RapidFit without providing a driving XML
	if ( !config->configFileNameFlag )
	{
		cerr << "No Input XML defined" << endl;
		return -33;
	}

	//	This Section of Code deals with the seed value that is to be used to start whatever your studying
	if( !config->RuntimeSeed.empty() && !config->UUID_Flag )
	{
		cout << "Setting Seed At Runtime to be: " << config->RuntimeSeed[0] << endl;
		config->xmlFile->SetSeed( unsigned(config->RuntimeSeed[0]) );
	} else if ( !config->RuntimeSeed.empty() && config->UUID_Flag )	{
		//  I have used the TRandom3 code as inspiration for 'Salting the Seed' of the UUID at runtime
		cout << "Using a Random seed of UUID x RuntimeSeed to be 'unique' for all machines"<<endl;
		TUUID uid;
		UChar_t uuid[16];
		uid.GetUUID(uuid);
		config->RuntimeSeed[0] = config->RuntimeSeed[0] * ( uuid[ 2*(config->RuntimeSeed[0]%8) ]*256 +uuid[ 2*(config->RuntimeSeed[0]%8) ] );
		config->xmlFile->SetSeed( unsigned(config->RuntimeSeed[0]) );
	}

	cout << "SEED being used is:\t" << config->xmlFile->GetSeed() << endl;

	//	Command line arguments are passed and interpreted within the parser to override what is read from the XMLFile
	if( config->parameterTemplateFlag )
	{
		config->argumentParameterSet = InputParsing::MakeParameterSet( config->parameterTemplates );
	}

	//Actually configure a fit: first configure the output
	config->makeOutput = config->xmlFile->GetOutputConfiguration();
	config->makeOutput->SetDebug( config->debug );

	//	Command line arguments override the config file

	//	If the Minimiser wasn't defined at runtime consult the XML
	if( !config->theMinimiserFlag )
	{
		config->theMinimiser = config->xmlFile->GetMinimiserConfiguration();
	}
	//	If the range for a Contour Scan was provided at runtime
	if( config->defineContourFlag )
	{
		for( unsigned short int i=0; i < config->Contour_X.size(); ++i)
		{
			config->makeOutput->ClearScanList();
			config->makeOutput->Clear2DScanList();
			config->makeOutput->AddContour( config->Contour_X[i], config->Contour_Y[i] );
		}
	}
	//	If the range for a Scan was provided at runtime
	if( config->defineScanFlag )
	{
		for( unsigned short int i=0; i < config->Scan_X.size(); ++i )
		{
			config->makeOutput->ClearScanList();
			config->makeOutput->Clear2DScanList();
			config->makeOutput->AddScan( config->Scan_X[i] );
		}
	}
	//	The Minimisation Function hasn't been defined yet, request one from the XML
	if( !config->theFunctionFlag )
	{
		config->theFunction = config->xmlFile->GetFitFunctionConfiguration();
		// If weights were specified then we need to let the output plotting know
		if( config->theFunction->GetWeightsWereUsed() ) config->makeOutput->SetWeightsWereUsed( config->theFunction->GetWeightName() ) ;
	}
	//	The number of repeats hasn't been defined at runtime, look in the XML
	if( !config->numberRepeatsFlag )
	{
		config->numberRepeats = config->xmlFile->GetNumberRepeats();
	}

	//      No Parameter Template had been provided, look in the XML
	if( !config->parameterTemplateFlag )
	{
		config->argumentParameterSet = config->xmlFile->GetFitParameters( config->CommandLineParamvector );
	}

	//	PDFs used in fits and the dataset to be fit to hasn't been defined at the command line
	if( config->pdfsAndData.size() == 0 )
	{
		//	Read in from XML
		config->pdfsAndData = config->xmlFile->GetPDFsAndData();
		//	If we are performing a scan we want to check for Data Generation instances and generate/store the data in a cache for future use
		if( config->doLLscanFlag || ( config->doLLcontourFlag || config->doFC_Flag ) )
		{
			//	Loop over all ToFits containing data
			for( unsigned int pdf_num=0; pdf_num< config->pdfsAndData.size(); ++pdf_num )
			{
				//	Get the DataSet in this PDFWithData
				vector<DataSetConfiguration*> DataConfigs = config->pdfsAndData[pdf_num]->GetAllDataSetConfigs();
				//	Loop over all DataSetConfiguration objects existing within this PDFWithData
				for( unsigned int config_num=0; config_num < DataConfigs.size(); ++config_num )
				{
					//	If it's a file the data is never changed/destroyed
					if( DataConfigs[config_num]->GetSource() != "File" )
					{
						cout << "SCAN REQUESTED, GENERATING AND CACHING DATA" << endl;
						config->pdfsAndData[pdf_num]->SetPhysicsParameters( config->xmlFile->GetFitParameters() );
						vector<IDataSet*> gen_data;
						gen_data.push_back( config->pdfsAndData[pdf_num]->GetDataSet() );		//	Generate a Foam DataSet
						config->pdfsAndData[pdf_num]->AddCachedData( gen_data );		//	Cache it for the lifetime of a scan
						config->pdfsAndData[pdf_num]->SetUseCache( true );
					}
				}
			}
		}
	}

	//	Do Pulls?
	if( config->doPullsFlag )
	{
		config->makeOutput->SetPullFileName( config->pullFileName );
	}
	else
	{
		config->doPullsFlag = config->makeOutput->DoPullPlots();
	}
	//	Do Plotting
	if( config->doPlottingFlag )
	{
		config->makeOutput->MakeAllPlots( config->plotFileName );
	}

	//	FINISHED CONFIGURING RAPIDFIT PARAMETERS
	return 0;
}

int BuildTemplateXML( RapidFitConfiguration* config )
{
	cout << "Building Template XML for:\t";

	for( vector<string>::iterator template_i = config->templatePDFs.begin(); template_i != config->templatePDFs.end(); ++template_i )
	{
		cout << (*template_i) << "\t";
	}

	cout << endl;
	cout << "This is ONLY a template XML and as such I can't give you parameter configurations, his needs to be configured BY YOU" << endl;

	cout << endl;

	PDFConfigurator* empty = new PDFConfigurator();
	vector<IPDF*> requestedPDFs;
	for( vector<string>::iterator template_i = config->templatePDFs.begin(); template_i != config->templatePDFs.end(); ++template_i )
	{
		requestedPDFs.push_back( ClassLookUp::LookUpPDFName( (*template_i), empty ) );
	}

	vector<string> WantedParameters;
	vector<PhaseSpaceBoundary*> RequiredObservables;

	for( vector<IPDF*>::iterator PDF_i = requestedPDFs.begin(); PDF_i != requestedPDFs.end(); ++PDF_i )
	{
		vector<string> this_datapoint = (*PDF_i)->GetPrototypeDataPoint();
		vector<string> this_parameterset = (*PDF_i)->GetPrototypeParameterSet();

		RequiredObservables.push_back( new PhaseSpaceBoundary(this_datapoint) );

		WantedParameters = StringProcessing::CombineUniques( WantedParameters, this_parameterset );
	}

	ParameterSet* templateParameterSet = new ParameterSet( WantedParameters );

	MinimiserConfiguration* Minuit = new MinimiserConfiguration( "Minuit" );
	Minuit->SetSteps( 1000000 );
	Minuit->SetTolerance( 0.001 );
	Minuit->SetQuality( 1 );

	FitFunctionConfiguration* Function = new FitFunctionConfiguration( "NegativeLogLikelihoodThreaded" );
	Function->SetThreads( 4 );

	stringstream xml;
	xml << setprecision(10) << endl;
	xml << "<RapidFit>" << endl << endl;
	xml << templateParameterSet->XML() << endl;
	xml << endl;
	xml << Minuit->XML() << endl;
	xml << endl;
	xml << Function->XML() << endl;
	xml << endl;


	vector<PhaseSpaceBoundary*>::iterator obs_i = RequiredObservables.begin();
	vector<IPDF*>::iterator PDF_i = requestedPDFs.begin();
	for( ; PDF_i != requestedPDFs.end(); ++PDF_i, ++obs_i )
	{
		DataSetConfiguration* set_i = new DataSetConfiguration( "Foam", 10000, string(), vector<string>(), vector<string>(), ClassLookUp::CopyPDF( *PDF_i ), *obs_i );
		vector<DataSetConfiguration*> input_vec_set( 1, set_i );
		PDFWithData* toFit = new PDFWithData( ClassLookUp::CopyPDF( *PDF_i ), *obs_i, input_vec_set );

		xml << endl;
		xml << toFit->XML() << endl;
		xml << endl;

		delete toFit;
	}
	delete empty;
	while( !RequiredObservables.empty() )
	{
		if( RequiredObservables.back() != NULL ) delete RequiredObservables.back();
		RequiredObservables.pop_back();
	}
	while( !requestedPDFs.empty() )
	{
		if( requestedPDFs.back() != NULL ) delete requestedPDFs.back();
		requestedPDFs.pop_back();
	}
	delete templateParameterSet;
	delete Minuit;
	delete Function;

	xml << "</RapidFit>" << endl;

	string xml_filename = "Template_";
	xml_filename.append( StringProcessing::TimeString() );
	xml_filename.append(".xml");

	ofstream output_xmlFile;
	output_xmlFile.open( xml_filename.c_str() );
	output_xmlFile << xml.str();

	output_xmlFile.close();

	cout << endl << "Template file stored at:\t" << xml_filename << endl << endl;

	return 0;
}


