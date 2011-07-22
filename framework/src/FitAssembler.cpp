/**
  @class FitAssembler

  The intention is for this class to formalise the process of assembling the components of a fit
  Ideally it will be a set of nested static methods, starting from more and more rudimentary components

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

//	RapidFit Headers
#include "FitAssembler.h"
#include "FitResult.h"
#include "ClassLookUp.h"
#include "ScanParam.h"
#include "FitResultVector.h"
#include "ResultFormatter.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using namespace std;


void FitAssembler::SafeMinimise( IMinimiser* Minimiser )
{
	try
	{
		// Try a Fit, it it converges, continue to elsewhere in the program
		Minimiser->Minimise();
	}
	//  If it didn't fit tell the user why!
	catch( int e)
	{
		if ( e == 10 )
		{
			cerr << "\nCaught exception : fit failed for these parameters..." << endl; 
		}
		else if ( e == 13 )
		{
			cerr << "\nIntegration Error: Fit Failed..." << endl;
		}
	}
	catch (...)
	{
		cerr << "\n\n\n\t\t\tCaught Unknown Exception, THIS IS SERIOUS!!!\n\n\n" << endl;
	}
}

//The final stage - do the minimisation
FitResult * FitAssembler::DoFit( IMinimiser * Minimiser, FitFunction * TheFunction )
{
	Minimiser->SetupFit( TheFunction );

	SafeMinimise( Minimiser );

	return Minimiser->GetFitResult();
}

//Create the minimiser and fit function
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, PhysicsBottle * Bottle )
{
	IMinimiser * minimiser = MinimiserConfig->GetMinimiser( int(Bottle->GetParameterSet()->GetAllNames().size()) );
	FitFunction * theFunction = FunctionConfig->GetFitFunction( Bottle );
	//theFunction->SetPhysicsBottle(Bottle);

	FitResult* result = DoFit( minimiser, theFunction );

	return result;
}

//Create the physics bottle
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, vector<ParameterSet*> BottleParameters, const vector< PDFWithData* > BottleData, const vector< ConstraintFunction* > BottleConstraints )
{
	PhysicsBottle * bottle = new PhysicsBottle( BottleParameters.back() );

	//Fill the bottle - data generation occurs in this step
	for ( unsigned int resultIndex = 0; resultIndex < BottleData.size(); ++resultIndex )
	{
		BottleData[resultIndex]->SetPhysicsParameters(BottleParameters);
		IPDF* Requested_PDF = BottleData[resultIndex]->GetPDF();
		IDataSet* Requested_DataSet = BottleData[resultIndex]->GetDataSet();

		//		Requested_DataSet->SortBy("time");

		bottle->AddResult( Requested_PDF, Requested_DataSet );
	}

	//Add the constraints
	for ( unsigned int constraintIndex = 0; constraintIndex < BottleConstraints.size(); ++constraintIndex )
	{
		bottle->AddConstraint( BottleConstraints[constraintIndex] );
	}

	bottle->Finalise();

	FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );

	//delete bottle;
	return result;
}

//Create the physics bottle with pre-made data
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, vector< ParameterSet* > BottleParameters, const vector< IPDF* > AllPDFs, const vector< IDataSet* > AllData, const vector< ConstraintFunction* > BottleConstraints )
{
	if ( AllPDFs.size() == AllData.size() )
	{
		PhysicsBottle * bottle = new PhysicsBottle(BottleParameters.back());

		//Fill the bottle - data already generated
		for ( unsigned int resultIndex = 0; resultIndex < AllData.size(); ++resultIndex )
		{
			AllPDFs[resultIndex]->SetPhysicsParameters(BottleParameters.back());
			bottle->AddResult( AllPDFs[resultIndex], AllData[resultIndex] );
		}

		//Add the constraints
		for ( unsigned int constraintIndex = 0; constraintIndex < BottleConstraints.size(); ++constraintIndex )
		{
			bottle->AddConstraint( BottleConstraints[constraintIndex] );
		}

		bottle->Finalise();
		FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );

		//delete bottle;
		return result;
	}
	else
	{
		cerr << "Mismatched number of PDFs and DataSets" << endl;
		exit(1);
	}
}

void FitAssembler::CheckParameterSet( FitResult* ReturnableFitResult, vector< ParameterSet* > BottleParameters )
{
	vector<string> already_found = ReturnableFitResult->GetResultParameterSet()->GetAllNames();

	for( unsigned int i=0; i< BottleParameters.back()->GetAllNames().size(); ++i )
	{
		int found = StringProcessing::VectorContains( &already_found, &(BottleParameters.back()->GetAllNames()[i]) );

		//	There was something in the ParameterSet not in the FitResult, i.e. an unclaimed object which can't have changed during the fit
		if( found == -1 )
		{
			double Value = BottleParameters.back()->GetPhysicsParameter( BottleParameters.back()->GetAllNames()[i] )->GetValue();
			double OriginalValue = BottleParameters.back()->GetPhysicsParameter( BottleParameters.back()->GetAllNames()[i] )->GetValue();
			double Error = 0;
			double Minimum = BottleParameters.back()->GetPhysicsParameter( BottleParameters.back()->GetAllNames()[i] )->GetMinimum();
			double Maximum = BottleParameters.back()->GetPhysicsParameter( BottleParameters.back()->GetAllNames()[i] )->GetMaximum();
			double StepSize = BottleParameters.back()->GetPhysicsParameter( BottleParameters.back()->GetAllNames()[i] )->GetStepSize();
			string Type = BottleParameters.back()->GetPhysicsParameter( BottleParameters.back()->GetAllNames()[i] )->GetType();
			string Unit = BottleParameters.back()->GetPhysicsParameter( BottleParameters.back()->GetAllNames()[i] )->GetUnit();
			bool added = ReturnableFitResult->GetResultParameterSet()->ForceNewResultParameter( BottleParameters.back()->GetAllNames()[i],  Value, OriginalValue, Error, Minimum, Maximum, StepSize, Type, Unit );
			if( !added )
			{
				cerr << "Error finalizing FitResultVector Object" << endl << endl;
				exit(-984);
			}
		}
	}
}

//  Perform a safer fit which is gauranteed to return something which you can use :D
FitResult * FitAssembler::DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, vector< ParameterSet* > BottleParameters, const vector< PDFWithData* > BottleData, const vector< ConstraintFunction* > BottleConstraints, const int OutputLevel )
{
	FitResult* ReturnableFitResult=NULL;

	if( FunctionConfig->GetStrategy() == "Petes" )
	{
		ReturnableFitResult = Petes_DoSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, OutputLevel );
	}
	else
	{
		ReturnableFitResult = DoSingleSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, OutputLevel );
	}

	return ReturnableFitResult;
}

FitResult * FitAssembler::DoSingleSafeFit(MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, vector< ParameterSet* > BottleParameters, const vector< PDFWithData* > BottleData, const vector< ConstraintFunction* > BottleConstraints, const int OutputLevel )
{
	streambuf *nullbuf=NULL, *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
	ofstream filestr;
	filestr.open ("/dev/null");
	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		nullbuf = filestr.rdbuf();
		cout.rdbuf(nullbuf);
		cerr.rdbuf(nullbuf);
		clog.rdbuf(nullbuf);
	}

	FitResult* ReturnableFitResult=NULL;
	MinimiserConfig->SetOutputLevel( OutputLevel );
	vector<string> other_params = BottleParameters.back()->GetAllFloatNames();      //      This better at least contain all in prototypeparamset!!!
	vector<double> truth;
	for( unsigned short int j=0; j < other_params.size(); ++j )
	{
		truth.push_back( BottleParameters.back()->GetPhysicsParameter( other_params[j] )->GetValue() );
	}

	// Try a Fit, it it converges, continue to elsewhere in the program
	ReturnableFitResult = FitAssembler::DoFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints );

	//      Reset Std Output Streams
	if( OutputLevel <= -1 )
	{
		cout.rdbuf(cout_bak);
		cerr.rdbuf(cerr_bak);
		clog.rdbuf(clog_bak);
	}

	if( ReturnableFitResult->GetFitStatus() != 3 )
	{
		cerr << "\n\n\t\tFit Did NOT Converge Correctly, CHECK YOUR RESULTS!\n\n";
		for( unsigned short int j=0; j < other_params.size(); ++j )
		{
			BottleParameters.back()->GetPhysicsParameter( other_params[j] )->SetValue( truth[j] );
		}
		int status = -1;
		vector<string> NewNamesList = BottleParameters.back()->GetAllNames();
		ResultParameterSet* DummyFitResults = new ResultParameterSet( NewNamesList );
		PhysicsBottle* Bad_Bottle = new PhysicsBottle( BottleParameters.back() );
		ReturnableFitResult = new FitResult( LLSCAN_FIT_FAILURE_VALUE, DummyFitResults, status, Bad_Bottle );
	}

	CheckParameterSet( ReturnableFitResult, BottleParameters );

	return ReturnableFitResult;
}


FitResult * FitAssembler::Petes_DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, vector< ParameterSet* > BottleParameters, const vector< PDFWithData* > BottleData, const vector< ConstraintFunction* > BottleConstraints, const int OutputLevel )
{
	cout << endl << "******* Result of Petes Double fit strategy*********" << endl ;
	cout << "Starting Fit1:" << endl;
	// Normal fit
	FitResult* res0 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, OutputLevel ) ;
	double LLmin0 = res0->GetMinimumValue() ;
	bool good_result_0 = res0->GetFitStatus() == 3;
	cout << "Finished Fit1." << endl;
	if( !good_result_0  ) cout << "Fit-1 failed" << endl;

	// Conjugate fit
	double deltaPara = BottleParameters.back()->GetPhysicsParameter( string("delta_para") )->GetBlindedValue() ;
	double deltaPerp = BottleParameters.back()->GetPhysicsParameter( string("delta_perp") )->GetBlindedValue() ;
	BottleParameters.back()->GetPhysicsParameter( string("delta_para") )->SetBlindedValue( -deltaPara ) ;
	BottleParameters.back()->GetPhysicsParameter( string("delta_perp") )->SetBlindedValue( 3.14159-deltaPerp ) ;
	cout << endl << "Starting Fit2:" << endl;
	FitResult* res1 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, OutputLevel ) ;
	double LLmin1 = res1->GetMinimumValue() ;
	bool good_result_1 = res1->GetFitStatus() == 3;

	// Apply strategy 
	FitResult * res ;
	bool swapped = false ;
	if( !good_result_0 && !good_result_1 ){ 
		res = res0;     //both fail so it doesnt matter 
	}
	else if( good_result_0 && !good_result_1 ){     
		res = res0;     
	}
	else if( !good_result_0 && good_result_1 ){     
		res = res1;     
		swapped = true;
	}
	else if( LLmin0 < LLmin1+0.001 ){ 
		res = res0 ;
	}
	else { 
		res = res1 ;
		swapped = true;
	}

	//MinimiserConfig->SetOutputLevel( 0 );
	cout << setprecision(9) ;
	if( !good_result_1  ) cout << "Fit-2 failed" << endl;
	cout << "Finished Fit2." << endl;
	cout << "The 2 LLs were "  << LLmin0 << "   <=>   "  << LLmin1 << endl ;
	if( swapped  ) cout << "Swapped to conjugate fit result " << endl;
	cout << "******* Result of Petes Double fit strategy*********" << endl ;
	//MinimiserConfig->SetOutputLevel( OutputLevel );

	// Set input parameters back to what they were 
	BottleParameters.back()->GetPhysicsParameter( string("delta_para") )->SetBlindedValue( deltaPara ) ;  
	BottleParameters.back()->GetPhysicsParameter( string("delta_perp") )->SetBlindedValue( deltaPerp ) ;  

	return res ;

}

//  Interface for internal calls
void FitAssembler::DoScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, vector< ParameterSet* > BottleParameters, const vector< PDFWithData* > BottleData, const vector< ConstraintFunction* > BottleConstraints, ScanParam* Wanted_Param, FitResultVector* output_interface, const int OutputLevel )
{

	double uplim = Wanted_Param->GetMax();
	double lolim = Wanted_Param->GetMin();
	double npoints = Wanted_Param->GetPoints();
	string scanName = Wanted_Param->GetName();

	//	cout << "Performing Scan for the parameter " << scanName << endl ;

	// Get a pointer to the physics parameter to be scanned and fix it	
	// CAREFUL:  this must be reset as it was at the end.
	PhysicsParameter * scanParameter = BottleParameters.back()->GetPhysicsParameter(scanName);
	double originalValue = scanParameter->GetBlindedValue( ) ;
	string originalType = scanParameter->GetType( ) ;
	scanParameter->SetType( "Fixed" ) ;

	// Need to set up a loop , fixing the scan parameter at each point
	double deltaScan;

	for( int si=0; si<int(npoints); ++si)
	{
		cout << "\n\nSINGLE SCAN NUMBER\t\t" << si+1 << "\t\tOF\t\t" <<int(npoints)<< endl<<endl;
		// Set scan parameter value
		if( int(npoints)!=1 ) deltaScan = (uplim-lolim) / (npoints-1.) ;
		else deltaScan=0;
		double scanVal = lolim+deltaScan*si;
		scanParameter->SetBlindedValue( scanVal ) ;

		output_interface->StartStopwatch();

		//	Use the SafeFit as this always returns something when a PDF has been written to throw not exit
		//	Do a scan point fit
		FitResult * scanStepResult = FitAssembler::DoSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, OutputLevel );

		cout << "Fit Finished!\n" <<endl;
		//  THIS IS ALWAYS TRUE BY DEFINITION OF THE SCAN

		string name = Wanted_Param->GetName();
		double StepSize = BottleParameters.back()->GetPhysicsParameter( name )->GetStepSize();
		string type = BottleParameters.back()->GetPhysicsParameter( name )->GetType();
		string unit = BottleParameters.back()->GetPhysicsParameter( name )->GetUnit();
		scanStepResult->GetResultParameterSet()->SetResultParameter( name, scanVal, 0, 0., scanVal, scanVal, StepSize, type, unit );

		vector<string> Fixed_List = BottleParameters.back()->GetAllFixedNames();
		vector<string> Fit_List = scanStepResult->GetResultParameterSet()->GetAllNames();
		for( unsigned short int i=-1; i < Fixed_List.size() ; ++i )
		{
			bool found=false;
			for( unsigned short int j=0; j < Fit_List.size(); ++j )
			{
				if( Fit_List[j] == Fixed_List[i] )
				{
					found = true;
				}
			}
			if( !found )
			{
				double fixed_StepSize = BottleParameters.back()->GetPhysicsParameter( Fixed_List[i] )->GetStepSize();
				string fixed_type = BottleParameters.back()->GetPhysicsParameter( Fixed_List[i] )->GetType();
				string fixed_unit = BottleParameters.back()->GetPhysicsParameter( Fixed_List[i] )->GetUnit();
				double fixed_value = BottleParameters.back()->GetPhysicsParameter( Fixed_List[i] )->GetValue();
				scanStepResult->GetResultParameterSet()->ForceNewResultParameter( Fixed_List[i], fixed_value, fixed_value, 0, fixed_value, fixed_value, fixed_StepSize, fixed_type, fixed_unit );
			}
		}

		ResultFormatter::ReviewOutput( scanStepResult );

		output_interface->AddFitResult( scanStepResult );
	}

	//Reset the parameter as it was
	scanParameter->SetType( originalType ) ;
	scanParameter->SetBlindedValue( originalValue ) ;

}

//  Interface for internal calls
void FitAssembler::DoScan2D( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, vector< ParameterSet* > BottleParameters, const vector< PDFWithData* > BottleData, const vector< ConstraintFunction* > BottleConstraints, const pair<ScanParam*, ScanParam*> Param_Set, vector<FitResultVector*>* output_interface, const int OutputLevel )
{

	//	vector<string> namez = BottleParameters.back()->GetAllNames();
	vector<string> result_names = BottleParameters.back()->GetAllNames();
	double uplim = Param_Set.first->GetMax();
	double lolim = Param_Set.first->GetMin();
	double npoints = Param_Set.first->GetPoints();

	string scanName = Param_Set.first->GetName();
	string scanName2 = Param_Set.second->GetName();


	// Get a pointer to the physics parameter to be scanned and fix it
	// CAREFUL:  this must be reset as it was at the end.
	PhysicsParameter * scanParameter = BottleParameters.back()->GetPhysicsParameter(scanName);
	double originalValue = scanParameter->GetBlindedValue( );
	string originalType = scanParameter->GetType( );
	scanParameter->SetType( "Fixed" );

	// Need to set up a loop , fixing the scan parameter at each point

	double deltaScan;
	if( int(npoints) !=1 ) deltaScan = (uplim-lolim) / (npoints-1.) ;
	else deltaScan=0.;

	for( int si=0; si < int(npoints); ++si) {

		cout << "\n\n2DSCAN OUTER NUMBER\t\t" << si+1 << "\t\tOF\t\t" << int(npoints) <<endl<<endl;
		FitResultVector* Returnable_Result = new FitResultVector( result_names );

		// Set scan parameter value
		double scanVal = lolim + si*deltaScan;
		scanParameter->SetBlindedValue( scanVal ) ;

		// Do a scan point fit
		FitAssembler::DoScan( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, Param_Set.second, Returnable_Result, OutputLevel );

		//  THIS IS ALWAYS TRUE BY DEFINITION OF THE SCAN
		string name = Param_Set.first->GetName();
		double step = BottleParameters.back()->GetPhysicsParameter( name )->GetStepSize();
		string type = BottleParameters.back()->GetPhysicsParameter( name )->GetType();
		string unit = BottleParameters.back()->GetPhysicsParameter( name )->GetUnit();

		for( short int i=0; i < Returnable_Result->NumberResults(); ++i )
		{
			Returnable_Result->GetFitResult( i )->GetResultParameterSet()->SetResultParameter( name, scanVal, 0, 0.0, scanVal, scanVal, step, type, unit );
		}

		output_interface->push_back( Returnable_Result );
	}

	//Reset the parameter as it was
	scanParameter->SetType( originalType ) ;
	scanParameter->SetBlindedValue( originalValue ) ;

}

// Interface for external calls
vector<FitResultVector*> FitAssembler::ContourScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, vector< ParameterSet* > BottleParameters, const vector< PDFWithData* > BottleData, const vector< ConstraintFunction* > BottleConstraints, OutputConfiguration* OutputConfig, const string scanName, const string scanName2, const int OutputLevel )
{
	vector<FitResultVector*>* Returnable_Result = new vector<FitResultVector*>;

	pair< ScanParam*, ScanParam* > Param_Set = OutputConfig->Get2DScanParams( scanName, scanName2 );

	DoScan2D( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, Param_Set, Returnable_Result, OutputLevel );

	return *Returnable_Result;
}

//  Interface for external calls
FitResultVector* FitAssembler::SingleScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, vector< ParameterSet*> BottleParameters, vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, OutputConfiguration* OutputConfig, const string scanName, const int OutputLevel )
{
	FitResultVector* Returnable_Result = new FitResultVector( BottleParameters.back()->GetAllNames() );

	ScanParam* local_param = OutputConfig->GetScanParam( scanName );

	DoScan( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, local_param, Returnable_Result, OutputLevel );

	return Returnable_Result;
}

//	This is the Main loop of Code which performs the Numerical work needed for the Feldman-Cousins Fit in RapidFit
//	The output of this code is to be analysed by Conor's plotting tool
//	I make no claims as to how FC plots are extracted from this I will just tell you what this code Does
//
//	Step1:		Find Global Minima
//	Step2:		Find Minima at point i,j
//	Step3:		Setup Framework to, and Generate Data for toy study with control paramater for i,j Fixed
//			You can either set nuisence parameters to their central values when fit at Step 1 or 2
//			At the time of writing the nuisence parameters are still under investigation
//	Step4:		Repeat Step 3 as many times as requested, or do 100 toys as default
//	Step5:		Format output into standard block format. This gives us 2 things.
FitResultVector* FitAssembler::FeldmanCousins( FitResultVector* GlobalResult, FitResultVector* _2DResultForFC, vector<unsigned int> numberRepeats, const unsigned int NuisenceModel, const bool FC_Debug_Flag, OutputConfiguration* makeOutput, MinimiserConfiguration* theMinimiser, FitFunctionConfiguration* theFunction, XMLConfigReader* xmlFile, vector< PDFWithData* > pdfsAndData, const int OutputLevel )
{
	double Random_Seed = pdfsAndData[0]->GetPDF()->GetRandomFunction()->Rndm();
	FitResult* GlobalFitResult = GlobalResult->GetFitResult( 0 );
	vector<FitResultVector*> AllResults;

	//		Want to loop over all points on a 2D Plot that have been allocated to this instance
	for( int iFC=0; iFC < _2DResultForFC->NumberResults(); ++iFC )
	{
		TRandom3* new_rand = new TRandom3(UInt_t(Random_Seed));
		//		GET INPUT Data from fit Results
		//
		//  Get a ParameterSet that contains the input parameters from the output of a fit
		vector<pair<string, string> > _2DLLscanList = makeOutput->Get2DScanList();
		string name1 = _2DLLscanList[0].first;
		string name2 = _2DLLscanList[0].second;

		vector<ParameterSet*> InputParamSet;
		vector<ParameterSet*> InputFreeSet;
		vector<ParameterSet*> ControlParamSet;

		//	True by definition:
		double lim1 = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetResultParameter( name1 )->GetValue();
		double lim2 = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetResultParameter( name2 )->GetValue();

		if( NuisenceModel == 1 )
		{
			//		Use the inputs from Step 1
			InputParamSet.push_back( GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet() );
			InputFreeSet.push_back( GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet() );
			ControlParamSet.push_back( GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet() );
		} else if( NuisenceModel == 2 ){
			//		Use the inputs from Step 2
			InputParamSet.push_back( _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet() );
			InputFreeSet.push_back( _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet() );
			ControlParamSet.push_back( _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet() );
		}


		//		Just for clarity
		//		also when using values from Step1 these are not set correctly for the study
		ControlParamSet.back()->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
		ControlParamSet.back()->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
		InputParamSet.back()->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
		InputParamSet.back()->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
		InputParamSet.back()->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
		InputParamSet.back()->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
		InputParamSet.back()->GetPhysicsParameter( name1 )->SetType( "Fixed" );
		InputParamSet.back()->GetPhysicsParameter( name2 )->SetType( "Fixed" );
		InputFreeSet.back()->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
		InputFreeSet.back()->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
		InputFreeSet.back()->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
		InputFreeSet.back()->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
		InputFreeSet.back()->GetPhysicsParameter( name1 )->SetType( "Free" );
		InputFreeSet.back()->GetPhysicsParameter( name2 )->SetType( "Free" );

		cout << "Here" << endl;

		//	Collect all of the relevent Data from the XML
		//	Note: most of these had to be written for FCscans
		MinimiserConfiguration * ToyStudyMinimiser = theMinimiser;
		theMinimiser->SetOutputLevel(-999);
		FitFunctionConfiguration* ToyStudyFunction = theFunction;
		vector<vector<IPrecalculator*> > ToyPrecalculators = xmlFile->GetPrecalculators();
		vector<PhaseSpaceBoundary*> PhaseSpaceForToys = xmlFile->GetPhaseSpaceBoundaries();
		vector<ConstraintFunction*> ConstraintsForToys = xmlFile->GetConstraints();

		//	Read the number of events from the existing DataSets
		//	I think there needs be an equivalent function drafted to use sWeights
		vector<unsigned int> EventsPerPDF;
		vector<double> sweight_error;

		for( unsigned short int pdf_num=0; pdf_num < pdfsAndData.size(); ++pdf_num )
		{
			string sWeightObs="Nsig_sw";
			bool sWeighted=false;
			vector<string> phase_names = PhaseSpaceForToys[pdf_num]->GetAllNames();
			int sweight_num = StringProcessing::VectorContains( &phase_names, &sWeightObs );
			if( sweight_num != -1 ) sWeighted=true;
			if( sWeighted )
			{
				for(unsigned short int i=0; i< PhaseSpaceForToys.size(); ++i )
				{
					vector<double> new_constraint(1,1.0);
					PhaseSpaceForToys[i]->SetConstraint( "Nsig_sw", new_constraint, " " );
				}
				IDataSet* file_input = pdfsAndData[pdf_num]->GetDataSet();
				int point_number = file_input->GetDataNumber();
				double data_num=0; double data_err=0;
				for( int i=0; i < point_number; ++i)
				{
					double sweight_val = file_input->GetDataPoint( i )->GetObservable(sWeightObs)->GetValue();
					data_num+=sweight_val;
					data_err+=sweight_val*sweight_val;
				}
				data_err=sqrt(data_err);
				sweight_error.push_back( unsigned(int(data_err)) );
				EventsPerPDF.push_back( unsigned(int(data_num)) );
			}
			else	{
				EventsPerPDF.push_back( unsigned(pdfsAndData[pdf_num]->GetDataSet()->GetDataNumber()) );
				sweight_error.push_back( 0. );
			}
		}



		//		Number of datasets wanted i.e. Number of Toys
		unsigned int wanted_number_of_toys = 100;
		if( !numberRepeats.empty() ) wanted_number_of_toys = unsigned(numberRepeats[0]);



		//		GENERATE DATA
		//	We want to now generate a new DataSet for a Toy Study
		//	I choose a Foam DataSet with no Cuts or extra arguments as these are not required
		//
		//	Data is now only generated and held in memory for as long as it's needed as this
		//	has a lighter memory footprint
		//
		//	Generate Once, Fit twice
		//	This was a bit of a pain to code up.
		//	Although I like the idea of coding it up into a GenerateToyData object in future to avoid this

		vector<PDFWithData*> PDFsWithDataForToys;
		vector<vector<IDataSet*> > Memory_Data;
		Memory_Data.resize( pdfsAndData.size() );

		//	NB:
		//	Generating data the first time is more complex as we need to construct the correct PDFWithData
		//	The object setup here will be used later in the code as a template with cached data being replaced
		//	This reduces the amount of code and makes this procedure easier to understand as a whole

		cout << "\n\n\t\tGenerating Data For First Toy\n" << endl;
		//	We may have multiple PDFs in the XML
		for( unsigned short int pdf_num=0; pdf_num < pdfsAndData.size(); ++pdf_num )
		{
			//	Generate the Data for this Study using the given the XML PhaseSpace and PDF
			IPDF* PDF_from_XML = pdfsAndData[pdf_num]->GetPDF();

			//	Create a DataSetConfiguration to make the datasets
			vector<string> empty_args;
			vector<string> empty_arg_names;
			unsigned int wanted_events = EventsPerPDF[pdf_num];
			if( sweight_error[pdf_num] > 0. )  wanted_events *= unsigned(int(new_rand->Gaus()*sweight_error[pdf_num]));
			DataSetConfiguration* Toy_Foam_DataSet= new DataSetConfiguration( "Foam", wanted_events, "", empty_args, empty_arg_names, PDF_from_XML );

			//	Set the Input Physics Parameters to be as defined above ControlSet
			//	This is to seperate changing bottles for Physics from that used at Generation
			//	as we want the data to have a wider scope than 1 fit with the data
			Toy_Foam_DataSet->SetPhysicsParameters( ControlParamSet );

			//	Let the user know what's going on
			cout << "generating data for pdf: " << pdf_num << "\n";

			//	Make the data
			IDataSet* new_dataset = Toy_Foam_DataSet->MakeDataSet( PhaseSpaceForToys[pdf_num], PDF_from_XML );
			//	Store the Data Object in Memory
			Memory_Data[pdf_num].push_back( new_dataset );

			//	Store the information for each PDFWithData
			vector<DataSetConfiguration*> DataSetConfigForToys;
			DataSetConfigForToys.push_back( Toy_Foam_DataSet );
			PDFWithData* ToyPDFWithData = new PDFWithData( PDF_from_XML, PhaseSpaceForToys[pdf_num], DataSetConfigForToys, ToyPrecalculators[pdf_num] );

			//	Store this PDFWithData
			ToyPDFWithData->SetPhysicsParameters( InputParamSet );
			PDFsWithDataForToys.push_back( ToyPDFWithData );
		}

		//		Result Vectors for the Fit for clarity and Output Formatting
		FitResultVector* study1Results = new FitResultVector( GlobalFitResult->GetResultParameterSet()->GetAllNames() );
		FitResultVector* study2Results = new FitResultVector( GlobalFitResult->GetResultParameterSet()->GetAllNames() );

		//	This Forms the main sequence of the fit

		cout << "\n\n\t\tPerforming Fits to Toys in FC\n"<<endl;
		//		Perform Fits
		//	This will record ONLY Data from working fits i.e. FitStatus==3 unless the Debug flag is on
		for( unsigned short int dataset_num=0; dataset_num < wanted_number_of_toys; ++dataset_num )
		{

			bool toy_failed = false;

			FitResult* fit1Result = NULL;
			FitResult* fit2Result = NULL;

			vector<ParameterSet*> LocalInputFreeSet;
			vector<ParameterSet*> LocalInputFixedSet;

			if( NuisenceModel == 1 )	{
				//	Assuming Nuisence Parameters set to Global Minima
				//	We need to set this for EVERY FIT in order to have correct generation/pull values
				LocalInputFreeSet.push_back( GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet() );
				LocalInputFixedSet.push_back( GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet() );
			} else if( NuisenceModel == 2 )	{
				//	Assuming Nuisence Parameters set to Local Minima
				//
				LocalInputFreeSet.push_back( _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet() );
				LocalInputFixedSet.push_back( _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet() );
			} else if( NuisenceModel == 3 )	{
				//	Generate a ParameterSet with random numbers based on
				//	Fit Results and assign it to the ControlParamSet which is used for generation
				//	Also need to Setup Local ParameterSets as with standard case
			}

			ControlParamSet.back()->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
			ControlParamSet.back()->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
			LocalInputFreeSet.back()->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
			LocalInputFreeSet.back()->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
			LocalInputFreeSet.back()->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
			LocalInputFreeSet.back()->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
			LocalInputFreeSet.back()->GetPhysicsParameter( name1 )->SetType( "Free" );
			LocalInputFreeSet.back()->GetPhysicsParameter( name2 )->SetType( "Free" );

			LocalInputFixedSet.back()->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
			LocalInputFixedSet.back()->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
			LocalInputFixedSet.back()->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
			LocalInputFixedSet.back()->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
			LocalInputFixedSet.back()->GetPhysicsParameter( name1 )->SetType( "Fixed" );
			LocalInputFixedSet.back()->GetPhysicsParameter( name2 )->SetType( "Fixed" );

			//	We need to set some factors before we perform the fit
			for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); ++pdf_num )
			{
				vector<IDataSet*> wanted_set;
				wanted_set.push_back( Memory_Data[pdf_num][dataset_num] );
				PDFsWithDataForToys[pdf_num]->AddCachedData( wanted_set );
				PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( LocalInputFreeSet );
			}

			cout << "\n\n\t\tPerforming Fit To Toy: "<< (dataset_num+1) <<" of " << wanted_number_of_toys << endl<<endl;

			streambuf *nullbuf=NULL, *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
			ofstream filestr;
			filestr.open ("/dev/null");
			cout << "\n\tStarting Fit:" <<endl;
			//	If the user wanted silence we point the Std Output Streams to /dev/null
			if( OutputLevel <= -1 )
			{
				cout_bak = cout.rdbuf();
				cerr_bak = cerr.rdbuf();
				clog_bak = clog.rdbuf();
				nullbuf = filestr.rdbuf();
				cout.rdbuf(nullbuf);
				cerr.rdbuf(nullbuf);
				clog.rdbuf(nullbuf);
			}
			//	Fit once with control parameters Free
			fit1Result = FitAssembler::DoSafeFit( ToyStudyMinimiser, ToyStudyFunction, LocalInputFreeSet, PDFsWithDataForToys, ConstraintsForToys, OutputLevel );
			//	Reset Std Output Streams
			if( OutputLevel <= -1 )
			{
				cout.rdbuf(cout_bak);
				cerr.rdbuf(cerr_bak);
				clog.rdbuf(clog_bak);
			}
			cout << "\tFit Finished!\n" <<endl;

			//	Only Fit again to this dataset if it fits well with +2 dof
			//	This has the obvious savings in CPU resources
			//	We may want to keep this information so leaving it configurable.
			if( ( fit1Result->GetFitStatus() == 3 ) || FC_Debug_Flag )
			{
				if( FC_Debug_Flag )	cout << "\n\nYou are aware your requesting all output?\n"<<endl;

				ResultFormatter::ReviewOutput( fit1Result );

				//  Fit secondGlobalResult with control parameters Fixed
				for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); ++pdf_num )
				{
					vector<IDataSet*> wanted_set;
					wanted_set.push_back( Memory_Data[pdf_num][dataset_num] );
					PDFsWithDataForToys[pdf_num]->AddCachedData( wanted_set );
					PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( LocalInputFixedSet );
				}

				cout << "\n\t\tFirst Fit Successful, Performing the Second Fit " << (dataset_num+1) << " of " << wanted_number_of_toys <<endl;

				cout << "\n\tStarting Fit:" <<endl;
				//	If the user wanted silence we point the Std Output Streams to /dev/null
				if( OutputLevel <= -1 )
				{
					cout_bak = cout.rdbuf();
					cerr_bak = cerr.rdbuf();
					clog_bak = clog.rdbuf();
					nullbuf = filestr.rdbuf();
					cout.rdbuf(nullbuf);
					cerr.rdbuf(nullbuf);
					clog.rdbuf(nullbuf);
				}
				//	Use the SafeFit as this always returns something when a PDF has been written to throw not exit
				fit2Result = FitAssembler::DoSafeFit( ToyStudyMinimiser, ToyStudyFunction, LocalInputFixedSet, PDFsWithDataForToys, ConstraintsForToys, OutputLevel );
				//	Reset Std Output Streams
				if( OutputLevel <= -1 )
				{
					cout.rdbuf(cout_bak);
					cerr.rdbuf(cerr_bak);
					clog.rdbuf(clog_bak);
				}
				cout << "\tFit Finished!\n" <<endl;

				//	If either Fit Failed we want to 'dump the results' and run an extra Fit.
				if( (fit1Result->GetFitStatus() != 3) || (fit2Result->GetFitStatus() != 3) ) toy_failed = true;
			}
			else	toy_failed = false;

			if( !toy_failed || FC_Debug_Flag )
				ResultFormatter::ReviewOutput( fit2Result );

			//	Do we want to store the Data OR run another toy to get a better Fit
			if( toy_failed )
			{
				cerr << "\n\t\tA Single Toy Study Failed... Requesting More Data for another pass.\n" << endl;
				//  Increment counter so we can guarantee we get 'wanted_number_of_toys' which have fitted correctly
				++wanted_number_of_toys;
			}

			//	In my opinion the fact this is so easy and doesn't cause root to throw up is a testiment to the other authors in RapidFit :D
			cout << "\n\n\t\tDeleting Used Data:\n"<<endl;
			for( unsigned int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); ++pdf_num )
			{
				cout << "deleting data for pdf: " << pdf_num << "\n";
				//delete Memory_Data[pdf_num].back();
				Memory_Data[pdf_num].back() = NULL;
			}

			// Only generate data if I'm going to fit to it
			if( dataset_num < (wanted_number_of_toys-1) )
			{
				cout << "\n\n\t\tGenerating Data For Next Toy:\n" <<endl;

				for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); ++pdf_num )
				{
					// See above for more detailed description
					cout << "generating data for pdf: " << pdf_num << "\n";
					PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( ControlParamSet );
					unsigned int wanted_events = EventsPerPDF[pdf_num];
					if( sweight_error[pdf_num] > 0. )  wanted_events *= unsigned(int(new_rand->Gaus()*sweight_error[pdf_num]));
					IDataSet* new_dataset = PDFsWithDataForToys[pdf_num]->GetDataSetConfig()->MakeDataSet( PhaseSpaceForToys[pdf_num], PDFsWithDataForToys[pdf_num]->GetPDF(), int(wanted_events) );
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
		FitResultVector* ThisStudy = new FitResultVector( GlobalFitResult->GetResultParameterSet()->GetAllNames() );
		ThisStudy->AddFitResult( _2DResultForFC->GetFitResult( iFC ), false );
		ThisStudy->AddCPUTimes( _2DResultForFC->GetAllCPUTimes() );
		ThisStudy->AddRealTimes( _2DResultForFC->GetAllRealTimes() );
		//	The Generated Value for the Global and Local fit are best defined as -9999 as a sensible default
		for( unsigned short int num=0; num < GlobalFitResult->GetResultParameterSet()->GetAllNames().size(); ++num )
		{
			string name = GlobalFitResult->GetResultParameterSet()->GetAllNames()[num];
			GlobalFitResult->GetResultParameterSet()->GetResultParameter( name )->ForcePullValue( -9999 );
			GlobalFitResult->GetResultParameterSet()->GetResultParameter( name )->ForceOriginalValue( -9999 );
			ThisStudy->GetFitResult(0)->GetResultParameterSet()->GetResultParameter( name )->ForcePullValue( -9999 );
			ThisStudy->GetFitResult(0)->GetResultParameterSet()->GetResultParameter( name )->ForceOriginalValue( -9999 );
		}
		AllResults.push_back( GlobalResult );
		AllResults.push_back( ThisStudy );
		AllResults.push_back( study1Results );
		AllResults.push_back( study2Results );

		delete new_rand;
	}

	cout << "\n\nNumerical part of FeldMan-Cousins Scan Complete,\n\n Now process the data through the plotting tool :D\n\n";
	return new FitResultVector( AllResults );
}
