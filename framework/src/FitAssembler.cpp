/**
  @class FitAssembler

  The intention is for this class to formalise the process of assembling the components of a fit
  Ideally it will be a set of nested static methods, starting from more and more rudimentary components

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

#include "FitAssembler.h"
#include "FitResult.h"
#include "ClassLookUp.h"
#include "ScanParam.h"
#include "ToyStudyResult.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

//The final stage - do the minimisation
FitResult * FitAssembler::DoFit( IMinimiser * Minimiser, FitFunction * TheFunction )
{
	Minimiser->Minimise(TheFunction);
	return Minimiser->GetFitResult();
}

//Create the minimiser and fit function
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, PhysicsBottle * Bottle )
{
	IMinimiser * minimiser = MinimiserConfig->GetMinimiser( Bottle->GetParameterSet()->GetAllNames().size() );
	FitFunction * theFunction = FunctionConfig->GetFitFunction();
	theFunction->SetPhysicsBottle(Bottle);

	FitResult * result = DoFit( minimiser, theFunction );

	delete theFunction;
	delete minimiser;
	return result;
}

//Create the physics bottle
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints )
{
	PhysicsBottle * bottle = new PhysicsBottle( BottleParameters );

	//Fill the bottle - data generation occurs in this step
	for ( int resultIndex = 0; resultIndex < BottleData.size(); resultIndex++ )
	{
		BottleData[resultIndex]->SetPhysicsParameters(BottleParameters);
		bottle->AddResult( BottleData[resultIndex]->GetPDF(), BottleData[resultIndex]->GetDataSet() );
	}

	//Add the constraints
	for ( int constraintIndex = 0; constraintIndex < BottleConstraints.size(); constraintIndex++ )
	{
		bottle->AddConstraint( BottleConstraints[constraintIndex] );
	}

	bottle->Finalise();
	FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );

	delete bottle;
	return result;
}

//Create the physics bottle with pre-made data
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,
		vector< IPDF* > AllPDFs, vector< IDataSet* > AllData, vector< ConstraintFunction* > BottleConstraints )
{
	if ( AllPDFs.size() == AllData.size() )
	{
		PhysicsBottle * bottle = new PhysicsBottle(BottleParameters);

		//Fill the bottle - data already generated
		for ( int resultIndex = 0; resultIndex < AllData.size(); resultIndex++ )
		{
			AllPDFs[resultIndex]->SetPhysicsParameters(BottleParameters);
			bottle->AddResult( AllPDFs[resultIndex], AllData[resultIndex] );
		}

		//Add the constraints
		for ( int constraintIndex = 0; constraintIndex < BottleConstraints.size(); constraintIndex++ )
		{
			bottle->AddConstraint( BottleConstraints[constraintIndex] );
		}  

		bottle->Finalise();
		FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );

		delete bottle;
		return result;
	}
	else
	{
		cerr << "Mismatched number of PDFs and DataSets" << endl;
		exit(1);
	}
}


//  Interface for internal calls
void FitAssembler::DoScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters, vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, ScanParam* Wanted_Param, ToyStudyResult* output_interface )
{

	double uplim = Wanted_Param->GetMax();
	double lolim = Wanted_Param->GetMin();
	double npoints = Wanted_Param->GetPoints();
	string scanName = Wanted_Param->GetName();

	cout << "Performing Scan for the parameter " << scanName << endl ;
	
	// Get a pointer to the physics parameter to be scanned and fix it	
	// CAREFUL:  this must be reset as it was at the end.
	PhysicsParameter * scanParameter = BottleParameters->GetPhysicsParameter(scanName);
	double originalValue = scanParameter->GetBlindedValue( ) ;
	string originalType = scanParameter->GetType( ) ;
	scanParameter->SetType( "Fixed" ) ;

	// Need to set up a loop , fixing the scan parameter at each point
	double deltaScan = (uplim - lolim) / (npoints-1.) ;

	for( int si=0; si<npoints; si++) {

		// Set scan parameter value
		double scanVal = lolim + si*deltaScan;
		//scanParameter->SetBlindedValue( scanVal ) ;

		output_interface->StartStopwatch();

		// Do a scan point fit
		FitResult * scanStepResult;


		double realscanVal=scanVal;
		//  fraction of the full step to dance around the fixed value to find a working one
		double retry_frac=10.;
		//  Did any fit converge?
		bool fit_try_failed=false;


		//  Here there be Monsters!!!...
		//
		//  From here I will attempt to perform a fit and catch any fatal errors which cause the fit to abort
		//  Should the fit abort I will shift the value I am fixing by some amount and try again
		//  I will repeat this 4 times until I have had 5 fit failiures then I will deam this to be a lost cause and the users problem
		//  This can be coded up into an even more generic for try catch but I feal the algorithm would be painful to code and even more painful to read so I think this has the best performance/transparency

		for( short int i=0; i<=4; i++ )
		{
			fit_try_failed = false;
			//  Give a Fit 3 chances to work, fit at origin, then fit a small step either side of a not working result to see if that helps
			if( i == 1 ) realscanVal = scanVal+deltaScan/retry_frac; 
			if( i == 2 ) realscanVal = scanVal-deltaScan/retry_frac;
			if( i == 3 ) realscanVal = scanVal+2.*(deltaScan/retry_frac);
			if( i == 4 ) realscanVal = scanVal-2.*(deltaScan/retry_frac);

			scanParameter->SetBlindedValue( realscanVal ) ;

			try{
				// Try a Fit, it it converges, continue to elsewhere in the program
				scanStepResult = FitAssembler::DoFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints );
				continue;
			}
			//  If it didn't converge tell the user why!
			catch( int e){
				if ( e == 10 ){
				cerr << "\nCaught exception : fit failed for these parameters..." << endl;  }
				else if ( e == 13 ){
				cerr << "\nIntegration Error: Fit Failed..." << endl;  }
				fit_try_failed = true;
			} catch (...) {
                        	cerr << "\n\n\n\t\t\tCaught Unknown Exception, THIS IS SERIOUS!!!\n\n\n" << endl;
				fit_try_failed = true;
			}
			//  Exit at 4 as we no longer are retrying anything
			if( i == 4 ) continue;

			//  Let the user know that I'm trying another fit close to the value they originally asked for.
			cerr << "Retrying the fit again with ";
			if( (i == 0) || (i == 1) ) cerr << "1";
			if( (i == 2) || (i == 3) ) cerr << "2";
			cerr << " small steps ";
			if( (i == 0) || (i == 2) ) cerr <<"forward";
			if( (i == 1) || (i == 3) ) cerr << "backward";
			cerr << " from the failing fixed value" << endl;
		}

		//  If the fit failed 5 times I will simply return a dummy fit result full of zerod objects. It is up to the user to watch for and remove these
		if( fit_try_failed ){
			cerr << "Nothing more I'm willing to do, considering a Fit here a lost cause..." <<endl;
			realscanVal = scanVal;
			scanParameter->SetBlindedValue( realscanVal ) ;
			int status = -1;
			vector<string> NewNamesList = BottleParameters->GetAllNames();
			ResultParameterSet* DummyFitResults = new ResultParameterSet( NewNamesList );
			scanStepResult = new FitResult( LLSCAN_FIT_FAILURE_VALUE, DummyFitResults, status, BottleParameters );
		}

		//  THIS IS ALWAYS TRUE BY DEFINITION OF THE SCAN
		string name = Wanted_Param->GetName();
		string type = BottleParameters->GetPhysicsParameter( name )->GetType();
		string unit = BottleParameters->GetPhysicsParameter( name )->GetUnit();
		scanStepResult->GetResultParameterSet()->SetResultParameter( name, realscanVal, realscanVal, 0., realscanVal, realscanVal, type, unit );

		vector<string> Fixed_List = BottleParameters->GetAllFixedNames();
		vector<string> Fit_List = scanStepResult->GetResultParameterSet()->GetAllNames();
		for( short int i=0; i < Fixed_List.size() ; i++ )
		{
			bool found=false;
			for( short int j=0; j < Fit_List.size(); j++ )
			{
				if( Fit_List[j] == Fixed_List[i] )
				{
					found = true;
				}
			}
			if( !found )
			{
				string fixed_type = BottleParameters->GetPhysicsParameter( Fixed_List[i] )->GetType();
				string fixed_unit = BottleParameters->GetPhysicsParameter( Fixed_List[i] )->GetUnit();
				double fixed_value = BottleParameters->GetPhysicsParameter( Fixed_List[i] )->GetValue();
				scanStepResult->GetResultParameterSet()->ForceNewResultParameter( Fixed_List[i], fixed_value, fixed_value, 0, fixed_value, fixed_value, fixed_type, fixed_unit );
			}
		}


		output_interface->AddFitResult( scanStepResult );
	}
	
	//Reset the parameter as it was
	scanParameter->SetType( originalType ) ;
	scanParameter->SetBlindedValue( originalValue ) ;

}

//  Interface for internal calls
void FitAssembler::DoScan2D( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters, vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, pair<ScanParam*, ScanParam*> Param_Set, vector<ToyStudyResult*>* output_interface )
{
//	vector<string> namez = BottleParameters->GetAllNames();
	vector<string> result_names = BottleParameters->GetAllNames();

	double uplim = Param_Set.first->GetMax();
	double lolim = Param_Set.first->GetMin();
	double npoints = Param_Set.first->GetPoints();

	string scanName = Param_Set.first->GetName();
	string scanName2 = Param_Set.second->GetName();


	// Get a pointer to the physics parameter to be scanned and fix it
	// CAREFUL:  this must be reset as it was at the end.
	PhysicsParameter * scanParameter = BottleParameters->GetPhysicsParameter(scanName);
	double originalValue = scanParameter->GetBlindedValue( );
	string originalType = scanParameter->GetType( );
	scanParameter->SetType( "Fixed" );
	
	// Need to set up a loop , fixing the scan parameter at each point

	double deltaScan = (uplim-lolim) / (npoints-1.) ;
	
	for( int si=0; si < npoints; si++) {
		ToyStudyResult* Returnable_Result = new ToyStudyResult( result_names );
		
		// Set scan parameter value
		double scanVal = lolim + si*deltaScan ;
		scanParameter->SetTrueValue( scanVal ) ;

		// Do a scan point fit
		FitAssembler::DoScan( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, Param_Set.second, Returnable_Result );

		//  THIS IS ALWAYS TRUE BY DEFINITION OF THE SCAN
                string name = Param_Set.first->GetName();
                string type = BottleParameters->GetPhysicsParameter( name )->GetType();
                string unit = BottleParameters->GetPhysicsParameter( name )->GetUnit();

		for( short int i=0; i < Returnable_Result->NumberResults(); i++ )
		{
			Returnable_Result->GetFitResult( i )->GetResultParameterSet()->SetResultParameter( name, scanVal, scanVal, 0.0, scanVal, scanVal, type, unit );
		}

		output_interface->push_back( Returnable_Result );
	}

	//Reset the parameter as it was
	scanParameter->SetType( originalType ) ;
	scanParameter->SetBlindedValue( originalValue ) ;

}

// Interface for external calls
vector<ToyStudyResult*> FitAssembler::ContourScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters, vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, OutputConfiguration* OutputConfig, string scanName, string scanName2 )
{
	vector<ToyStudyResult*>* Returnable_Result = new vector<ToyStudyResult*>;

	pair< ScanParam*, ScanParam* > Param_Set = OutputConfig->Get2DScanParams( scanName, scanName2 );

	DoScan2D( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, Param_Set, Returnable_Result );

	return *Returnable_Result;
}

//  Interface for external calls
ToyStudyResult* FitAssembler::SingleScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters, vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, OutputConfiguration* OutputConfig, string scanName )
{
	ToyStudyResult* Returnable_Result = new ToyStudyResult( BottleParameters->GetAllNames() );

	ScanParam* local_param = OutputConfig->GetScanParam( scanName );

	DoScan( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, local_param, Returnable_Result );

	return Returnable_Result;
}

