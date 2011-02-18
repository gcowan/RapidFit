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
#include "LLscanResult2D.h"
#include "LLscanResult.h"
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


vector<LLscanResult*> FitAssembler::DoCVScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, OutputConfiguration* OutputConfig, string scanName )
{
	ScanParam* temp_param = OutputConfig->GetScanParam( scanName, "CVscan" );

	vector<LLscanResult*> Returnable_Scan_Results;


	PhysicsParameter * ScanParameter = BottleParameters->GetPhysicsParameter(scanName);
	vector<string> all_parameters = BottleParameters->GetAllNames();
	vector<string> parameters_of_interest;
	vector<TString> Plotting_string;
	for( short int i=0; i < all_parameters.size(); i++ )
	{
		if( ( ScanParameter->GetType() != "Fixed" ) && ( all_parameters[i] != scanName ) )
		{
			parameters_of_interest.push_back( all_parameters[i] );
			TString local_Plotting_String( all_parameters[i] );
			local_Plotting_String.Append(":");
			local_Plotting_String.Append( scanName );
		}
	}

	vector<FitResult*> ScanResults = FitAssembler::DoScan( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, temp_param );

	//  Create containers for persistant results here
	vector<double> axis_x;

	for( short int i=0; i < parameters_of_interest.size(); i++ )
	{
		//  Create containers for non-persistant results here
		vector<double> axis_y;

		for( short int j=0; j < ScanResults.size(); j++ )
		{
			//  Only need to fill X axis once as this is the controlled parameter i.e. same in all plots
			if( axis_x.size() != ScanResults.size() ) {
			       axis_x.push_back( ScanResults[j]->GetResultParameterSet()->GetResultParameter( scanName )->GetValue() ); }//  X Axis
			
			//  Fill the Y axis parameters for all results
			axis_y.push_back( ScanResults[j]->GetResultParameterSet()->GetResultParameter( parameters_of_interest[i] )->GetValue() );//  Y Axis
		}

		LLscanResult * local_result = new LLscanResult( scanName, axis_x, axis_y );

		Returnable_Scan_Results.push_back( local_result );
	}

	return Returnable_Scan_Results;
}

//====================================================================================================
//Do a likelihood scan
LLscanResult * FitAssembler::DoLLScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, OutputConfiguration* OutputConfig, string scanName )
{
	ScanParam* temp_param = OutputConfig->GetScanParam( scanName, "LLscan" );
	return FitAssembler::DoLLScan( MinimiserConfig,  FunctionConfig, BottleParameters,  BottleData,  BottleConstraints, temp_param );
}


LLscanResult * FitAssembler::DoLLScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, ScanParam* Wanted_Param)
{
  
	string scanName = Wanted_Param->GetName();

	cout << "Constructing LLscan for parameter " << scanName << endl ;
	
	// Need to set up a loop , fixing the scan parameter at each point
	vector<double> scanParameterValues ;
	vector<double> scanLLValues ;

	vector<FitResult*> Fitting_Results = FitAssembler::DoScan( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, Wanted_Param );

	for( int si=0; si<Fitting_Results.size(); si++) {

		scanParameterValues.push_back( Fitting_Results[si]->GetResultParameterSet()->GetResultParameter(scanName)->GetValue() );

		if( Fitting_Results[si]->GetFitStatus() == 3 ) 
		{
			scanLLValues.push_back( Fitting_Results[si]->GetMinimumValue() );
		}
		else{
			scanLLValues.push_back( LLSCAN_FIT_FAILURE_VALUE );
		}
	}

	LLscanResult * result = new LLscanResult( scanName, scanParameterValues, scanLLValues ) ;
	result->print() ; //PELC
	
	return result;
}


//====================================================================================================
//Do a likelihood scan
vector<FitResult*> FitAssembler::DoScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters, vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, ScanParam* Wanted_Param )
{

	double uplim = Wanted_Param->GetMax();
	double lolim = Wanted_Param->GetMin();
	double npoints = Wanted_Param->GetPoints();
	string scanName = Wanted_Param->GetName();

	vector<FitResult*> Returnable_Fit_Results;

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
		double scanVal = lolim + si*deltaScan ;
		scanParameter->SetBlindedValue( scanVal ) ;

		// Do a scan point fit
		FitResult * scanStepResult;

		try{
			scanStepResult= FitAssembler::DoFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints );
		}
		catch( int e){
			cerr << "Caught exception : fit failed for these parameters - continuing to next scan value" << endl;
			int status = -1;
			vector<string> NewNamesList = BottleParameters->GetAllNames();
			ResultParameterSet* DummyFitResults = new ResultParameterSet( NewNamesList );
			scanStepResult = new FitResult( LLSCAN_FIT_FAILURE_VALUE, DummyFitResults, status, BottleParameters );
		}

		//  THIS IS ALWAYS TRUE BY DEFINITION OF THE SCAN
		string name = Wanted_Param->GetName();
		string type = BottleParameters->GetPhysicsParameter( name )->GetType();
		string unit = BottleParameters->GetPhysicsParameter( name )->GetUnit();
		scanStepResult->GetResultParameterSet()->SetResultParameter( name, scanVal, scanVal, 0., scanVal, scanVal, type, unit );
		Returnable_Fit_Results.push_back( scanStepResult );

	}
	
	//Reset the parameter as it was
	scanParameter->SetType( originalType ) ;
	scanParameter->SetBlindedValue( originalValue ) ;

	return Returnable_Fit_Results;

}

LLscanResult2D * FitAssembler::DoLLScan2D( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters, vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, OutputConfiguration* OutputConfig, string scanName, string scanName2 )
{
  
	vector<vector<FitResult*> > Fit_Results = DoScan2D( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, OutputConfig, scanName, scanName2, "LLscan", "LLscan" );
	
	cout << "Scan Over" << endl;
	vector<LLscanResult* > LLScanResults;
	vector<double> scanParameterValues;
	vector<double> scanParameterValues3;

	for( int si=0; si < Fit_Results.size(); si++) {

		vector<double> scanLLValues;
		vector<double> scanParameterValues2;

		for( int si2=0; si2 < Fit_Results[si].size(); si2++) {

			scanParameterValues2.push_back( Fit_Results[si][si2]->GetResultParameterSet()->GetResultParameter(scanName2)->GetValue() );

			if( Fit_Results[si][si2]->GetFitStatus() == 3 ) 
			{
				scanLLValues.push_back( Fit_Results[si][si2]->GetMinimumValue() );
			}
			else{
				scanLLValues.push_back( LLSCAN_FIT_FAILURE_VALUE );
			}
			if( si == 0 ) scanParameterValues3.push_back( scanParameterValues2.back() );
		}

		LLscanResult * _1D_temp_result = new LLscanResult( scanName2, scanParameterValues2, scanLLValues ) ;
		LLScanResults.push_back( _1D_temp_result );

		scanParameterValues.push_back( Fit_Results[si][0]->GetResultParameterSet()->GetResultParameter(scanName)->GetValue() );
	}

	LLscanResult2D * result = new LLscanResult2D( scanName, scanParameterValues, scanName2, scanParameterValues3, LLScanResults ) ;

	return result;
}

//===========================================================================================
//Do a double likelihood scan
vector<vector<FitResult* > > FitAssembler::DoScan2D( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,
									vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, 
									OutputConfiguration* OutputConfig, string scanName, string scanName2, string scanType, string scanType2 )
{

	pair< ScanParam*, ScanParam* > Param_Set = OutputConfig->Get2DScanParams( scanName, scanName2, scanType, scanType2 );

	double uplim = Param_Set.first->GetMax();
	double lolim = Param_Set.first->GetMin();
	double npoints = Param_Set.first->GetPoints();

	cout << "Performing LL scan for parameters " << scanName << "  and  " << scanName2 << endl ;
	
	// Get a pointer to the physics parameter to be scanned and fix it	
	// CAREFUL:  this must be reset as it was at the end.
	PhysicsParameter * scanParameter = BottleParameters->GetPhysicsParameter(scanName) ;
	double originalValue = scanParameter->GetBlindedValue( ) ;
	string originalType = scanParameter->GetType( ) ;
	scanParameter->SetType( "Fixed" ) ;
	
	// Need to set up a loop , fixing the scan parameter at each point
	vector<vector<FitResult*> > Scan_Results;

	double deltaScan = (uplim-lolim) / (npoints-1.) ;
	
	for( int si=0; si < npoints; si++) {
		
		// Set scan parameter value
		double scanVal = lolim + si*deltaScan ;
		scanParameter->SetBlindedValue( scanVal ) ;

		// Do a scan point fit
		vector<FitResult*> oneScan = FitAssembler::DoScan( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, Param_Set.second );

		//  THIS IS ALWAYS TRUE BY DEFINITION OF THE SCAN
                string name = Param_Set.first->GetName();
                string type = BottleParameters->GetPhysicsParameter( name )->GetType();
                string unit = BottleParameters->GetPhysicsParameter( name )->GetUnit();

                for( short int i=0; i < oneScan.size(); i++ )
		{
			oneScan[i]->GetResultParameterSet()->SetResultParameter( name, scanVal, scanVal, 0.0, scanVal, scanVal, type, unit );
                }

		Scan_Results.push_back( oneScan );
	}

	//Reset the parameter as it was
	scanParameter->SetType( originalType ) ;
	scanParameter->SetBlindedValue( originalValue ) ;

	return Scan_Results;
}

vector<FitResult* > FitAssembler::Linearize( vector<vector<FitResult*> > Input_Results )
{
	vector<FitResult*> Output_Results;
	for( short int i=0; i< Input_Results.size(); i++)
	{
		for( short int j=0; j< Input_Results[i].size(); j++ )
		{
			Output_Results.push_back( Input_Results[i][j] );
		}
	}
	return Output_Results;
}
