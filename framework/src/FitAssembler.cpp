/**
  @class FitAssembler

  The intention is for this class to formalise the process of assembling the components of a fit
  Ideally it will be a set of nested static methods, starting from more and more rudimentary components

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

#include "FitAssembler.h"
#include "ClassLookUp.h"
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
	PhysicsBottle * bottle = new PhysicsBottle(BottleParameters);

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



//====================================================================================================
//Do a likelihood scan
LLscanResult * FitAssembler::DoScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, string scanName, int npoints, double uplim, double lolim )
{
	
	cout << "Performing LL scan for parameter " << scanName << endl ;
	
	// Get a pointer to the physics parameter to be scanned and fix it	
	// CAREFUL:  this must be reset as it was at the end.
	PhysicsParameter * scanParameter = BottleParameters->GetPhysicsParameter(scanName) ;
	double originalValue = scanParameter->GetBlindedValue( ) ;
	string originalType = scanParameter->GetType( ) ;
	scanParameter->SetType( "Fixed" ) ;

	// Need to set up a loop , fixing the scan parameter at each point
	vector<double> scanParameterValues ;
	vector<double> scanLLValues ;
	double deltaScan = (uplim - lolim) / (npoints-1.) ;
		
	for( int si=0; si<npoints; si++) {
			
		// Set scan parameter value
		double scanVal = lolim + si*deltaScan ;
		cout << "Scan value: " << scanVal << endl;
		scanParameter->SetBlindedValue( scanVal ) ;

		cout << "Calling DoFit" << endl;
		// Do a scan point fit
		FitResult * scanStepResult = FitAssembler::DoFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints );
		cout << "Called" << endl;
					
		scanParameterValues.push_back( scanVal ) ;
		if( scanStepResult->GetFitStatus() == 3 ) 
		{
			scanLLValues.push_back( scanStepResult->GetMinimumValue() );
		}
		else{
			scanLLValues.push_back( LLSCAN_FIT_FAILURE_VALUE );
		}			
	
	}

	LLscanResult * result = new LLscanResult( scanName, scanParameterValues, scanLLValues ) ;
	result->print() ; //PELC		
	
	//Reset the parameter as it was
	scanParameter->SetType( originalType ) ;
	scanParameter->SetBlindedValue( originalValue ) ;
	
	return result;
								
}


//===========================================================================================
//Do a double likelihood scan
LLscanResult2D * FitAssembler::DoScan2D( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,
									vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, 
									string scanName, string scanName2, int npoints, double uplim, double lolim,  double uplim2, double lolim2 )
{
	
	cout << "Performing LL scan for parameters " << scanName << "  and  " << scanName2 << endl ;
	
	// Get a pointer to the physics parameter to be scanned and fix it	
	// CAREFUL:  this must be reset as it was at the end.
	PhysicsParameter * scanParameter = BottleParameters->GetPhysicsParameter(scanName) ;
	double originalValue = scanParameter->GetBlindedValue( ) ;
	string originalType = scanParameter->GetType( ) ;
	scanParameter->SetType( "Fixed" ) ;
	
	// Need to set up a loop , fixing the scan parameter at each point
	vector<double> scanParameterValues ;
	vector<LLscanResult*> listOfLLscans ;
	double deltaScan = (uplim-lolim) / (npoints-1.) ;
	
	for( int si=0; si<npoints; si++) {
		
		// Set scan parameter value
		double scanVal = lolim + si*deltaScan ;
		scanParameter->SetBlindedValue( scanVal ) ;
		
		// Do a scan point fit
		LLscanResult * oneScan = FitAssembler::DoScan( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, scanName2, npoints, uplim2, lolim2 );
		
		scanParameterValues.push_back( scanVal ) ;
		listOfLLscans.push_back( oneScan );
		
	}
	
	//To create theresult object, we need to get a list of the scan parameters for scan name 2
	vector<double> scanParameterValues2 ;
	if( npoints > 0 ) scanParameterValues2 = listOfLLscans[0]->GetParameterValues();
	
	LLscanResult2D * result = new LLscanResult2D( scanName, scanParameterValues, scanName2, scanParameterValues2, listOfLLscans ) ;
	
	//Reset the parameter as it was
	scanParameter->SetType( originalType ) ;
	scanParameter->SetBlindedValue( originalValue ) ;
	
	return result;
	
}

