/*!
 * @class FitAssembler
 *
 * The intention is for this class to formalise the process of assembling the components of a fit
 * Ideally it will be a set of nested static methods, starting from more and more rudimentary components
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */


///	RapidFit Headers
#include "FitAssembler.h"
#include "FitResult.h"
#include "ClassLookUp.h"
#include "ScanParam.h"
#include "FitResultVector.h"
#include "ResultFormatter.h"
#include "StringProcessing.h"
///	System Headers
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using namespace::std;


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

	cout << "\nStarting Fit!" << endl;

	SafeMinimise( Minimiser );

	cout << "\nMinimised!\n" << endl;

	FitResult* final_result = Minimiser->GetFitResult();

	return final_result;
}

//Create the minimiser and fit function
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, PhysicsBottle * Bottle )
{

	IMinimiser * minimiser = MinimiserConfig->GetMinimiser( int(Bottle->GetParameterSet()->GetAllNames().size()) );

	FitFunction * theFunction = FunctionConfig->GetFitFunction();
	theFunction->SetPhysicsBottle(Bottle);

	FitResult* result = DoFit( minimiser, theFunction );

	return result;
}

//Create the physics bottle
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet* BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints )
{
	vector<IPDF*> allPDFs;
	for( unsigned int i=0; i< BottleData.size(); ++i )
	{
		allPDFs.push_back( BottleData[i]->GetPDF() );
	}

	ParameterSet* checkedBottleParameters = CheckInputParams( BottleParameters, allPDFs );

	PhysicsBottle * bottle = new PhysicsBottle( checkedBottleParameters );

	//Fill the bottle - data generation occurs in this step
	for ( unsigned int resultIndex = 0; resultIndex < BottleData.size(); ++resultIndex )
	{
		BottleData[resultIndex]->SetPhysicsParameters( checkedBottleParameters );
		IPDF* Requested_PDF = BottleData[resultIndex]->GetPDF();
		IDataSet* Requested_DataSet = BottleData[resultIndex]->GetDataSet();

		bottle->AddResult( Requested_PDF, Requested_DataSet );
	}

	//Add the constraints
	for ( unsigned int constraintIndex = 0; constraintIndex < BottleConstraints.size(); ++constraintIndex )
	{
		bottle->AddConstraint( BottleConstraints[constraintIndex] );
	}

	delete checkedBottleParameters;

	FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );

	return result;
}

//Check that the provided ParameterSet only Contains the Parameters claimed by the PDFs to protect the Minimiser from runtime mistakes
ParameterSet* FitAssembler::CheckInputParams( ParameterSet* givenParams, vector<IPDF*> allPDFs )
{
	vector<string> param_names;
	for( unsigned int i=0; i< allPDFs.size(); ++i )
	{
		vector<string> input_from_PDF = allPDFs[i]->GetPrototypeParameterSet();
		param_names = StringProcessing::CombineUniques( param_names, input_from_PDF );
	}

	if( param_names.size() != givenParams->GetAllNames().size() )
	{
		vector<string> temp = givenParams->GetAllNames();
		for( unsigned int i=0; i< temp.size(); ++i )
		{
			if( StringProcessing::VectorContains( &param_names, &(temp[i]) ) == -1 )
			{
				cout << "Physics Parameter: " << temp[i] << " Provided, but not requested by a PDF, Ignoring" << endl;
			}
		}
	}

	ParameterSet* wantedParameterSet = ( new ParameterSet( param_names ) );
	for( unsigned int i=0; i< param_names.size(); ++i )
	{
		PhysicsParameter* phys_param = givenParams->GetPhysicsParameter( param_names[i] );
		wantedParameterSet->SetPhysicsParameter( param_names[i], new PhysicsParameter(*phys_param) );
	}

	return wantedParameterSet;
}

//Create the physics bottle with pre-made data
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet* BottleParameters,
		vector< IPDF* > AllPDFs, vector< IDataSet* > AllData, vector< ConstraintFunction* > BottleConstraints )
{
	ParameterSet* internalBottleParameters = FitAssembler::CheckInputParams( BottleParameters, AllPDFs );

	if ( AllPDFs.size() == AllData.size() )
	{
		PhysicsBottle * bottle = new PhysicsBottle( internalBottleParameters );

		//Fill the bottle - data already generated
		for ( unsigned int resultIndex = 0; resultIndex < AllData.size(); ++resultIndex )
		{
			AllPDFs[resultIndex]->UpdatePhysicsParameters(internalBottleParameters);
			bottle->AddResult( AllPDFs[resultIndex], AllData[resultIndex] );
		}

		//Add the constraints
		for ( unsigned int constraintIndex = 0; constraintIndex < BottleConstraints.size(); ++constraintIndex )
		{
			bottle->AddConstraint( BottleConstraints[constraintIndex] );
		}

		//bottle->Finalise();
		FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );


		delete internalBottleParameters;

		return result;
	}
	else
	{
		cerr << "Mismatched number of PDFs and DataSets" << endl;
		exit(1);
	}
}

//     This checks the ParameterSet in the Physics Parameters against all of the parameters given in the XML and then adds any missing parameters the user requested
void FitAssembler::CheckParameterSet( FitResult* ReturnableFitResult, ParameterSet* BottleParameters )
{
	vector<string> already_found = ReturnableFitResult->GetResultParameterSet()->GetAllNames();

	for( unsigned int i=0; i< BottleParameters->GetAllNames().size(); ++i )
	{
		int found = StringProcessing::VectorContains( &already_found, &(BottleParameters->GetAllNames()[i]) );

		//	There was something in the ParameterSet not in the FitResult, i.e. an unclaimed object which can't have changed during the fit
		if( found == -1 )
		{
			cout << "ALERT:\t" << "Parameter " << BottleParameters->GetAllNames()[i] << " was not claimed by any PDF in the fit and was NOT passed to the Minimiser!!!" << endl;
			double Value = BottleParameters->GetPhysicsParameter( BottleParameters->GetAllNames()[i] )->GetValue();
			double OriginalValue = BottleParameters->GetPhysicsParameter( BottleParameters->GetAllNames()[i] )->GetValue();
			double Minimum = BottleParameters->GetPhysicsParameter( BottleParameters->GetAllNames()[i] )->GetMinimum();
			double Maximum = BottleParameters->GetPhysicsParameter( BottleParameters->GetAllNames()[i] )->GetMaximum();
			double Error = BottleParameters->GetPhysicsParameter( BottleParameters->GetAllNames()[i] )->GetStepSize();
			string Type = BottleParameters->GetPhysicsParameter( BottleParameters->GetAllNames()[i] )->GetType();
			string Unit = BottleParameters->GetPhysicsParameter( BottleParameters->GetAllNames()[i] )->GetUnit();
			bool added = ReturnableFitResult->GetResultParameterSet()->ForceNewResultParameter( BottleParameters->GetAllNames()[i],  Value, OriginalValue, Error, Minimum, Maximum, Type, Unit );
			if( !added )
			{
				cerr << "Error finalizing FitResultVector Object" << endl << endl;
				exit(-984);
			}
		}
	}
}

//  Perform a safer fit which should always to return something which you can use :D
FitResult * FitAssembler::DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet* BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, int OutputLevel )
{
	FitResult* ReturnableFitResult=NULL;

	ParameterSet* internalBottleParameters = new ParameterSet( *BottleParameters );

	//	Decide what 'Strategy' to use By default just perform the fit
	if( FunctionConfig->GetStrategy() == "Petes" )
	{
		ReturnableFitResult = Petes_DoSafeFit( MinimiserConfig, FunctionConfig, internalBottleParameters, BottleData, BottleConstraints, OutputLevel );
	}
	else if( FunctionConfig->GetStrategy() == "PetesGamma" )
	{
		ReturnableFitResult = PetesGamma_DoSafeFit( MinimiserConfig, FunctionConfig, internalBottleParameters, BottleData, BottleConstraints, OutputLevel );
	}
	else if( FunctionConfig->GetStrategy() == "Robs" )
	{
		ReturnableFitResult = Robs_DoSafeFit( MinimiserConfig, FunctionConfig, internalBottleParameters, BottleData, BottleConstraints, OutputLevel );
	}
	else  //	Default!
	{
		ReturnableFitResult = DoSingleSafeFit( MinimiserConfig, FunctionConfig, internalBottleParameters, BottleData, BottleConstraints, OutputLevel );
	}

	/*!	Check that the ParameterSet in the output contains the provided ParameterSet regardless of what was passed to the Minimiser	*/
	CheckParameterSet( ReturnableFitResult, internalBottleParameters );

	return ReturnableFitResult;
}

FitResult * FitAssembler::Robs_DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, const ParameterSet* BottleParameters,
		const vector< PDFWithData* > BottleData, const vector< ConstraintFunction* > BottleConstraints, const int OutputLevel )
{
	cout << "Starting Fit1:" << endl;
	double deltaPara = BottleParameters->GetPhysicsParameter( string("delta_para") )->GetBlindedValue();
	double deltaPerp = BottleParameters->GetPhysicsParameter( string("delta_perp") )->GetBlindedValue();
	double deltaGamma = BottleParameters->GetPhysicsParameter( string("deltaGamma") )->GetBlindedValue();
	double phi_s = BottleParameters->GetPhysicsParameter( string("Phi_s") )->GetBlindedValue();

	ParameterSet* set1 = new ParameterSet( *BottleParameters );

	// Normal fit
	FitResult* res0 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, set1, BottleData, BottleConstraints, OutputLevel ) ;
	bool good_result_0 = res0->GetFitStatus() == 3;
	cout << "Finished Fit1." << endl;
	if( !good_result_0  ) cout << "Fit-1 failed" << endl;

	ParameterSet* set2 = new ParameterSet( *BottleParameters );

	// Conjugate fit
	//double deltaS = BottleParameters->GetPhysicsParameter( string("delta_s") )->GetBlindedValue() ;   
	if( set2->GetPhysicsParameter( string("delta_para") )->GetType() != "Fixed" ) set2->GetPhysicsParameter( string("delta_para") )->SetBlindedValue( -deltaPara ) ;
	if( set2->GetPhysicsParameter( string("delta_perp") )->GetType() != "Fixed" ) set2->GetPhysicsParameter( string("delta_perp") )->SetBlindedValue( 3.14159-deltaPerp ) ;
	//BottleParameters->GetPhysicsParameter( string("delta_s") )->SetBlindedValue( -deltaS );
	cout << endl << "Starting Fit2:" << endl;
	FitResult* res1 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, set2, BottleData, BottleConstraints, OutputLevel ) ;
	cout << "Finished Fit2." << endl;
	bool good_result_1 = res1->GetFitStatus() == 3;

	ParameterSet* set3 = new ParameterSet( *BottleParameters );

	if( set3->GetPhysicsParameter( string("delta_para") )->GetType() != "Fixed" ) set3->GetPhysicsParameter( string("delta_para") )->SetBlindedValue( deltaPara ) ;
	if( set3->GetPhysicsParameter( string("delta_perp") )->GetType() != "Fixed" ) set3->GetPhysicsParameter( string("delta_perp") )->SetBlindedValue( deltaPerp ) ;
	if( set3->GetPhysicsParameter( string("deltaGamma") )->GetType() != "Fixed" ) set3->GetPhysicsParameter( string("deltaGamma") )->SetBlindedValue( -deltaGamma ) ;
	if( set3->GetPhysicsParameter( string("Phi_s") )->GetType() != "Fixed" ) set3->GetPhysicsParameter( string("Phi_s") )->SetBlindedValue( 3.14159-phi_s ) ;
	cout << endl << "Starting Fit3:" << endl;
	FitResult* res2 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, set3, BottleData, BottleConstraints, OutputLevel ) ;
	cout << "Finished Fit3." << endl;
	bool good_result_2 = res2->GetFitStatus() == 3;

	ParameterSet* set4 = new ParameterSet( *BottleParameters );

	if( set4->GetPhysicsParameter( string("delta_para") )->GetType() != "Fixed" ) set4->GetPhysicsParameter( string("delta_para") )->SetBlindedValue( -deltaPara ) ;
	if( set4->GetPhysicsParameter( string("delta_perp") )->GetType() != "Fixed" ) set4->GetPhysicsParameter( string("delta_perp") )->SetBlindedValue( 3.14159-deltaPerp ); 
	if( set4->GetPhysicsParameter( string("deltaGamma") )->GetType() != "Fixed" ) set4->GetPhysicsParameter( string("deltaGamma") )->SetBlindedValue( -deltaGamma ) ;
	if( set4->GetPhysicsParameter( string("Phi_s") )->GetType() != "Fixed" ) set4->GetPhysicsParameter( string("Phi_s") )->SetBlindedValue( 3.14159-phi_s ) ;
	cout << endl << "Starting Fit4:" << endl;
	FitResult* res3 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, set4, BottleData, BottleConstraints, OutputLevel ) ;
	cout << "Finished Fit4." << endl;
	bool good_result_3 = res3->GetFitStatus() == 3;


	// Apply strategy 
	FitResult * res = res0;
	bool swapped = false ;

	vector<FitResult*> results;
	if( good_result_0 ) results.push_back( res0 ); else results.push_back( NULL );
	if( good_result_1 ) results.push_back( res1 ); else results.push_back( NULL );
	if( good_result_2 ) results.push_back( res2 ); else results.push_back( NULL );
	if( good_result_3 ) results.push_back( res3 ); else results.push_back( NULL );

	for( unsigned int i=0; i< results.size(); ++i )
	{
		if( results[i] != NULL )
		{
			if( results[i]->GetMinimumValue() < res->GetMinimumValue() )
			{
				res = results[i];
				cout << "Replacing Output With FitResult " << i << endl;
			}
		}
	}

	swapped = (res != res0);

	cout << endl << "Fit Finished!" << endl;

	//MinimiserConfig->SetOutputLevel( 0 );
	if( swapped  ) cout << "Swapped to some conjugate fit result " << endl;

	delete set1; delete set2; delete set3; delete set4;

	return res;
}

FitResult * FitAssembler::DoSingleSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet* BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, int OutputLevel )
{
	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;

	//      If the user wanted silence we point the Std Output Streams to the oblivion of NULL
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		cout.rdbuf(0);
		cerr.rdbuf(0);
		clog.rdbuf(0);
	}

	FitResult* ReturnableFitResult=NULL;
	MinimiserConfig->SetOutputLevel( OutputLevel );
	vector<string> other_params = BottleParameters->GetAllFloatNames();      //      This better at least contain all in prototypeparamset!!!
	vector<double> truth;
	for( unsigned short int j=0; j < other_params.size(); ++j )
	{
		truth.push_back( BottleParameters->GetPhysicsParameter( other_params[j] )->GetValue() );
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

	bool bad_fit=false;

	//	Have to protect against null objects being returned from the Wrappers
	if( ReturnableFitResult != NULL )
	{
		if( ReturnableFitResult->GetFitStatus() != 3 )
		{
			bad_fit = true;
		}
	} else {
		bad_fit = true;
	}

	if( bad_fit )
	{
		cerr << "\n\n\t\tFit Did NOT Converge Correctly, CHECK YOUR RESULTS!\n\n";
		for( unsigned short int j=0; j < other_params.size(); ++j )
		{
			BottleParameters->GetPhysicsParameter( other_params[j] )->SetValue( truth[j] );
		}
		int status = -1;
		vector<string> NewNamesList = BottleParameters->GetAllNames();
		ResultParameterSet* DummyFitResults = new ResultParameterSet( NewNamesList );
		for( unsigned int j=0; j< NewNamesList.size(); ++j )
		{
			double val = -999.;
			double origval = -999.;
			double err = -999.;
			//double gen = -999.;
			double min = -999.;
			double max = -999.;
			string type = BottleParameters->GetPhysicsParameter( NewNamesList[j] )->GetType();
			string unit = BottleParameters->GetPhysicsParameter( NewNamesList[j] )->GetUnit();
			DummyFitResults->SetResultParameter( NewNamesList[j], val, origval, err, min, max, type, unit );
		}

		PhysicsBottle* Bad_Bottle = new PhysicsBottle( BottleParameters );
		ReturnableFitResult = new FitResult( LLSCAN_FIT_FAILURE_VALUE, DummyFitResults, status, Bad_Bottle );
	}

	return ReturnableFitResult;
}


FitResult * FitAssembler::Petes_DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet* BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, int OutputLevel )
{
	cout << endl << "******* Result of Petes Double fit strategy  (NOT flipping delta_s as well) *********" << endl ;
	cout << "Starting Fit1:" << endl;
	// Normal fit
	FitResult* res0 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, OutputLevel ) ;
	double LLmin0 = res0->GetMinimumValue() ;
	bool good_result_0 = res0->GetFitStatus() == 3;
	cout << "Finished Fit1." << endl;
	if( !good_result_0  ) cout << "Fit-1 failed" << endl;

	// Conjugate fit
	double deltaPara = BottleParameters->GetPhysicsParameter( string("delta_para") )->GetBlindedValue();
	double deltaPerp = BottleParameters->GetPhysicsParameter( string("delta_perp") )->GetBlindedValue();
	//double deltaS = BottleParameters->GetPhysicsParameter( string("delta_s") )->GetBlindedValue();
	BottleParameters->GetPhysicsParameter( string("delta_para") )->SetBlindedValue( -deltaPara );
	BottleParameters->GetPhysicsParameter( string("delta_perp") )->SetBlindedValue( 3.14159-deltaPerp );
	//BottleParameters->GetPhysicsParameter( string("delta_s") )->SetBlindedValue( -deltaS );
	cout << endl << "Starting Fit2:" << endl;
	FitResult* res1 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, OutputLevel );
	double LLmin1 = res1->GetMinimumValue();
	bool good_result_1 = res1->GetFitStatus() == 3;

	// Apply strategy 
	FitResult * res=NULL;
	bool swapped = false;
	if( !good_result_0 && !good_result_1 )
	{
		res = res0;     //both fail so it doesnt matter 
	}
	else if( good_result_0 && !good_result_1 )
	{  
		res = res0;
	}
	else if( !good_result_0 && good_result_1 )
	{
		res = res1;
		swapped = true;
	}
	else if( LLmin0 < LLmin1+0.001 )
	{
		res = res0;
	}
	else
	{
		res = res1;
		swapped = true;
	}

	//MinimiserConfig->SetOutputLevel( 0 );
	cout << setprecision(9);
	if( !good_result_1  ) cout << "Fit-2 failed" << endl;
	cout << "Finished Fit2." << endl;
	cout << "The 2 LLs were "  << LLmin0 << "   <=>   "  << LLmin1 << endl ;
	if( swapped  ) cout << "Swapped to conjugate fit result " << endl;
	cout << "******* Result of Petes Double fit strategy*********" << endl ;
	//MinimiserConfig->SetOutputLevel( OutputLevel );

	// Set input parameters back to what they were 
	BottleParameters->GetPhysicsParameter( string("delta_para") )->SetBlindedValue( deltaPara );
	BottleParameters->GetPhysicsParameter( string("delta_perp") )->SetBlindedValue( deltaPerp );

	return res;
}

FitResult * FitAssembler::PetesGamma_DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet* BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, int OutputLevel )
{
	cout << endl << "******* Result of Petes Double fit strategy for flipping deltaGamma *********" << endl ;
	cout << "Starting Fit1:" << endl;
	// Normal fit
	FitResult* res0 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, OutputLevel ) ;
	double LLmin0 = res0->GetMinimumValue() ;
	bool good_result_0 = res0->GetFitStatus() == 3;
	cout << "Finished Fit1." << endl;
	if( !good_result_0  ) cout << "Fit-1 failed" << endl;
	
	// Conjugate fit
	double deltaGamma = BottleParameters->GetPhysicsParameter( string("deltaGamma") )->GetBlindedValue() ;
	BottleParameters->GetPhysicsParameter( string("deltaGamma") )->SetBlindedValue( -deltaGamma ) ;
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
	cout << "******* Result of Petes flipping deltaGamma strategy*********" << endl ;
	//MinimiserConfig->SetOutputLevel( OutputLevel );
	
	// Set input parameters back to what they were 
	BottleParameters->GetPhysicsParameter( string("deltaGamma") )->SetBlindedValue( deltaGamma ) ;  
	
	return res ;
	
}

