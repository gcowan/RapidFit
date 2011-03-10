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
	IMinimiser * minimiser = MinimiserConfig->GetMinimiser( int(Bottle->GetParameterSet()->GetAllNames().size()) );
	FitFunction * theFunction = FunctionConfig->GetFitFunction();
	theFunction->SetPhysicsBottle(Bottle);

	FitResult * result = DoFit( minimiser, theFunction );

	delete theFunction;
	delete minimiser;
	return result;
}

//Create the physics bottle
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints )
{
	PhysicsBottle * bottle = new PhysicsBottle( BottleParameters );

	//Fill the bottle - data generation occurs in this step
	for ( unsigned int resultIndex = 0; resultIndex < BottleData.size(); resultIndex++ )
	{
		BottleData[resultIndex]->SetPhysicsParameters(BottleParameters);
		IPDF* Requested_PDF = BottleData[resultIndex]->GetPDF();
		IDataSet* Requested_DataSet = BottleData[resultIndex]->GetDataSet();
		bottle->AddResult( Requested_PDF, Requested_DataSet );
	}

	//Add the constraints
	for ( unsigned int constraintIndex = 0; constraintIndex < BottleConstraints.size(); constraintIndex++ )
	{
		bottle->AddConstraint( BottleConstraints[constraintIndex] );
	}

	bottle->Finalise();
	FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );

	delete bottle;
	return result;
}

//Create the physics bottle with pre-made data
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters, vector< IPDF* > AllPDFs, vector< IDataSet* > AllData, vector< ConstraintFunction* > BottleConstraints )
{
	if ( AllPDFs.size() == AllData.size() )
	{
		PhysicsBottle * bottle = new PhysicsBottle(BottleParameters);

		//Fill the bottle - data already generated
		for ( unsigned int resultIndex = 0; resultIndex < AllData.size(); resultIndex++ )
		{
			AllPDFs[resultIndex]->SetPhysicsParameters(BottleParameters);
			bottle->AddResult( AllPDFs[resultIndex], AllData[resultIndex] );
		}

		//Add the constraints
		for ( unsigned int constraintIndex = 0; constraintIndex < BottleConstraints.size(); constraintIndex++ )
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

//void FitAssembler::ShakeBottle( ParameterSet* BottleParameters, vector< PDFWithData* > BottleData, unsigned int some_number )
//{
//	// To be written in a 'safe' way
//}
//  Perform a safer fit which is gauranteed to return something which you can use :D
FitResult * FitAssembler::DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints )
{
	vector<string> other_params = BottleParameters->GetAllFloatNames();
	vector<double> truth;
	for( unsigned short int j=0; j < other_params.size(); j++ )
	{
		truth.push_back( BottleParameters->GetPhysicsParameter(other_params[j])->GetTrueValue() );
	}


	bool fit_fail_status=false;
	FitResult* ReturnableFitResult;
	//  Try to fit 5 times and then abort
	for( unsigned short int i=0; i<=0; i++ )
	{
		//  Left in for correctness
		fit_fail_status = false;

		try{
			// Try a Fit, it it converges, continue to elsewhere in the program
			ReturnableFitResult = FitAssembler::DoFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints );
			if( ReturnableFitResult->GetFitStatus() != 3 )
			{
				cerr << "\n\n\t\tFit Did NOT Converge Correctly, CHECK YOUR RESULTS!\n\n";
			}
			return ReturnableFitResult;
		}
		//  If it didn't fit tell the user why!
		catch( int e){
			if ( e == 10 ){
				cerr << "\nCaught exception : fit failed for these parameters..." << endl;  }
			else if ( e == 13 ){
				cerr << "\nIntegration Error: Fit Failed..." << endl;  }
			fit_fail_status = true;
		} catch (...) {
                        cerr << "\n\n\n\t\t\tCaught Unknown Exception, THIS IS SERIOUS!!!\n\n\n" << endl;
			fit_fail_status = true;
		}
		//cerr << "Fit Did Not converge, shaking the bottle and starting again!" <<endl;
		//cerr << "This is Retry " << i+1 << "of 4"<<endl;
		//  Give the physics bottle a bit of a shake and see if it works this time
		//ShakeBottle( BottleParameters, BottleData, (i+5) );
	}

	//  If the fit failed 5 times I will simply return a dummy fit result full of zerod objects. It is up to the user to watch for and remove these
	if( fit_fail_status ){
		cerr << "Nothing more I'm willing to do, considering a Fit here a lost cause..." <<endl;

		for( unsigned short int j=0; j < other_params.size(); j++ )
		{
			BottleParameters->GetPhysicsParameter( other_params[j] )->SetTrueValue( truth[j] );
		}
		int status = -1;
		vector<string> NewNamesList = BottleParameters->GetAllNames();
		ResultParameterSet* DummyFitResults = new ResultParameterSet( NewNamesList );
		ReturnableFitResult = new FitResult( LLSCAN_FIT_FAILURE_VALUE, DummyFitResults, status, BottleParameters );
	}
	
	return ReturnableFitResult;
}


//  Interface for internal calls
void FitAssembler::DoScan( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters, vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, ScanParam* Wanted_Param, ToyStudyResult* output_interface )
{

	double uplim = Wanted_Param->GetMax();
	double lolim = Wanted_Param->GetMin();
	double npoints = Wanted_Param->GetPoints();
	string scanName = Wanted_Param->GetName();

//	cout << "Performing Scan for the parameter " << scanName << endl ;
	
	// Get a pointer to the physics parameter to be scanned and fix it	
	// CAREFUL:  this must be reset as it was at the end.
	PhysicsParameter * scanParameter = BottleParameters->GetPhysicsParameter(scanName);
	double originalValue = scanParameter->GetBlindedValue( ) ;
	string originalType = scanParameter->GetType( ) ;
	scanParameter->SetType( "Fixed" ) ;

	// Need to set up a loop , fixing the scan parameter at each point
	double deltaScan;

	for( int si=0; si<int(npoints); si++) {
		cout << "\n\nSINGLE SCAN NUMBER\t\t" << si+1 << "\t\tOF\t\t" <<int(npoints)<< endl<<endl;
		// Set scan parameter value
		if( int(npoints)!=1 ) deltaScan = (uplim-lolim) / (npoints-1.) ;
		else deltaScan=0;
		double scanVal = lolim+deltaScan*si;
		scanParameter->SetBlindedValue( scanVal ) ;

		output_interface->StartStopwatch();

		// Do a scan point fit
		FitResult * scanStepResult = FitAssembler::DoSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints );

		//  THIS IS ALWAYS TRUE BY DEFINITION OF THE SCAN
		string name = Wanted_Param->GetName();
		string type = BottleParameters->GetPhysicsParameter( name )->GetType();
		string unit = BottleParameters->GetPhysicsParameter( name )->GetUnit();
		scanStepResult->GetResultParameterSet()->SetResultParameter( name, scanVal, 0, 0., scanVal, scanVal, type, unit );

		vector<string> Fixed_List = BottleParameters->GetAllFixedNames();
		vector<string> Fit_List = scanStepResult->GetResultParameterSet()->GetAllNames();
		for( unsigned short int i=0; i < Fixed_List.size() ; i++ )
		{
			bool found=false;
			for( unsigned short int j=0; j < Fit_List.size(); j++ )
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

	double deltaScan;
	if( int(npoints) !=1 ) deltaScan = (uplim-lolim) / (npoints-1.) ;
	else deltaScan=0.;

	for( int si=0; si < int(npoints); si++) {
	  
		cout << "\n\n2DSCAN OUTER NUMBER\t\t" << si+1 << "\t\tOF\t\t" << int(npoints) <<endl<<endl;
		ToyStudyResult* Returnable_Result = new ToyStudyResult( result_names );
		
		// Set scan parameter value
		double scanVal = lolim + si*deltaScan;
		scanParameter->SetBlindedValue( scanVal ) ;

		// Do a scan point fit
		FitAssembler::DoScan( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, Param_Set.second, Returnable_Result );

		//  THIS IS ALWAYS TRUE BY DEFINITION OF THE SCAN
                string name = Param_Set.first->GetName();
                string type = BottleParameters->GetPhysicsParameter( name )->GetType();
                string unit = BottleParameters->GetPhysicsParameter( name )->GetUnit();

		for( short int i=0; i < Returnable_Result->NumberResults(); i++ )
		{
			Returnable_Result->GetFitResult( i )->GetResultParameterSet()->SetResultParameter( name, scanVal, 0, 0.0, scanVal, scanVal, type, unit );
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
ToyStudyResult* FitAssembler::FeldmanCousins( ToyStudyResult* GlobalResult, ToyStudyResult* _2DResultForFC, vector<unsigned short int> numberRepeats, unsigned short int NuisenceModel, bool FC_Debug_Flag, OutputConfiguration* makeOutput, MinimiserConfiguration* theMinimiser, FitFunctionConfiguration* theFunction, XMLConfigReader* xmlFile, vector< PDFWithData* > pdfsAndData )
{
	FitResult* GlobalFitResult = GlobalResult->GetFitResult( 0 );
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

		ParameterSet* InputParamSet = NULL;
		ParameterSet* InputFreeSet = NULL;
		ParameterSet* ControlParamSet = NULL;

		//	True by definition:
		double lim1 = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetResultParameter( name1 )->GetValue();
		double lim2 = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetResultParameter( name2 )->GetValue();

		if( NuisenceModel == 1 )
		{
			//		Use the inputs from Step 1
			InputParamSet = GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet();
			InputFreeSet = GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet();
			ControlParamSet = GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet();
		} else if( NuisenceModel == 2 ){
			//		Use the inputs from Step 2
			InputParamSet = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet();
			InputFreeSet = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet();
			ControlParamSet = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet();
		}


		//		Just for clarity
		//		also when using values from Step1 these are not set correctly for the study
		ControlParamSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
		ControlParamSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
		InputParamSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
		InputParamSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
		InputParamSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
		InputParamSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
		InputParamSet->GetPhysicsParameter( name1 )->SetType( "Fixed" );
		InputParamSet->GetPhysicsParameter( name2 )->SetType( "Fixed" );
		InputFreeSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
		InputFreeSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
		InputFreeSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
		InputFreeSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
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
		unsigned short int wanted_number_of_toys = 100;
		if( !numberRepeats.empty() ) wanted_number_of_toys = numberRepeats[0];



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
		ToyStudyResult* study1Results = new ToyStudyResult( GlobalFitResult->GetResultParameterSet()->GetAllNames() );
		ToyStudyResult* study2Results = new ToyStudyResult( GlobalFitResult->GetResultParameterSet()->GetAllNames() );


		//	Try to make minuit as quiet as possible here... does someone have a way to gag it's output?!?
		ToyStudyMinimiser->GetMinimiser( int(InputParamSet->GetAllNames().size()) )->SetOutputLevel(-1);


		//	This Forms the main sequence of the fit

		cout << "\n\n\t\tPerforming Fits to Toys in FC\n"<<endl;
		//		Perform Fits
		//	This will record ONLY Data from working fits i.e. FitStatus==3 unless the Debug flag is on
		for( unsigned short int dataset_num=0; dataset_num < wanted_number_of_toys; dataset_num++ )
		{

			bool toy_failed = false;

			FitResult* fit1Result = NULL;
			FitResult* fit2Result = NULL;

			ParameterSet* LocalInputFreeSet=NULL;
			ParameterSet* LocalInputFixedSet=NULL;

			if( NuisenceModel == 1 )	{
				//	Assuming Nuisence Parameters set to Global Minima
				//	We need to set this for EVERY FIT in order to have correct generation/pull values
				LocalInputFreeSet = GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet();
				LocalInputFixedSet = GlobalFitResult->GetResultParameterSet()->GetDummyParameterSet();
			} else if( NuisenceModel == 2 )	{
				//	Assuming Nuisence Parameters set to Local Minima
				//
				LocalInputFreeSet = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet();
				LocalInputFixedSet = _2DResultForFC->GetFitResult( iFC )->GetResultParameterSet()->GetDummyParameterSet();
			} else if( NuisenceModel == 3 )	{
				//	Generate a ParameterSet with random numbers based on
				//	Fit Results and assign it to the ControlParamSet which is used for generation
				//	Also need to Setup Local ParameterSets as with standard case
			}

			ControlParamSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
			ControlParamSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
			LocalInputFreeSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
			LocalInputFreeSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
			LocalInputFreeSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
			LocalInputFreeSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
			LocalInputFreeSet->GetPhysicsParameter( name1 )->SetType( "Free" );
			LocalInputFreeSet->GetPhysicsParameter( name2 )->SetType( "Free" );

			LocalInputFixedSet->GetPhysicsParameter( name1 )->SetBlindedValue( lim1 );
			LocalInputFixedSet->GetPhysicsParameter( name1 )->ForceOriginalValue( lim1 );
			LocalInputFixedSet->GetPhysicsParameter( name2 )->SetBlindedValue( lim2 );
			LocalInputFixedSet->GetPhysicsParameter( name2 )->ForceOriginalValue( lim2 );
			LocalInputFixedSet->GetPhysicsParameter( name1 )->SetType( "Fixed" );
			LocalInputFixedSet->GetPhysicsParameter( name2 )->SetType( "Fixed" );

			//	We need to set some factors before we perform the fit
			for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); pdf_num++ )
			{
				vector<IDataSet*> wanted_set;
				wanted_set.push_back( Memory_Data[pdf_num][dataset_num] );
				PDFsWithDataForToys[pdf_num]->AddCachedData( wanted_set );
				PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( LocalInputFreeSet );
			}

			cout << "\n\n\t\tPerforming Fit To Toy: "<< (dataset_num+1) <<" of " << wanted_number_of_toys << endl<<endl;
			//	Fit once with control parameters Free
			fit1Result = FitAssembler::DoSafeFit( ToyStudyMinimiser, ToyStudyFunction, LocalInputFreeSet, PDFsWithDataForToys, ConstraintsForToys );

			//	Only Fit again to this dataset if it fits well with +2 dof
			//	This has the obvious savings in CPU resources
			//	We may want to keep this information so leaving it configurable.
			if( ( fit1Result->GetFitStatus() == 3 ) || FC_Debug_Flag )
			{
				if( FC_Debug_Flag )	cout << "\n\nYou are aware your requesting all output?\n"<<endl;

				//  Fit secondGlobalResult with control parameters Fixed
				for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); pdf_num++ )
				{
					vector<IDataSet*> wanted_set;
					wanted_set.push_back( Memory_Data[pdf_num][dataset_num] );
					PDFsWithDataForToys[pdf_num]->AddCachedData( wanted_set );
					PDFsWithDataForToys[pdf_num]->SetPhysicsParameters( LocalInputFixedSet );
				}

				cout << "\n\n\t\tFirst Fit Successful, Performing the Second Fit " << (dataset_num+1) << " of " << wanted_number_of_toys <<endl;
				//	Use the SafeFit as this always returns something when a PDF has been written to throw not exit
				fit2Result = FitAssembler::DoSafeFit( ToyStudyMinimiser, ToyStudyFunction, LocalInputFixedSet, PDFsWithDataForToys, ConstraintsForToys );

				//	If either Fit Failed we want to 'dump the results' and run an extra Fit.
				if( (fit2Result->GetFitStatus() != 3) || (fit2Result->GetFitStatus() != 3) ) toy_failed = true;
			}
			else	toy_failed = false;

			//	Do we want to store the Data OR run another toy to get a better Fit
			if( toy_failed )
			{
				cerr << "\n\n\t\tA Single Toy Study Failed... Requesting More Data for another pass.\n" << endl;
				//  Increment counter so we can guarantee we get 'wanted_number_of_toys' which have fitted correctly
				wanted_number_of_toys++;
			}

			//	In my opinion the fact this is so easy and doesn't cause root to throw up is a testiment to the other authors in RapidFit :D
			cout << "\n\n\t\tDeleting Used Data:\n"<<endl;
			for( unsigned int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); pdf_num++ )
			{
				cout << "deleting data for pdf: " << pdf_num << "\n";
				delete Memory_Data[pdf_num].back();
				Memory_Data[pdf_num].back() = NULL;
			}

			// Only generate data if I'm going to fit to it
			if( dataset_num < (wanted_number_of_toys-1) )
			{
				cout << "\n\n\t\tGenerating Data For Next Toy:\n" <<endl;

				for( unsigned short int pdf_num=0; pdf_num < PDFsWithDataForToys.size(); pdf_num++ )
				{
					// See above for more detailed description
					cout << "generating data for pdf: " << pdf_num << "\n";
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
		ToyStudyResult* ThisStudy = new ToyStudyResult( GlobalFitResult->GetResultParameterSet()->GetAllNames() );
		ThisStudy->AddFitResult( _2DResultForFC->GetFitResult( iFC ), false );
		ThisStudy->AddCPUTimes( _2DResultForFC->GetAllCPUTimes() );
		ThisStudy->AddRealTimes( _2DResultForFC->GetAllRealTimes() );
		//	The Generated Value for the Global and Local fit are best defined as -9999 as a sensible default
		for( unsigned short int num=0; num < GlobalFitResult->GetResultParameterSet()->GetAllNames().size(); num++ )
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
	}

	return new ToyStudyResult( AllResults );
}
