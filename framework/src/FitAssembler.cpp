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
#include "PhysicsBottle.h"
///	System Headers
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <signal.h>

#ifdef __USE_VALGRIND
#include <valgrind/callgrind.h>
#endif

using namespace::std;

//	We will catch any throw statments internal to the minimisation process
void FitAssembler::SafeMinimise( IMinimiser* Minimiser )
{
#ifdef __USE_VALGRIND
	CALLGRIND_START_INSTRUMENTATION;
#endif
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
#ifdef __USE_VALGRIND
	CALLGRIND_STOP_INSTRUMENTATION;
	CALLGRIND_DUMP_STATS;
#endif
}

//The final stage - do the minimisation
FitResult * FitAssembler::DoFit( IMinimiser * Minimiser, IFitFunction * TheFunction )
{
	Minimiser->SetupFit( TheFunction );

	cout << "Using Fixed Parameters:" << endl;
	ParameterSet* fittingSet = TheFunction->GetParameterSet();
	vector<string> FixedNames = fittingSet->GetAllFixedNames();
	for( unsigned int i=0; i< FixedNames.size(); ++i )
	{
		cout << setw(25) << left << FixedNames[i] << fittingSet->GetPhysicsParameter( FixedNames[i] )->GetValue() << endl;
	}

	cout << "\nStarting Fit!" << endl;

	SafeMinimise( Minimiser );

	cout << "\nMinimised!\n" << endl;

	FitResult* final_result = Minimiser->GetFitResult();

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Returning FitResult" << endl;
	}

	/*
	//      Now to remove any caches created from the current PDFs on data
	PhysicsBottle* thisBottle = TheFunction->GetPhysicsBottle();
	for( int i=0; i< thisBottle->NumberResults(); ++i )
	{
		IDataSet* thisDataSet = thisBottle->GetResultDataSet( i );
	}
	*/
	return final_result;
}

//Create the minimiser and fit function
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, PhysicsBottle * Bottle )
{

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Testing if Fixed ParameterSet" << endl;
	}

	if( Bottle->GetParameterSet()->GetAllFloatNames().size() == 0 )
	{
		vector<string> ParamNames = Bottle->GetParameterSet()->GetAllFixedNames();
		ResultParameterSet* fixed_set = new ResultParameterSet( ParamNames );
		for( unsigned int i=0; i< ParamNames.size(); ++i )
		{
			PhysicsParameter* param = Bottle->GetParameterSet()->GetPhysicsParameter( ParamNames[i] );
			fixed_set->SetResultParameter( ParamNames[i], param->GetValue(), param->GetValue(), 0., 0., 0., param->GetType(), param->GetUnit() );
		}

		if( DebugClass::DebugThisClass( "FitAssembler" ) )
		{
			cout << "FitAssembler: Setting Fixed Physics Bottle" << endl;
		}

		IFitFunction* theFunction = FunctionConfig->GetFitFunction();

		vector<IDataSet*> allDataSets = Bottle->GetAllDataSets();
		for( unsigned int i=0; i< allDataSets.size(); ++i )
		{
			unsigned int num= (unsigned)allDataSets[i]->GetDataNumber();
			for( unsigned int j=0; j< num; ++j )
			{
				allDataSets[i]->GetDataPoint( (int)j )->ClearDiscreteIndexMap();
			}
		}

		theFunction->SetPhysicsBottle(Bottle);

		if( DebugClass::DebugThisClass( "FitAssembler" ) )
		{
			cout << "FitAssembler: Evaluating DataSet" << endl;
		}

		double someTest = theFunction->Evaluate();
		cout << someTest << endl;
		//exit(0);

		if( DebugClass::DebugThisClass( "FitAssembler" ) )
		{
			cout << "FitAssembler: Constructing Fixed Result" << endl;
		}

		FitResult* fixed_result = new FitResult( 0., fixed_set, 3, Bottle );

		if( DebugClass::DebugThisClass( "FitAssembler" ) )
		{
			cout << "FitAssembler: Setting Result Minimum Value" << endl;
		}

		delete fixed_set;
		fixed_result->SetMinimumValue(theFunction->Evaluate());
		return fixed_result;
	}
	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Constructing Minimiser" << endl;
	}

	IMinimiser * minimiser = MinimiserConfig->GetMinimiser( int(Bottle->GetParameterSet()->GetAllNames().size()) );
	cout << endl;

	/*
	   vector<IPDF*> pdfs = Bottle->GetAllPDFs();
	   for( unsigned int i=0; i< pdfs.size(); ++i )
	   {
	   pdfs[i]->GetPhysicsParameters()->Print();
	   }
	 */

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Constructing FitFunction" << endl;
	}

	IFitFunction* theFunction = FunctionConfig->GetFitFunction();

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler:  About to pass PhysicsBottle to FitFunction and Test Integrator" << endl;
	}

	theFunction->SetPhysicsBottle(Bottle);

	cout << endl;

	FitResult* result = FitAssembler::DoFit( minimiser, theFunction );

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Destroying FitFunction and returning FitResult" << endl;
	}
	delete theFunction;

	return result;
}

//Create the physics bottle
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet* BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints )
{
	//double someVal;
	vector<IPDF*> allPDFs;
	for( unsigned int i=0; i< BottleData.size(); ++i )
	{
		allPDFs.push_back( BottleData[i]->GetPDF() );
	}

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: For PDFs:" << endl;
		for( unsigned int i=0; i< allPDFs.size(); ++i ) cout << allPDFs[i]->GetLabel() << endl;
		cout << "FitAssembler: Constructing PhaseSpace:" << endl;
	}

	vector<int> allDataNum;
	for( unsigned int i=0; i< BottleData.size(); ++i )	allDataNum.push_back( BottleData[i]->GetDataSetConfig()->GetDataSetSize() );

	ParameterSet* checkedBottleParameters = FitAssembler::CheckInputParams( BottleParameters, allPDFs, allDataNum );

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		checkedBottleParameters->Print();
		cout << "FitAssembler: Sorting Parameters to get Floated Parameters first" << endl;
	}

	vector<string> Wanted_List_OF_Floated_Parameters = FunctionConfig->GetFloatedParameterList();

	checkedBottleParameters->FloatedFirst( Wanted_List_OF_Floated_Parameters );

	ParameterSet* checkedGenerationParameters = GenerationParameters( checkedBottleParameters, BottleParameters );

	PhysicsBottle * bottle = new PhysicsBottle( checkedBottleParameters );

	//Fill the bottle - data generation occurs in this step
	for( unsigned int resultIndex = 0; resultIndex < BottleData.size(); ++resultIndex )
	{
		//	Use the Raw input as the Generation PDF may require Parameters not involved in the main fit
		if( DebugClass::DebugThisClass( "FitAssembler" ) )
		{
			cout << "Setting Physics Parameters in Bottle" << endl;
		}
		BottleData[resultIndex]->SetPhysicsParameters( checkedGenerationParameters );

		IPDF* Requested_PDF = BottleData[resultIndex]->GetPDF();

		IPDF* genPDF = BottleData[resultIndex]->GetDataSetConfig()->GetGenerationPDF();

		if( genPDF != NULL )
		{
			if( DebugClass::DebugThisClass( "FitAssembler" ) )
			{
				cout << "FitAssembler: Checking for all required Generation PDFs" << endl;
			}

			ParameterSet* checkedSet = CheckInputParams( checkedGenerationParameters, vector<IPDF*>(1,genPDF), vector<int>(1,1) );

			genPDF->UpdatePhysicsParameters( checkedSet );

			delete checkedSet;
		}

		if( DebugClass::DebugThisClass( "FitAssembler" ) )
		{
			cout << "FitAssembler: Generating DataSet: " << resultIndex+1 << " of: " << BottleData.size() << endl;
		}

		IDataSet* Requested_DataSet = BottleData[resultIndex]->GetDataSet();

		if( FunctionConfig->GetWeightsWereUsed() )
		{
			vector<IPDF*> thisPDF = vector<IPDF*>( 1, BottleData[resultIndex]->GetPDF() );
			vector<IDataSet*> thisData = vector<IDataSet*>( 1, Requested_DataSet );
			string weightName, alphaName;
			if( FunctionConfig->GetWeightsWereUsed() )	weightName = FunctionConfig->GetWeightName();
			if( FunctionConfig->GetUseCustomAlpha() )	alphaName = FunctionConfig->GetAlphaName();
			FitAssembler::CheckInputObs( thisPDF, thisData, weightName, alphaName );
			Requested_DataSet->UseEventWeights( FunctionConfig->GetWeightName() );
			Requested_DataSet->PrintYield();
		}

		if( Requested_DataSet->GetDataNumber() > 0 )
		{

			if( FunctionConfig->GetNormaliseWeights() && Requested_DataSet->GetWeightsWereUsed() ) Requested_DataSet->NormaliseWeights();

			if( DebugClass::DebugThisClass( "FitAssembler" ) )
			{
				cout << "FitAssembler: Adding PDF & DataSet to Bottle: " << resultIndex+1 << " of: " << BottleData.size() << endl;
			}

			bottle->AddResult( Requested_PDF, Requested_DataSet );
		}
		else
		{
			if( DebugClass::DebugThisClass( "FitAssembler" ) )
			{
				cout << "FitAssembler: ***NOT*** Adding PDF & DataSet to Bottle: " << resultIndex+1 << " of: " << BottleData.size() << endl;
				cout << "FitAssembler: DataSet Size is <= 0!!" << endl;
			}
		}
	}

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Checking Weights and Normalisation" << endl;
	}

	if( FunctionConfig->GetSingleNormaliseWeights() )
	{
		double sum_total=0.;
		double sum_sq_total=0.;
		for( unsigned int i=0; i< (unsigned)bottle->NumberResults(); ++i )
		{
			IDataSet* Requested_DataSet = bottle->GetResultDataSet( (int)i );
			if( Requested_DataSet->GetWeightsWereUsed() )
			{
				sum_total += Requested_DataSet->GetSumWeights();
				sum_sq_total += Requested_DataSet->GetSumWeightsSq();
			}
		}
		for( unsigned int i=0; i< BottleData.size(); ++i )
		{
			IDataSet* Requested_DataSet = bottle->GetResultDataSet( (int)i );
			if( Requested_DataSet->GetWeightsWereUsed() )
			{
				Requested_DataSet->ApplyAlpha( sum_total, sum_sq_total );
				//Requested_DataSet->SortBy( "time" );
			}
			//Requested_DataSet->SortBy( "time" );
		}
	}

	if( FunctionConfig->GetUseCustomAlpha() )
	{
		for( unsigned int i=0; i< BottleData.size(); ++i )
		{
			IDataSet* Requested_DataSet = bottle->GetResultDataSet( (int)i );
			Requested_DataSet->ApplyExternalAlpha( FunctionConfig->GetAlphaName() );
		}
	}

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Adding Constraints" << endl;
	}

	//Add the constraints
	for ( unsigned int constraintIndex = 0; constraintIndex < BottleConstraints.size(); ++constraintIndex )
	{
		bottle->AddConstraint( BottleConstraints[constraintIndex] );
	}

	delete checkedBottleParameters;

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Passing MinimiserConfig, FunctionConfig and Bottle" << endl;
	}

	vector<IPDF*> thesePDFs = bottle->GetAllPDFs();
	vector<IDataSet*> theseDataSets = bottle->GetAllDataSets();

	string weightName, alphaName;
	if( FunctionConfig->GetWeightsWereUsed() ) weightName = FunctionConfig->GetWeightName();
	if( FunctionConfig->GetUseCustomAlpha() ) alphaName = FunctionConfig->GetAlphaName();
	FitAssembler::CheckInputObs( thesePDFs, theseDataSets, weightName, alphaName );

	FitResult * result = FitAssembler::DoFit( MinimiserConfig, FunctionConfig, bottle );

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Adding all PhysicsParameters from XML to output" << endl;
	}

	BottleParameters->AddPhysicsParameters( bottle->GetParameterSet(), true );

	vector<string> fixed_names = checkedGenerationParameters->GetAllFixedNames();

	for( vector<string>::iterator name_i = fixed_names.begin(); name_i != fixed_names.end(); ++name_i )
	{
		PhysicsParameter* thisParam = checkedGenerationParameters->GetPhysicsParameter( *name_i );

		ResultParameter* thisResult = new ResultParameter( thisParam );

		bool wasAdded = result->GetResultParameterSet()->ForceNewResultParameter( thisResult );

		if( !wasAdded ) result->GetResultParameterSet()->SetResultParameter( *name_i, thisResult );

		if( wasAdded )
		{
			if( result->GetResultParameterSet()->GetResultParameter( *name_i )->GetType() == "Fixed" )
			{
				result->GetResultParameterSet()->GetResultParameter( *name_i )->SetError(-2.);
			}
			else
			{
				result->GetResultParameterSet()->GetResultParameter( *name_i )->SetError(-1.);
			}
		}

		delete thisResult;
	}

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Destroying Bottle and returning FitResult" << endl;
	}

	delete checkedGenerationParameters;

	delete bottle;

	return result;
}

void FitAssembler::CheckInputObs( const vector<IPDF*> AllPDFs, const vector<IDataSet*> allData, const string WeightName, const string AlphaName )
{
	if( AllPDFs.size() != allData.size() )
	{
		cerr << "FitAssembler: Number of User PDFs NOT the Same as the Number of DataSets being used to minimise" << endl;
		cerr << endl;
		exit(0);
	}
	else
	{
		bool missing_any_all=false;
		for( unsigned int i=0; i< allData.size(); ++i )
		{

			bool missing_any=false;

			ObservableRef WeightRef(WeightName);
			ObservableRef AlphaRef(AlphaName);

			try
			{
				if( allData[i]->GetDataNumber() != 0 )
				{

					if( !WeightName.empty() )
					{
						Observable* WeightObs = allData[i]->GetDataPoint( 0 )->GetObservable( WeightRef );
						if( WeightObs == NULL || WeightRef.GetIndex() < 0 )
						{
							missing_any = true;
							cerr << "FitAssembler: Couldn't find Observable " << WeightName << " in DataPoint 0!" << endl;
							cerr << "FitAssembler: Please check it is in the DataSet for the PDF to access it!" << endl;
							cerr << "FitAssembler: This is likely being requested by your FitFunctionConfiguration!" << endl;
							cerr << endl;
						}
					}

					if( !AlphaName.empty() )
					{
						Observable* AlphaObs = allData[i]->GetDataPoint( 0 )->GetObservable( AlphaRef );
						if( AlphaObs == NULL || AlphaRef.GetIndex() < 0 )
						{
							missing_any = true;
							cerr << "FitAssembler: Couldn't find Observable " << AlphaName << " in DataPoint 0!" << endl;
							cerr << "FitAssembler: Please check it is in the DataSet for the PDF to access it!" << endl;
							cerr << "FitAssembler: This is likely being requested by your FitFunctionConfiguration!" << endl;
							cerr << endl;
						}
					}
				}
			}
			catch(...)
			{
				missing_any = true;
				cerr << "FitAssembler: Error Getting WeightName or AlphaName: \"" << WeightName << "\" , \"" << AlphaName << "\"" << endl;
				cerr << "FitAssembler: This is likely due to your FitFunction Requirements!" << endl;
				cerr << endl;
			}

			if( allData[i] == NULL )
			{
				cerr << "FitAssembler: NULL Pointer to DataSet Please Fix!!!" << endl;
				exit(0);
			}
			if( AllPDFs[i] == NULL )
			{
				cerr << "FitAssembler: NULL Pointer to PDF Please Fix!!!" << endl;
				exit(0);
			}

			vector<ObservableRef> wantedObservables;
			for( unsigned int j=0; j< AllPDFs[i]->GetPrototypeDataPoint().size(); ++j )
			{
				wantedObservables.push_back( ObservableRef( AllPDFs[i]->GetPrototypeDataPoint()[j] ) );
			}
			for( unsigned int j=0; j< wantedObservables.size(); ++j )
			{

				try
				{
					if( allData[i]->GetDataNumber() != 0 )
					{
						Observable* thisObs = allData[i]->GetDataPoint( 0 )->GetObservable( wantedObservables[j] );

						if( thisObs == NULL || wantedObservables[j].GetIndex() < 0 )
						{
							missing_any = true;
							cerr << "FitAssembler: Couldn't find Observable " << string( wantedObservables[j] ) << " in DataPoint 0!" << endl;
							cerr << "FitAssembler: Please check it is in the DataSet for the PDF to access it!" << endl;
							cerr << endl;
						}
					}
				}
				catch(...)
				{
					missing_any = true;
					cerr << "FitAssembler: Fell over looking for Observable  " << string( wantedObservables[j] ) << " in DataPoint 0!" << endl;
					cerr << "FitAssembler: Please check it is in the DataSet for the PDF to access it!" << endl;
					cerr << endl;
				}

			}
			for( unsigned int j=0; j< wantedObservables.size(); ++j )
			{
				try
				{
					PhaseSpaceBoundary* thisBoundary = allData[i]->GetBoundary();
					IConstraint* thisConstraint = thisBoundary->GetConstraint( wantedObservables[j] );

					if( thisConstraint == NULL || wantedObservables[j].GetIndex() < 0 )
					{
						missing_any = true;
						cerr << "FitAssembler: Couldn't find a PhaseSpace Constraint Observable " << string( wantedObservables[j] ) << "!" << endl;
						cerr << "FitAssembler: Please check it is in the PhaseSpaceBoundary!" << endl;
						cerr << endl;
					}
				}
				catch(...)
				{
					missing_any = true;
					cerr << "FitAssembler: Fell over while finding a PhaseSpace Constraint Observable " << string( wantedObservables[j] ) << "!" << endl;
					cerr << "FitAssembler: Please check it is in the PhaseSpaceBoundary!" << endl;
					cerr << endl;
				}
			}
			if( missing_any ) cerr << "Failure in DataSet: " << i+1 << " of " << allData.size() << endl;

			missing_any_all = missing_any_all||missing_any;
		}
		if( missing_any_all ) exit(0);
	}
}

//Check that the provided ParameterSet only Contains the Parameters claimed by the PDFs to protect the Minimiser from runtime mistakes
ParameterSet* FitAssembler::CheckInputParams( const ParameterSet* givenParams, const vector<IPDF*> allPDFs, const vector<int> allDataNum )
{
	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Checking ParameterSet from XML" << endl;
	}
	vector<string> param_names;
	for( unsigned int i=0; i< allPDFs.size(); ++i )
	{
		if( allDataNum[i] > 0 )
		{
			vector<string> input_from_PDF = allPDFs[i]->GetPrototypeParameterSet();
			param_names = StringProcessing::CombineUniques( param_names, input_from_PDF );
		}
	}

	if( param_names.size() != givenParams->GetAllNames().size() )
	{
		vector<string> temp = givenParams->GetAllNames();
	}

	bool caught = false;
	ParameterSet* wantedParameterSet = ( new ParameterSet( param_names ) );
	for( unsigned int i=0; i< param_names.size(); ++i )
	{
		PhysicsParameter* phys_param = NULL;
		try
		{
			phys_param = givenParams->GetPhysicsParameter( param_names[i] );
		}
		catch(...)
		{
			cout << endl;
			if( param_names[i] == "UnNamed" )
			{
				cout << "You haven't Named this Parameter, please check all PhysicsParameter-s are named correctly in your PDF Constructor" << endl << endl;
			}
			else
			{
				cout << "FitAssembler: You are missing the Parameter: " << param_names[i] << " in your XML File. Please add it to the ParameterSet!" << endl << endl;
				cout << "Eg:" << endl;
				PhysicsParameter* temp = new PhysicsParameter( param_names[i] );
				string xmlStr = temp->XML();
				delete temp;
				cout << xmlStr << endl << endl;
			}
			caught = true;
		}
		if( !caught ) wantedParameterSet->SetPhysicsParameter( param_names[i], new PhysicsParameter(*phys_param) );
	}

	if( caught ) exit(0);

	return wantedParameterSet;
}

ParameterSet* FitAssembler::GenerationParameters( const ParameterSet* checkedBottleParameters, const ParameterSet* BottleParameters )
{
	vector<string> allForGeneration = BottleParameters->GetAllNames();

	ParameterSet* GenerationParameterSet = new ParameterSet( *checkedBottleParameters );

	for( vector<string>::iterator param_i = allForGeneration.begin(); param_i != allForGeneration.end(); ++param_i )
	{
		PhysicsParameter* thisParam = BottleParameters->GetPhysicsParameter( *param_i );
		GenerationParameterSet->AddPhysicsParameter( thisParam );
	}

	return GenerationParameterSet;
}

//Create the physics bottle with pre-made data
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet* BottleParameters,
		vector< IPDF* > AllPDFs, vector< IDataSet* > AllData, vector< ConstraintFunction* > BottleConstraints )
{
	vector<int> allDataNum;
	for( unsigned int i=0; i< AllData.size(); ++i )	allDataNum.push_back( AllData[i]->GetDataNumber() );

	ParameterSet* internalBottleParameters = FitAssembler::CheckInputParams( BottleParameters, AllPDFs, allDataNum );

	FitAssembler::CheckInputObs( AllPDFs, AllData, FunctionConfig->GetWeightName(), FunctionConfig->GetAlphaName() );

	if ( AllPDFs.size() == AllData.size() )
	{
		PhysicsBottle * bottle = new PhysicsBottle( internalBottleParameters );

		//Fill the bottle - data already generated
		for ( unsigned int resultIndex = 0; resultIndex < AllData.size(); ++resultIndex )
		{
			if( allDataNum[resultIndex] > 0 )
			{
				AllPDFs[resultIndex]->UpdatePhysicsParameters(internalBottleParameters);
				bottle->AddResult( AllPDFs[resultIndex], AllData[resultIndex] );
			}
		}

		//Add the constraints
		for ( unsigned int constraintIndex = 0; constraintIndex < BottleConstraints.size(); ++constraintIndex )
		{
			bottle->AddConstraint( BottleConstraints[constraintIndex] );
		}

		//bottle->Finalise();
		FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );

		delete bottle;
		delete internalBottleParameters;

		BottleParameters->AddPhysicsParameters( bottle->GetParameterSet(), true );

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
	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Checking Result ParameterSet" << endl;
	}
	vector<string> already_found = ReturnableFitResult->GetResultParameterSet()->GetAllNames();
	vector<string> bottle_names = BottleParameters->GetAllNames();

	for( unsigned int i=0; i< bottle_names.size(); ++i )
	{
		int found = StringProcessing::VectorContains( &already_found, &(bottle_names[i]) );

		PhysicsParameter* Phys_Param = BottleParameters->GetPhysicsParameter( bottle_names[i] );

		//	There was something in the ParameterSet not in the FitResult, i.e. an unclaimed object which can't have changed during the fit
		if( found == -1 )
		{
			//cout << "ALERT:\t" << "Parameter " << BottleParameters->GetAllNames()[i] << " was not claimed by any PDF in the fit and was NOT passed to the Minimiser!!!" << endl;
			double Value = Phys_Param->GetValue();
			double OriginalValue = Phys_Param->GetValue();
			double Minimum = Phys_Param->GetMinimum();
			double Maximum = Phys_Param->GetMaximum();
			//double Error = Phys_Param->GetStepSize();
			string Type = Phys_Param->GetType();
			string Unit = Phys_Param->GetUnit();
			bool added = ReturnableFitResult->GetResultParameterSet()->ForceNewResultParameter( bottle_names[i],  Value, OriginalValue, -1., Minimum, Maximum, Type, Unit );
			if( !added )
			{
				cerr << "Error finalizing FitResultVector Object" << endl << endl;
				exit(-984);
			}
		}
		else
		{
			if( Phys_Param->GetType() != "Fixed" )
			{
				ResultParameter* result_Param = ReturnableFitResult->GetResultParameterSet()->GetResultParameter( already_found[(unsigned)found] );
				double result_val = result_Param->GetValue();
				double result_err = result_Param->GetError();
				double param_min = Phys_Param->GetMinimum();
				double param_max = Phys_Param->GetMaximum();
				bool nolim = fabs( param_min - param_max ) > 1E-6;
				if( (fabs( result_val - param_min ) <= (0.5 * result_err)) && nolim )
				{
					cout << "Caution: Parameter " << bottle_names[i] << " is within 0.5 sigma of a limit!!!" << endl;
					cout << "\t\t" << bottle_names[i] << ": " << result_val << " Min: " << param_min << endl << endl;
				}
				if( (fabs( result_val - param_max ) <= (0.5 * result_err)) && nolim )
				{
					cout << "Caution: Parameter " << bottle_names[i] << " is within 0.5 sigma of a limit!!!" << endl;
					cout << "\t\t" << bottle_names[i] << ": " << result_val << " Max: " << param_max << endl << endl;
				}
			}
		}
	}
}

//  Perform a safer fit which should always to return something which you can use :D
FitResult * FitAssembler::DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet* BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, bool forceContinue, int OutputLevel )
{
	FitResult* ReturnableFitResult=NULL;
	ParameterSet* internalBottleParameters = new ParameterSet( *BottleParameters );

	//	Decide what 'Strategy' to use By default just perform the fit
	if( FunctionConfig->GetStrategy() == "Petes" )
	{
		ReturnableFitResult = Petes_DoSafeFit( MinimiserConfig, FunctionConfig, internalBottleParameters, BottleData, BottleConstraints, forceContinue, OutputLevel );
	}
	else if( FunctionConfig->GetStrategy() == "PetesGamma" )
	{
		ReturnableFitResult = PetesGamma_DoSafeFit( MinimiserConfig, FunctionConfig, internalBottleParameters, BottleData, BottleConstraints, forceContinue, OutputLevel );
	}
	else if( FunctionConfig->GetStrategy() == "Robs" )
	{
		ReturnableFitResult = Robs_DoSafeFit( MinimiserConfig, FunctionConfig, internalBottleParameters, BottleData, BottleConstraints, forceContinue, OutputLevel );
	}
	else  //	Default!
	{
		ReturnableFitResult = DoSingleSafeFit( MinimiserConfig, FunctionConfig, internalBottleParameters, BottleData, BottleConstraints, forceContinue, OutputLevel );
	}

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler:: Finished Passing Back, checking Ourput FitResult contains all input Parameters" << endl;
	}

	/*!	Check that the ParameterSet in the output contains the provided ParameterSet regardless of what was passed to the Minimiser	*/
	FitAssembler::CheckParameterSet( ReturnableFitResult, internalBottleParameters );

	BottleParameters->AddPhysicsParameters( internalBottleParameters, true );

	delete internalBottleParameters;

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler:: Returning FitResult" << endl;
	}

	return ReturnableFitResult;
}

FitResult * FitAssembler::Robs_DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, const ParameterSet* BottleParameters,
		const vector< PDFWithData* > BottleData, const vector< ConstraintFunction* > BottleConstraints, bool forceContinue, const int OutputLevel )
{
	cout << "Starting Fit1:" << endl;
	double deltaPara = BottleParameters->GetPhysicsParameter( string("delta_para") )->GetBlindedValue();
	double deltaPerp = BottleParameters->GetPhysicsParameter( string("delta_perp") )->GetBlindedValue();
	double deltaGamma = BottleParameters->GetPhysicsParameter( string("deltaGamma") )->GetBlindedValue();
	double phi_s = BottleParameters->GetPhysicsParameter( string("Phi_s") )->GetBlindedValue();

	ParameterSet* set1 = new ParameterSet( *BottleParameters );

	// Normal fit
	FitResult* res0 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, set1, BottleData, BottleConstraints, forceContinue, OutputLevel );
	bool good_result_0 = res0->GetFitStatus() == 3;
	cout << "Finished Fit1." << endl;
	if( !good_result_0  ) cout << "Fit-1 failed" << endl;

	ParameterSet* set2 = new ParameterSet( *BottleParameters );

	// Conjugate fit
	//double deltaS = BottleParameters->GetPhysicsParameter( string("delta_s") )->GetBlindedValue();
	if( set2->GetPhysicsParameter( string("delta_para") )->GetType() != "Fixed" ) set2->GetPhysicsParameter( string("delta_para") )->SetBlindedValue( -deltaPara );
	if( set2->GetPhysicsParameter( string("delta_perp") )->GetType() != "Fixed" ) set2->GetPhysicsParameter( string("delta_perp") )->SetBlindedValue( 3.14159-deltaPerp );
	//BottleParameters->GetPhysicsParameter( string("delta_s") )->SetBlindedValue( -deltaS );
	cout << endl << "Starting Fit2:" << endl;
	FitResult* res1 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, set2, BottleData, BottleConstraints, forceContinue, OutputLevel );
	cout << "Finished Fit2." << endl;
	bool good_result_1 = res1->GetFitStatus() == 3;

	ParameterSet* set3 = new ParameterSet( *BottleParameters );

	if( set3->GetPhysicsParameter( string("delta_para") )->GetType() != "Fixed" ) set3->GetPhysicsParameter( string("delta_para") )->SetBlindedValue( deltaPara );
	if( set3->GetPhysicsParameter( string("delta_perp") )->GetType() != "Fixed" ) set3->GetPhysicsParameter( string("delta_perp") )->SetBlindedValue( deltaPerp );
	if( set3->GetPhysicsParameter( string("deltaGamma") )->GetType() != "Fixed" ) set3->GetPhysicsParameter( string("deltaGamma") )->SetBlindedValue( -deltaGamma );
	if( set3->GetPhysicsParameter( string("Phi_s") )->GetType() != "Fixed" ) set3->GetPhysicsParameter( string("Phi_s") )->SetBlindedValue( 3.14159-phi_s );
	cout << endl << "Starting Fit3:" << endl;
	FitResult* res2 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, set3, BottleData, BottleConstraints, forceContinue, OutputLevel );
	cout << "Finished Fit3." << endl;
	bool good_result_2 = res2->GetFitStatus() == 3;

	ParameterSet* set4 = new ParameterSet( *BottleParameters );

	if( set4->GetPhysicsParameter( string("delta_para") )->GetType() != "Fixed" ) set4->GetPhysicsParameter( string("delta_para") )->SetBlindedValue( -deltaPara );
	if( set4->GetPhysicsParameter( string("delta_perp") )->GetType() != "Fixed" ) set4->GetPhysicsParameter( string("delta_perp") )->SetBlindedValue( 3.14159-deltaPerp );
	if( set4->GetPhysicsParameter( string("deltaGamma") )->GetType() != "Fixed" ) set4->GetPhysicsParameter( string("deltaGamma") )->SetBlindedValue( -deltaGamma );
	if( set4->GetPhysicsParameter( string("Phi_s") )->GetType() != "Fixed" ) set4->GetPhysicsParameter( string("Phi_s") )->SetBlindedValue( 3.14159-phi_s );
	cout << endl << "Starting Fit4:" << endl;
	FitResult* res3 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, set4, BottleData, BottleConstraints, forceContinue, OutputLevel );
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
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, bool forceContinue, int OutputLevel )
{
	streambuf *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;

	cout_bak = cout.rdbuf();
	cerr_bak = cerr.rdbuf();
	clog_bak = clog.rdbuf();

	if( DebugClass::DebugThisClass( "FitAssembler" ) )
	{
		cout << "FitAssembler: Debugging Start" << endl;
	}

	//      If the user wanted silence we point the Std Output Streams to the oblivion of NULL
	if( OutputLevel <= -1 )
	{
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
		if( ReturnableFitResult->GetFitStatus() < 2 )
		{
			bad_fit = true;
		}
	} else {
		bad_fit = true;
	}

	if( !forceContinue && bad_fit )
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
		delete Bad_Bottle;
	}

	return ReturnableFitResult;
}


FitResult * FitAssembler::Petes_DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration* FunctionConfig, ParameterSet* BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, bool forceContinue, int OutputLevel )
{
	cout << endl << "******* Result of Petes Double fit strategy  (NOT flipping delta_s as well) *********" << endl ;
	cout << "Starting Fit1:" << endl;
	// Normal fit
	FitResult* res0 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, forceContinue, OutputLevel ) ;
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
	FitResult* res1 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, forceContinue, OutputLevel );
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

FitResult * FitAssembler::PetesGamma_DoSafeFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration* FunctionConfig, ParameterSet* BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints, bool forceContinue, int OutputLevel )
{
	cout << endl << "******* Result of Petes Double fit strategy for flipping deltaGamma *********" << endl ;
	cout << "Starting Fit1:" << endl;
	// Normal fit
	FitResult* res0 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, forceContinue, OutputLevel );
	double LLmin0 = res0->GetMinimumValue() ;
	bool good_result_0 = res0->GetFitStatus() == 3;
	cout << "Finished Fit1." << endl;
	if( !good_result_0  ) cout << "Fit-1 failed" << endl;

	// Conjugate fit
	double deltaGamma = BottleParameters->GetPhysicsParameter( string("deltaGamma") )->GetBlindedValue() ;
	BottleParameters->GetPhysicsParameter( string("deltaGamma") )->SetBlindedValue( -deltaGamma ) ;
	cout << endl << "Starting Fit2:" << endl;
	FitResult* res1 = DoSingleSafeFit( MinimiserConfig, FunctionConfig, BottleParameters, BottleData, BottleConstraints, forceContinue, OutputLevel );
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

