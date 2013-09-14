/**
  @class FitFunction

  Parent class for the function to minimise
  Overload the evaluate methods and UP value for Chi2, NLL, etc.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	ROOT Headers
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#ifdef RAPIDFIT_USETGLTIMER
#include "TGLStopwatch.h"
#endif
//	RapidFit Headers
#include "FitFunction.h"
#include "Threading.h"
#include "ClassLookUp.h"
#include "RapidFitIntegrator.h"
#include "StringProcessing.h"
#include "MemoryDataSet.h"
#include "ProdPDF.h"
//	System Headers
#include <iostream>
#include <iomanip>
#include <float.h>
#include <cstdlib>
#include <ctime>

using namespace::std;

//Default constructor
FitFunction::FitFunction() :
	Name("Unknown"), allData(), testDouble(), useWeights(false), weightObservableName(), Fit_File(NULL), Fit_Tree(NULL), branch_objects(), branch_names(), fit_calls(0),
	Threads(-1), stored_pdfs(), StoredBoundary(), StoredDataSubSet(), StoredIntegrals(), finalised(false), fit_thread_data(NULL), testIntegrator( true ), weightsSquared( false ),
	debug(new DebugClass(false) ), traceNum(0), step_time(-1), callNum(0), integrationConfig(new RapidFitIntegratorConfig())
{
}

//Destructor
FitFunction::~FitFunction()
{
	//cout << "Hello from FitFunction destructor" << endl;

	//	Close any open files...
	//	common sence and OO says call the destructors too... ROOT says not to and I'm too fed up to argue!
	if( Fit_File != NULL )
	{
		Fit_Tree->Write();
		Fit_File->Close();
	}
	if( fit_thread_data != NULL ) delete [] fit_thread_data;
	//if( allData != NULL ) delete allData;
	while( !StoredBoundary.empty() )
	{
		if( StoredBoundary.back() != NULL ) delete StoredBoundary.back();
		StoredBoundary.pop_back();
	}
	while( !stored_pdfs.empty() )
	{
		if( stored_pdfs.back() != NULL ) delete stored_pdfs.back();
		stored_pdfs.pop_back();
	}
	while( !StoredIntegrals.empty() )
	{
		if( StoredIntegrals.back() != NULL ) delete StoredIntegrals.back();
		StoredIntegrals.pop_back();
	}

	if( debug != NULL ) delete debug;
	if( integrationConfig != NULL ) delete integrationConfig;
}

void FitFunction::SetupTrace( const TString FileName, const int inputTraceNum )
{
	//	Create the output file
	Fit_File = new TFile( FileName, "UPDATE" );
	traceNum = inputTraceNum;
}

void FitFunction::SetIntegratorConfig( const RapidFitIntegratorConfig* gsl )
{
	if( integrationConfig != NULL ) delete integrationConfig;
	integrationConfig = new RapidFitIntegratorConfig( *gsl );
}

void FitFunction::SetupTraceTree()
{
	Fit_File->cd();
	TString TraceName("Trace_");
	TraceName+=traceNum;
	Fit_Tree = new TTree( TraceName, TraceName );

	//	Yes I could point the FitFunction to the address of the objects in memory in RapidFit...
	//	However that seems INCREADIBLY DANGEROUS

	//	Initialize the branch structure within the TTree
	for( unsigned int i=0; i< allData->GetParameterSet()->GetAllNames().size(); ++i )
	{
		branch_objects.push_back( 0 );
		branch_names.push_back( allData->GetParameterSet()->GetAllNames()[i] );
		TString Branch_Name = string( branch_names.back() ).c_str();
		TString Branch_Name_2( string( branch_names.back() ).c_str()); Branch_Name_2.Append("/D");
		Fit_Tree->Branch( Branch_Name, &(branch_objects.back()), Branch_Name_2 );
	}

	branch_objects.push_back( 0 );
	Fit_Tree->Branch( "NLL", &(branch_objects.back()), "NLL/D" );
	Fit_Tree->Branch( "Call", &(fit_calls), "Call/D" );
	Fit_Tree->Branch( "time", &(step_time), "time/D" );
}

//Set the physics bottle to fit with
void FitFunction::SetPhysicsBottle( const PhysicsBottle * NewBottle )
{
	allData = new PhysicsBottle( *NewBottle );
	if( Fit_File != NULL ) this->SetupTraceTree();

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "FitFunction" ) )
		{
			cout << "FitFunction: I am Performing a fit with " << NewBottle->NumberResults() << " seperate NLLs" << endl;
			cout << "FitFunction: I am using the " << Name << " Fit Function to Evaluate" << endl;
			cout << "FitFunction: I have been asked to use " << Threads << " parallel threads" << endl;
		}
	}

	//Initialise the integrators
	for( int resultIndex = 0; resultIndex < NewBottle->NumberResults(); ++resultIndex )
	{
		//	Update Internal ParameterSet in PDF
		NewBottle->GetResultPDF(resultIndex)->UpdatePhysicsParameters( allData->GetParameterSet() );

		if( debug != NULL )
		{
			if( debug->DebugThisClass( "FitFunction" ) )
			{
				cout << "FitFunction: Constructing Integrator Object for ToFit " << resultIndex+1 << endl;
			}
		}

		//RapidFitIntegrator * resultIntegrator =  new RapidFitIntegrator( NewBottle->GetResultPDF(resultIndex), false, gslIntegrator );
		//resultIntegrator->SetDebug( debug );

		if( debug != NULL )
		{
			if( debug->DebugThisClass( "FitFunction" ) )
			{
				if( testIntegrator )
				{
					cout << "FitFunction: Performing Integrator test" << endl;
				}
				else
				{
					cout << "FitFunction: NOT Performing Integrator test" << endl;
				}
			}
		}

		if( debug != NULL )
		{
			if( debug->DebugThisClass( "FitFunction" ) )
			{
				cout << "FitFunction: Performing Integration Test" << endl;
			}
		}

		double someVal=0.;
		NewBottle->GetResultPDF(resultIndex)->GetPDFIntegrator()->ForceTestStatus( false );
		NewBottle->GetResultPDF(resultIndex)->SetUpIntegrator( integrationConfig );
		allData->GetResultPDF(resultIndex)->SetUpIntegrator( integrationConfig );

		if( integrationConfig->useGSLIntegrator ) cout << "Using GSL!" << endl;

		if( NewBottle->GetResultDataSet(resultIndex)->GetDataNumber() > 0 )
		{
			if( testIntegrator )
			{
				allData->GetResultPDF(resultIndex)->GetPDFIntegrator()->ForceTestStatus( false );
				allData->GetResultPDF(resultIndex)->GetPDFIntegrator()->SetDebug( debug );
				someVal = allData->GetResultPDF(resultIndex)->GetPDFIntegrator()->Integral(
						allData->GetResultDataSet(resultIndex)->GetDataPoint( 0 ),
						allData->GetResultDataSet(resultIndex)->GetBoundary() );
			}
		}
		(void) someVal;

		if( debug != NULL )
		{
			if( debug->DebugThisClass( "FitFunction" ) )
			{
				cout << "FitFunction: Finished Performing Integration Test" << endl;
			}
		}

		if( Threads > 0 )
		{
			//      Create simple data subsets. We no longer care about the handles that IDataSet takes care of
			if( debug != NULL )
			{
				if( debug->DebugThisClass( "FitFunction" ) )
				{
					cout << "FitFunction: Splitting DataSet" << endl;
				}
			}
			StoredDataSubSet.push_back( Threading::divideData( NewBottle->GetResultDataSet(resultIndex), Threads ) );
			vector<IDataSet*> sets;
			for( unsigned int i=0; i<  StoredDataSubSet.back().size(); ++i )
			{
				sets.push_back( new MemoryDataSet( NewBottle->GetResultDataSet(resultIndex)->GetBoundary(), StoredDataSubSet.back()[i] ) );
			}
			stored_datasets.push_back( sets );
			for( int i=0; i< Threads; ++i )
			{
				if( debug != NULL )
				{
					if( debug->DebugThisClass( "FitFunction" ) )
					{
						cout << "FitFunction: Cloning PhaseSpaceBoundary" << endl;
					}
				}
				StoredBoundary.push_back( new PhaseSpaceBoundary( *(NewBottle->GetResultDataSet(resultIndex)->GetBoundary()) ) );
				if( debug != NULL )
				{
					if( debug->DebugThisClass( "FitFunction" ) )
					{
						cout << "FitFunction: CopyingPdf " << NewBottle->GetResultPDF( resultIndex )->GetLabel() << endl;
					}
				}
				stored_pdfs.push_back( ClassLookUp::CopyPDF( NewBottle->GetResultPDF( resultIndex ) ) );
				stored_pdfs.back()->SetDebug( debug );
				stored_pdfs.back()->SetUpIntegrator( integrationConfig );
			}
		}
	}

	if( Threads > 0 )
	{
		fit_thread_data = new Fitting_Thread[ (unsigned) Threads ];
	}

	if( debug != NULL )
	{
		if( debug->DebugThisClass( "FitFunction" ) )
		{
			cout << "FitFunction: PhysicsBottle Set" << endl;
		}
	}
}

//Return the physics bottle
PhysicsBottle* FitFunction::GetPhysicsBottle() const
{
	return allData;
}

// Get and set the fit parameters
void FitFunction::SetParameterSet( const ParameterSet * NewParameters )
{
	allData->SetParameterSet(NewParameters);

	//Initialise the integrators
	for ( int resultIndex = 0; resultIndex < allData->NumberResults(); ++resultIndex )
	{
		allData->SetParameterSet( NewParameters );

		allData->GetResultPDF( resultIndex )->UpdatePhysicsParameters( allData->GetParameterSet() );

		for( int i=0; i< Threads; ++i )
		{
			stored_pdfs[ (unsigned)(i + resultIndex*Threads) ]->UpdatePhysicsParameters( allData->GetParameterSet() );
			//stored_pdfs[ (unsigned)(i + resultIndex*Threads) ]->UnsetCache();
		}
	}
}

ParameterSet * FitFunction::GetParameterSet() const
{
	return allData->GetParameterSet();
}

//Return the value to minimise
double FitFunction::Evaluate()
{
	++callNum;
	//time_t start, end;
	//time(&start);

#ifdef RAPIDFIT_USETGLTIMER
	TGLStopwatch* thisWatch = NULL;
	if( Fit_Tree !=NULL )
	{
		thisWatch = new TGLStopwatch();
		thisWatch->Start();
	}
#endif
	double minimiseValue = 0.0;
	double temp=0.;
	//Calculate the function value for each PDF-DataSet pair
	for( int resultIndex = 0; resultIndex < allData->NumberResults(); ++resultIndex )
	{
		//cout << endl << resultIndex << ": " << allData->GetResultDataSet( resultIndex )->GetDataNumber() << endl;
		if( allData->GetResultDataSet( resultIndex )->GetDataNumber() == 0 )
		{
			if( callNum < 6 )
			{
				cerr << "Are you aware DataSet: " << resultIndex+1 << " has zero size?" << endl;
			}
			continue;
		}
		if( allData->GetResultDataSet( resultIndex )->GetDataNumber() < 1 )
		{
			temp = 0.;
		}
		else
		{
			//cout << "Eval Set: " << allData->GetResultDataSet( resultIndex ) << "\t" << resultIndex << endl;
			temp = this->EvaluateDataSet( allData->GetResultPDF( resultIndex ), allData->GetResultDataSet( resultIndex ), resultIndex );
			//cout << "Result: " << temp << endl;
		}
		if( abs(temp) >= DBL_MAX )
		{
			return DBL_MAX;
		}
		else
		{
			minimiseValue+=temp;
		}
		//cout << "temp: " << minimiseValue << endl;
	}

	//Calculate the value of each constraint
	vector< ConstraintFunction* > constraints = allData->GetConstraints();
	for (unsigned int constraintIndex = 0; constraintIndex < constraints.size(); ++constraintIndex )
	{
		if( minimiseValue < DBL_MAX )
		{
			minimiseValue += constraints[constraintIndex]->Evaluate( allData->GetParameterSet() );
		}
		else
		{
			return DBL_MAX;
		}
	}

	//time(&end);
	++fit_calls;
	//step_time = difftime( end, start );
	if( Fit_Tree !=NULL )
	{
#ifdef RAPIDFIT_USETGLTIMER
		step_time = thisWatch->End();
#else
		step_time = -1.;
#endif
		//thisWatch->Stop();
		//step_time = thisWatch->CpuTime();
		//delete thisWatch;

		for(unsigned int i=0; i< allData->GetParameterSet()->GetAllNames().size(); ++i )
		{
			branch_objects[i] = (Double_t) allData->GetParameterSet()->GetPhysicsParameter( branch_names[i] )->GetBlindedValue();
			//cout << (Double_t) branch_objects[i] << "\t" ;
			Fit_Tree->SetBranchAddress( string(branch_names[i]).c_str(), &(branch_objects[i]) );
		}
		branch_objects[branch_objects.size()] = (Double_t) minimiseValue;
		Fit_Tree->SetBranchAddress( "NLL", &(branch_objects[branch_objects.size()]) );
		Fit_Tree->SetBranchAddress( "Call", &(fit_calls) );

		Fit_Tree->SetBranchAddress( "time", &(step_time) );
		//cout << endl;
		Fit_Tree->Fill();
		Fit_Tree->Write("",TObject::kOverwrite);
		Fit_File->Write("",TObject::kOverwrite);
	}


	if( std::isnan(minimiseValue) )
	{
		this->GetParameterSet()->Print();
		minimiseValue = DBL_MAX;
	}
	if( debug != NULL )
	{
		if( debug->DebugThisClass( "FitFunction" ) )
		{
			cout << endl;
		}
	}

	return minimiseValue;
}

//Return the value to minimise for a given PDF/DataSet pair
double FitFunction::EvaluateDataSet( IPDF * TestPDF, IDataSet * TestDataSet, int number )
{
	(void)TestPDF;
	(void)TestDataSet;
	(void)number;

	return 1.0;
}

//Return the Up value for error calculation
double FitFunction::UpErrorValue( const int Sigma )
{
	(void)Sigma;
	return 1.0;
}

//Set the FitFunction to use per-event weights
void FitFunction::UseEventWeights( const string WeightName )
{
	useWeights = true;
	weightObservableName = WeightName;
}

void FitFunction::SetThreads( const int input )
{
	Threads = input;
	//      Get the number of cores on the compile machine
	unsigned int num_cores = (unsigned)Threading::numCores();

	if( input < 0 )
	{
		Threads = (int)num_cores;
	}
}

int FitFunction::GetThreads() const
{
	return Threads;
}

void FitFunction::SetIntegratorTest( const bool input )
{
	testIntegrator = input;
}

void FitFunction::SetUseWeightsSquared( const bool Input )
{
	weightsSquared = Input;
}

string FitFunction::GetWeightName() const
{
	return weightObservableName;
}

bool FitFunction::GetWeightsWereUsed() const
{
	return useWeights;
}

vector<string> FitFunction::ConstrainedParameter() const
{
	vector<string> allparams;
	vector<ConstraintFunction*> allconstraints = allData->GetConstraints();
	for( unsigned int i=0; i< allconstraints.size(); ++i )
	{
		allparams = StringProcessing::CombineUniques( allparams, allconstraints[i]->ConstrainedParameter() );
	}
	return allparams;
}

void FitFunction::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
}

unsigned int FitFunction::GetCallNum()
{
	return callNum;
}

