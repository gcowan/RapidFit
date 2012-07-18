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
//	RapidFit Headers
#include "FitFunction.h"
#include "Threading.h"
#include "ClassLookUp.h"
#include "RapidFitIntegrator.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <iomanip>
#include <float.h>
#include <cstdlib>

using namespace::std;

//Default constructor
FitFunction::FitFunction() :
	allData(), allIntegrators(), testDouble(), useWeights(false), weightObservableName(), Fit_File(NULL), Fit_Tree(NULL), branch_objects(), branch_names(), fit_calls(0),
	Threads(-1), stored_pdfs(), StoredBoundary(), StoredDataSubSet(), StoredIntegrals(), finalised(false), fit_thread_data(NULL), testIntegrator( true ), weightsSquared( false ),
	debug(new DebugClass(false) )
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
	while( !allIntegrators.empty() )
	{
		if( allIntegrators.back() != NULL ) delete allIntegrators.back();
		allIntegrators.pop_back();
	}

	if( debug != NULL ) delete debug;
}

void FitFunction::SetupTrace( const TString FileName, const int inputTraceNum )
{
	//	Create the output file
	Fit_File = new TFile( FileName, "UPDATE" );
	traceNum = inputTraceNum;
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
}

//Set the physics bottle to fit with
void FitFunction::SetPhysicsBottle( const PhysicsBottle * NewBottle )
{
	allData = new PhysicsBottle( *NewBottle );

	if( Fit_File != NULL ) this->SetupTraceTree();

	//Initialise the integrators
	for ( int resultIndex = 0; resultIndex < NewBottle->NumberResults(); ++resultIndex )
	{
		NewBottle->GetResultPDF(resultIndex)->UpdatePhysicsParameters( allData->GetParameterSet() );

		RapidFitIntegrator * resultIntegrator =  new RapidFitIntegrator( NewBottle->GetResultPDF(resultIndex) );
		resultIntegrator->SetDebug( debug );
		allIntegrators.push_back( resultIntegrator );

		allIntegrators.back()->ForceTestStatus( false );
		double someval=0.;
		if( testIntegrator == false ) allIntegrators.back()->ForceTestStatus( true );
		else someval = allIntegrators.back()->Integral( NewBottle->GetResultDataSet(resultIndex)->GetDataPoint(0), NewBottle->GetResultDataSet(resultIndex)->GetBoundary() );
		(void) someval;
		allIntegrators.back()->ForceTestStatus( true );

		if( Threads > 0 )
		{
			//      Create simple data subsets. We no longer care about the handles that IDataSet takes care of
			StoredDataSubSet.push_back( Threading::divideData( NewBottle->GetResultDataSet(resultIndex), Threads ) );
			for( int i=0; i< Threads; ++i )
			{
				StoredBoundary.push_back( new PhaseSpaceBoundary( *(NewBottle->GetResultDataSet(resultIndex)->GetBoundary()) ) );
				stored_pdfs.push_back( ClassLookUp::CopyPDF( NewBottle->GetResultPDF( resultIndex ) ) );
				stored_pdfs.back()->SetDebug( debug );
				StoredIntegrals.push_back( new RapidFitIntegrator( stored_pdfs.back() ) );
				StoredIntegrals.back()->SetDebug( debug );

				// Give the Integral Testing code a Sharp Kick!!!
				StoredIntegrals.back()->ForceTestStatus( true );
			}
		}
	}

	if( Threads > 0 )
	{
		fit_thread_data = new Fitting_Thread[ (unsigned) Threads ];
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
			stored_pdfs[ (unsigned)(i + resultIndex*Threads) ]->UnsetCache();
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
	double minimiseValue = 0.0;
	double temp=0.;
	//Calculate the function value for each PDF-DataSet pair
	for( int resultIndex = 0; resultIndex < allData->NumberResults(); ++resultIndex )
	{
		if( allData->GetResultDataSet( resultIndex )->GetBoundary() == 0 )
		{
			cerr << "Are you aware DataSet: " << resultIndex+1 << " has zero size?" << endl;
			continue;
		}
		temp = this->EvaluateDataSet( allData->GetResultPDF( resultIndex ), allData->GetResultDataSet( resultIndex ), allIntegrators[unsigned(resultIndex)], resultIndex );
		if( temp >= DBL_MAX )
		{
			minimiseValue=DBL_MAX;
			break;
		}
		else
		{
			minimiseValue+=temp;
		}
	}

	//Calculate the value of each constraint
	vector< ConstraintFunction* > constraints = allData->GetConstraints();
	for (unsigned int constraintIndex = 0; constraintIndex < constraints.size(); ++constraintIndex )
	{
		if( minimiseValue < DBL_MAX )
		{
			minimiseValue += constraints[constraintIndex]->Evaluate( allData->GetParameterSet() );
		}
	}

	++fit_calls;

	if( Fit_Tree !=NULL )
	{
		for(unsigned int i=0; i< allData->GetParameterSet()->GetAllNames().size(); ++i )
		{
			branch_objects[i] = (Double_t) allData->GetParameterSet()->GetPhysicsParameter( branch_names[i] )->GetBlindedValue();
			//cout << (Double_t) branch_objects[i] << "\t" ;
			Fit_Tree->SetBranchAddress( string(branch_names[i]).c_str(), &(branch_objects[i]) );
		}
		branch_objects[branch_objects.size()] = (Double_t) minimiseValue;
		Fit_Tree->SetBranchAddress( "NLL", &(branch_objects[branch_objects.size()]) );
		Fit_Tree->SetBranchAddress( "Call", &(fit_calls) );
		//cout << endl;
		Fit_Tree->Fill();
	}
	cout << "NLL: " << setprecision(10) << minimiseValue << "\b\b\b\b\r\r\r\r" << flush;
	if( isnan(minimiseValue) )
	{
		minimiseValue = DBL_MAX;
	}
	//cout << endl;
	return minimiseValue;
}

//Return the value to minimise for a given PDF/DataSet pair
double FitFunction::EvaluateDataSet( IPDF * TestPDF, IDataSet * TestDataSet, RapidFitIntegrator * ResultIntegrator, int number )
{
	(void)TestPDF;
	(void)TestDataSet;
	(void)ResultIntegrator;
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
	if( input_debug != NULL )
	{
		if( debug != NULL ) delete debug;
		debug = new DebugClass( *input_debug );

		if( debug->DebugThisClass("FitFunction") )
		{
			debug->SetStatus(true);
			cout << "FitFunction: Debugging Enabled!" << endl;
		}
		else
		{
			debug->SetStatus(false);
		}
	}
	else
	{
		if( debug != NULL ) delete debug;
		debug = new DebugClass(false);
	}
}

