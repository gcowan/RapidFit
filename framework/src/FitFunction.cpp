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
//	System Headers
#include <iostream>
#include <iomanip>

//Default constructor
FitFunction::FitFunction() : allData(), allIntegrators(), testDouble(), useWeights(false), weightObservableName(), Fit_File(NULL), Fit_Tree(NULL), branch_objects(), branch_names(), fit_calls(0), Threads(0)
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
}

void FitFunction::SetupTrace( TString FileName, int traceNum )
{
	//	Create the output file
	Fit_File = new TFile( FileName, "UPDATE" );
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
void FitFunction::SetPhysicsBottle( PhysicsBottle * NewBottle )
{
	NewBottle->Finalise();
	allData = NewBottle;

	//Initialise the integrators
	for ( int resultIndex = 0; resultIndex < NewBottle->NumberResults(); ++resultIndex )
	{
		RapidFitIntegrator * resultIntegrator = new RapidFitIntegrator( NewBottle->GetResultPDF(resultIndex) );
		allIntegrators.push_back( resultIntegrator );
	}
}

//Return the physics bottle
PhysicsBottle* FitFunction::GetPhysicsBottle()
{
	return allData;
}

// Get and set the fit parameters
bool FitFunction::SetParameterSet( ParameterSet * NewParameters )
{
	return allData->SetParameterSet(NewParameters);
}
ParameterSet * FitFunction::GetParameterSet()
{
	return allData->GetParameterSet();
}

//Return the value to minimise
double FitFunction::Evaluate()
{
	double minimiseValue = 0.0;

	//Calculate the function value for each PDF-DataSet pair
	for (int resultIndex = 0; resultIndex < allData->NumberResults(); ++resultIndex)
	{
		minimiseValue += EvaluateDataSet( allData->GetResultPDF( resultIndex ), allData->GetResultDataSet( resultIndex ), allIntegrators[unsigned(resultIndex)] );
	}

	//Calculate the value of each constraint
	vector< ConstraintFunction* > constraints = allData->GetConstraints();
	for (unsigned int constraintIndex = 0; constraintIndex < constraints.size(); ++constraintIndex )
	{
		minimiseValue += constraints[constraintIndex]->Evaluate( allData->GetParameterSet() );
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
	cout << "NLL: " << setprecision(20) << minimiseValue << "\r\r\r" << flush;
	return minimiseValue;
}

//Return the value to minimise for a given PDF/DataSet pair
double FitFunction::EvaluateDataSet( IPDF * TestPDF, IDataSet * TestDataSet, RapidFitIntegrator * ResultIntegrator )
{
	//	Stupid gcc
	(void)TestPDF;
	(void)TestDataSet;
	(void)ResultIntegrator;

	return 1.0;
}

//Return the Up value for error calculation
double FitFunction::UpErrorValue( int Sigma )
{
	//	Stupid gcc
	(void)Sigma;
	return 1.0;
}

//Finalise the PhysicsBottle
void FitFunction::Finalise()
{
	allData->Finalise();
}

//Set the FitFunction to use per-event weights
void FitFunction::UseEventWeights( string WeightName )
{
	useWeights = true;
	weightObservableName = WeightName;
}

void FitFunction::SetThreads( int input )
{
	Threads = input;
}

int FitFunction::GetThreads()
{
	return Threads;
}
