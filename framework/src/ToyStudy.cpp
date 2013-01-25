/**
  @class ToyStudy

  A class that automates a whole toy study

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

//	ROOT Headers
#include "TString.h"
//	RapidFit Headers
#include "ToyStudy.h"
#include "FitAssembler.h"
#include "XMLConfigReader.h"
#include "StringProcessing.h"
#include "ResultFormatter.h"
//	System Headers
#include <iostream>

using namespace::std;

//Constructor with correct arguments
ToyStudy::ToyStudy( MinimiserConfiguration * TheMinimiser, FitFunctionConfiguration * TheFunction, ParameterSet* StudyParameters,
		vector< PDFWithData* > PDFsAndData, vector< ConstraintFunction* > InputConstraints, int NumberStudies ) :
		IStudy(), fixedNumToys(false), saveAllToys(false)
{
	pdfsAndData = PDFsAndData;
	studyParameters = StudyParameters;
	theMinimiser = TheMinimiser;
	theFunction = TheFunction;
	theFunction->SetIntegratorTest( false );
	allResults = NULL;
	numberStudies = NumberStudies;
	allConstraints = InputConstraints;

	if ( numberStudies < 1 )
	{
		cerr << "Erroneous number of studies: defaulting to 1" << endl;
		numberStudies = 1;
	}
	else if ( numberStudies == 1 )
	{
		cout << "Doing a single toy study. Check this is what you expect." << endl;
	}
}

//Destructor
ToyStudy::~ToyStudy()
{
	if( allResults!=NULL ) delete allResults;
}

void ToyStudy::SetFixedNumberToys()
{
	fixedNumToys = true;
}

void ToyStudy::setSaveAllToys()
{
	saveAllToys = true;
}

//Automate the toy study
void ToyStudy::DoWholeStudy( int OutputLevel )
{
	//Make a vector of unique parameter names
	vector<string> uniqueNames;
	for (unsigned int pdfIndex = 0; pdfIndex < pdfsAndData.size(); ++pdfIndex )
	{
		//This is not strictly necessary, but suppresses a warning message
		pdfsAndData[pdfIndex]->SetPhysicsParameters(studyParameters);

		uniqueNames = StringProcessing::CombineUniques( uniqueNames, pdfsAndData[pdfIndex]->GetPDF()->GetPrototypeParameterSet() );
	}
	allResults = new FitResultVector(uniqueNames);

	//	Do NOT want PDFWithData to cache the dataset as this will just keep running over the same data
	for( unsigned int i=0; i< pdfsAndData.size(); ++i )
	{
		pdfsAndData[i]->SetUseCache( false );
	}

	TString filename="filename_";

	//Loop over all studies
	for ( int studyIndex = 0; studyIndex < numberStudies; ++studyIndex )
	{
		cout << "\n\n\t\tStarting ToyStudy\t\t" << studyIndex+1 << "\tof:\t" << numberStudies << endl;
		allResults->StartStopwatch();

		ParameterSet* thisSet = new ParameterSet( *studyParameters );

		FitResult* new_result = FitAssembler::DoSafeFit( theMinimiser, theFunction, thisSet, pdfsAndData, allConstraints, false, OutputLevel );

		delete thisSet;

		//cout << "\n\nFinished Study" << endl;

		TString this_filename = filename;
		this_filename.Append("_S");
		this_filename+=studyIndex; 

		//	Have to explicitly call this to request the data be deleted between runs, not normally an issue but causes problems with large (pb) datasets
		for( unsigned int i=0; i< pdfsAndData.size(); ++i )
		{
			TString this_filename2=this_filename;
			this_filename2.Append("_D");
			this_filename2+=i;
			this_filename2.Append(".root");
			if( saveAllToys ) ResultFormatter::MakeRootDataFile( this_filename2.Data(), vector<IDataSet*>(1, pdfsAndData[i]->GetDataSet()) );
			//cout << "Removing DataSet Number: " << i << endl;
			pdfsAndData[i]->ClearCache();
		}

		if( new_result->GetFitStatus() != 3 )
		{
			cerr << "Fit fell over!\t Requesting another fit." << endl;
			if( !fixedNumToys ) ++numberStudies;
		}

		allResults->AddFitResult( new_result );
	}
}

//Get the result of the toy study
FitResultVector* ToyStudy::GetStudyResult()
{
	return allResults;
}


//	Set the number of repeats
void ToyStudy::SetNumRepeats( int new_num_repeats )
{
	numberStudies = new_num_repeats;
}

//	Set the command line physics parameters
//	Not actually implemented in the ToyStudy class ATM
void ToyStudy::SetCommandLineParams( vector<string> new_CommandLineParams )
{
	(void) new_CommandLineParams;
}

