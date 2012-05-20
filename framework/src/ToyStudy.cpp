/**
  @class ToyStudy

  A class that automates a whole toy study

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	RapidFit Headers
#include "ToyStudy.h"
#include "FitAssembler.h"
#include "XMLConfigReader.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>

using namespace::std;

//Constructor using an XML config file directly
ToyStudy::ToyStudy( string FileName ) : IStudy()
{
	XMLConfigReader * xml = new XMLConfigReader( FileName );

	theMinimiser = xml->GetMinimiserConfiguration();
	theFunction = xml->GetFitFunctionConfiguration();
	theFunction->SetIntegratorTest( false );
	studyParameters = xml->GetFitParameters();
	pdfsAndData = xml->GetPDFsAndData();
	numberStudies = xml->GetNumberRepeats();
	allConstraints = xml->GetConstraints();
	delete_objects = true;

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

//Constructor with correct arguments
ToyStudy::ToyStudy( MinimiserConfiguration * TheMinimiser, FitFunctionConfiguration * TheFunction, ParameterSet* StudyParameters, vector< PDFWithData* > PDFsAndData, vector< ConstraintFunction* > InputConstraints, int NumberStudies ) : IStudy()
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

	//Loop over all studies
	for ( int studyIndex = 0; studyIndex < numberStudies; ++studyIndex )
	{
		cout << "\n\n\t\tStarting ToyStudy\t\t" << studyIndex+1 << "\tof:\t" << numberStudies << endl;
		allResults->StartStopwatch();

		FitResult* new_result = FitAssembler::DoSafeFit( theMinimiser, theFunction, studyParameters, pdfsAndData, allConstraints, OutputLevel );

		//cout << "\n\nFinished Study" << endl;

		//	Have to explicitly call this to request the data be deleted between runs, not normally an issue but causes problems with large (pb) datasets
		for( unsigned int i=0; i< pdfsAndData.size(); ++i )
		{
			//cout << "Removing DataSet Number: " << i << endl;
			pdfsAndData[i]->ClearCache();
		}

		if( new_result->GetFitStatus() != 3 )
		{
			cerr << "Fit fell over!\t Requesting another fit." << endl;
			++numberStudies;
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

