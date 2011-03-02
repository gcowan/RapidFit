/**
  @class ToyStudy

  A class that automates a whole toy study

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

#include "ToyStudy.h"
#include "FitAssembler.h"
#include "XMLConfigReader.h"
#include "StringProcessing.h"
#include <iostream>

using namespace std;

//Default constructor
ToyStudy::ToyStudy()
{
}

//Constructor using an XML config file directly
ToyStudy::ToyStudy( string FileName )
{
	XMLConfigReader * xml = new XMLConfigReader(FileName);
	if ( xml->IsLoaded() )
	{
		theMinimiser = xml->GetMinimiserConfiguration();
		theFunction = xml->GetFitFunctionConfiguration();
		studyParameters = xml->GetFitParameters();
		pdfsAndData = xml->GetPDFsAndData();
		numberStudies = xml->GetNumberRepeats();
		allConstraints = xml->GetConstraints();

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
	else
	{
		cerr << "XML file (" << FileName << ") not found" << endl;
	}
}

//Constructor with correct arguments
ToyStudy::ToyStudy( MinimiserConfiguration * TheMinimiser, FitFunctionConfiguration * TheFunction, ParameterSet * StudyParameters, vector< PDFWithData* > PDFsAndData, vector< ConstraintFunction* > InputConstraints, int NumberStudies )
: pdfsAndData(PDFsAndData), studyParameters(StudyParameters), theMinimiser(TheMinimiser), theFunction(TheFunction), numberStudies(NumberStudies), allConstraints(InputConstraints)
{
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
	delete allResults;
}

//Automate the toy study
ToyStudyResult * ToyStudy::DoWholeStudy( bool force_n_studies )
{
	//Make a vector of unique parameter names
	vector<string> uniqueNames;
	for (unsigned int pdfIndex = 0; pdfIndex < pdfsAndData.size(); pdfIndex++ )
	{
		//This is not strictly necessary, but suppresses a warning message
		pdfsAndData[pdfIndex]->SetPhysicsParameters(studyParameters);

		uniqueNames = StringProcessing::CombineUniques( uniqueNames, pdfsAndData[pdfIndex]->GetPDF()->GetPrototypeParameterSet() );
	}
	allResults = new ToyStudyResult(uniqueNames);

	//Loop over all studies
	for ( int studyIndex = 0; studyIndex < numberStudies; studyIndex++ )
	{
		cout << "\n\n\t\tStarting ToyStudy\t\t" << studyIndex+1 << "\tof:\t" << numberStudies << endl;
		allResults->StartStopwatch();
		allResults->AddFitResult( FitAssembler::DoSafeFit( theMinimiser, theFunction, studyParameters, pdfsAndData, allConstraints ) );
		if( ( allResults->GetFitResult( studyIndex )->GetFitStatus() < 0 ) && force_n_studies )
		{
			cout << "Fit fell over!\t Requesting another fit but keeping this result." << endl;
			numberStudies++;
		}
	}

	return allResults;
}

//Get the result of the toy study
ToyStudyResult * ToyStudy::GetToyStudyResult()
{
	return allResults;
}
