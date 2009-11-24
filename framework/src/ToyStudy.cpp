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
		minimiserName = xml->GetMinimiserName();
		functionName = xml->GetFitFunctionName();
		studyParameters = xml->GetFitParameters();
		pdfsAndData = xml->GetPDFsAndData();
		numberStudies = xml->GetNumberRepeats();

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
ToyStudy::ToyStudy( string MinimiserName, string FunctionName, ParameterSet * StudyParameters, vector< PDFWithData* > PDFsAndData, int NumberStudies )
	: minimiserName(MinimiserName), functionName(FunctionName), studyParameters(StudyParameters), pdfsAndData(PDFsAndData), numberStudies(NumberStudies)
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
ToyStudyResult * ToyStudy::DoWholeStudy()
{
	//Make a vector of unique parameter names
	vector<string> uniqueNames;
	for ( int pdfIndex = 0; pdfIndex < pdfsAndData.size(); pdfIndex++ )
	{
		//This is not strictly necessary, but suppresses a warning message
		pdfsAndData[pdfIndex]->SetPhysicsParameters(studyParameters);

		uniqueNames = StringProcessing::CombineUniques( uniqueNames, pdfsAndData[pdfIndex]->GetPDF()->GetPrototypeParameterSet() );
	}
	allResults = new ToyStudyResult(uniqueNames);

	//Loop over all studies
	for ( int studyIndex = 0; studyIndex < numberStudies; studyIndex++ )
	{
		cout << "Starting study #" << studyIndex << endl;
		allResults->StartStopwatch();
		allResults->AddFitResult( FitAssembler::DoFit( minimiserName, functionName, studyParameters, pdfsAndData ) );
	}

	return allResults;
}

//Get the result of the toy study
ToyStudyResult * ToyStudy::GetToyStudyResult()
{
	return allResults;
}
