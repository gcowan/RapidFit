/**
	@class PDFWithData

	A class for creating/storing a PDF and its associated data set

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-5
*/

#include "PDFWithData.h"
#include "DataFileLoader.h"
#include "ClassLookUp.h"
#include "InvalidObject.h"
#include <iostream>

using namespace std;

//Default constructor
PDFWithData::PDFWithData() : parametersAreSet(false)
{
}

//Constructor with correct aruments
PDFWithData::PDFWithData( IPDF * InputPDF, string DataSource, long DataAmount, vector<string> DataArguments, PhaseSpaceBoundary * InputBoundary )
	: generatePDF(InputPDF), fitPDF(InputPDF), dataSource(DataSource), dataAmount(DataAmount), dataArguments(DataArguments), inputBoundary(InputBoundary), parametersAreSet(false)
{
}

//Constructor with correct aruments
PDFWithData::PDFWithData( IPDF * GeneratePDF, IPDF * FitPDF, string DataSource, long DataAmount, vector<string> DataArguments, PhaseSpaceBoundary * InputBoundary )
	: generatePDF(GeneratePDF), fitPDF(FitPDF), dataSource(DataSource), dataAmount(DataAmount), dataArguments(DataArguments),
	inputBoundary(InputBoundary), parametersAreSet(false)
{
}
//Destructor
PDFWithData::~PDFWithData()
{
}

//Return the PDF
IPDF * PDFWithData::GetPDF()
{
	if (!parametersAreSet)
	{
		cout << "Warning: PDF parameters have not yet been set" << endl;
	}
	return fitPDF;
}

//Return the data set associated with the PDF
IDataSet * PDFWithData::GetDataSet()
{
	//Some kind of decision about what kind of data set to use?
	IDataSet * newDataSet;

	if ( dataSource == "File" )
	{
		DataFileLoader * dataLoader = new DataFileLoader( dataArguments[0], inputBoundary, dataAmount );
		newDataSet = dataLoader->GetDataSet();
		delete dataLoader;
	}
	else
	{
		//PDF parameters must be set before data can be generated
		if (parametersAreSet)
		{
			//Assume it's an accept/reject generator, or some child of it
			IDataGenerator * dataGenerator = ClassLookUp::LookUpDataGenerator( dataSource, inputBoundary, generatePDF );
			dataGenerator->GenerateData(dataAmount);
			newDataSet = dataGenerator->GetDataSet();
			delete dataGenerator;
		}
		else
		{
			cerr << "PDF parameters must be set before data can be generated" << endl;
			newDataSet = new InvalidObject("PDF parameters must be set before data can be generated");
		}
	}

	return newDataSet;
}

//Set the physics parameters of the PDF
bool PDFWithData::SetPhysicsParameters( ParameterSet * NewParameters )
{
	if ( fitPDF->SetPhysicsParameters(NewParameters) && generatePDF->SetPhysicsParameters(NewParameters) )
	{
		parametersAreSet = true;
		return true;
	}
	else
	{
		cerr << "Failed to set PDF physics parameters" << endl;
		return false;
	}
}
