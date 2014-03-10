/*!
 * @class PDFWithData
 *
 * A class for creating/storing a PDF and its associated data set
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

//	RapidFit Headers
#include "PDFWithData.h"
#include "ClassLookUp.h"
//	System Headers
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace::std;

//Constructor with correct arguments
PDFWithData::PDFWithData( IPDF * InputPDF, PhaseSpaceBoundary * InputBoundary, DataSetConfiguration* DataConfig ) :
	fitPDF(ClassLookUp::CopyPDF(InputPDF)), inputBoundary(new PhaseSpaceBoundary(*InputBoundary)),  parametersAreSet(false),
	dataSetMaker(), cached_data(NULL), useCache(false)
{
	dataSetMaker = new DataSetConfiguration( *DataConfig );
}

//Destructor
PDFWithData::~PDFWithData()
{
	fitPDF->Can_Remove_Cache( true );
	if( fitPDF != NULL ) delete fitPDF;
	if( inputBoundary != NULL ) delete inputBoundary;
	if( dataSetMaker != NULL ) delete dataSetMaker;
	this->ClearCache();
}

//Return the PDF
IPDF * PDFWithData::GetPDF() const
{
	return fitPDF;
}

void PDFWithData::AddCachedData( IDataSet* input_cache )
{
	cached_data = input_cache;
}

DataSetConfiguration* PDFWithData::GetDataSetConfig()
{
	return dataSetMaker;
}

void PDFWithData::ClearCache()
{
	if( cached_data != NULL ) delete cached_data;
	cached_data = NULL;
}

void PDFWithData::SetUseCache( bool input )
{
	useCache = input;
}

bool PDFWithData::GetUseCache() const
{
	return useCache;
}

//Return the data set associated with the PDF
IDataSet * PDFWithData::GetDataSet() const
{
	//Combine all data sources
	IDataSet * newDataSet=NULL;

	//Right, the behaviour of this class has changed vastly from what  used to happen, but the external use if still the same.

	//For any dataset NOT a File, calling this class generates a new dataset and passes a pointer

	//For a FILE based dataset the dataset is loaded once, and then kept in memory until it's removed

	/*!
	 * This makes the Implcit assumption you don't mix File and Toy dataTypes!
	 *
	 * This may not be strictly true!
	 */
	if( dataSetMaker->GetSource() == "File" ) useCache = true;

	if( cached_data == NULL || !useCache )
	{
		newDataSet = dataSetMaker->MakeDataSet( inputBoundary, fitPDF );
		cached_data = newDataSet;
	}
	else
	{
		newDataSet = cached_data;
	}

	newDataSet->Print();

	cout << endl;
	//cout << "Providing DataSet At: " << newDataSet << endl;

	return newDataSet;
}

//Set the physics parameters of the PDF
bool PDFWithData::SetPhysicsParameters( ParameterSet* NewParameters )
{
	//	I am in the process of adding more flexibility to the RapidFit structure
	//	As such this requires that the Paramaters from the XML be passed in vectors
	//	I will update this function in time to stop refering to just the first element it sees

	//Set the parameters for the stored PDF and all data set makers
	fitPDF->UpdatePhysicsParameters( NewParameters );
	dataSetMaker->SetPhysicsParameters( NewParameters );

	return true;
}

void PDFWithData::Print() const
{
	cout << "PDFWithData information:" << endl;
	cout << "This is to be coded up when needed" << endl;
}

string PDFWithData::XML() const
{
	stringstream xml;

	xml << "<ToFit>" << endl;

	xml << fitPDF->XML() << endl;

	xml << endl;
	xml << dataSetMaker->XML() << endl;
	xml << endl;

	xml << "</ToFit>" << endl;

	return xml.str();
}


