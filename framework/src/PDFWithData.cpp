/**
  @class PDFWithData

  A class for creating/storing a PDF and its associated data set

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-5
 */

#include "PDFWithData.h"
#include "DataFileLoader.h"
#include "ClassLookUp.h"
#include <stdlib.h>
#include <iostream>

using namespace std;

//Default constructor
PDFWithData::PDFWithData() : parametersAreSet(false)
{
}

//Constructor with correct aruments
PDFWithData::PDFWithData( IPDF * InputPDF, PhaseSpaceBoundary * InputBoundary, vector< DataSetConfiguration* > DataConfig, vector< IPrecalculator* > InputPrecalculators ) : fitPDF(InputPDF),
	inputBoundary(InputBoundary), dataSetMakers(DataConfig), parametersAreSet(false), dataProcessors(InputPrecalculators)
{
	if ( DataConfig.size() < 1 )
	{
		cerr << "No data sets configured" << endl;
		exit(1);
	}       
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
	//Combine all data sources
	IDataSet * newDataSet = dataSetMakers[0]->MakeDataSet( inputBoundary, fitPDF );
	for ( int sourceIndex = 1; sourceIndex < dataSetMakers.size(); sourceIndex++ )
	{
		IDataSet * extraData = dataSetMakers[sourceIndex]->MakeDataSet( inputBoundary, fitPDF );
		for ( int dataIndex = 0; dataIndex < extraData->GetDataNumber(); dataIndex++ )
		{
			newDataSet->AddDataPoint( extraData->GetDataPoint(dataIndex) );
		}
		delete extraData;
	}

	//Precalculation, if required
	for ( int precalculatorIndex = 0; precalculatorIndex < dataProcessors.size(); precalculatorIndex++ )
	{
		IDataSet * oldDataSet = newDataSet;
		newDataSet = dataProcessors[precalculatorIndex]->ProcessDataSet(oldDataSet);
		delete oldDataSet;
	}

	cout << "DataSet contains " << newDataSet->GetDataNumber() << " events" << endl;
	return newDataSet;
}

//Set the physics parameters of the PDF
bool PDFWithData::SetPhysicsParameters( ParameterSet * NewParameters )
{
	//Set the parameters for the stored PDF and all data set makers
	bool success = fitPDF->SetPhysicsParameters(NewParameters);
	for ( int dataIndex = 0; dataIndex < dataSetMakers.size(); dataIndex++ )
	{
		success &= dataSetMakers[dataIndex]->SetPhysicsParameters(NewParameters);
	}

	if (success)
	{
		parametersAreSet = true;
		return true;
	}
	else
	{
		cerr << "Failed to set PDF parameters in initialisation" << endl;
		exit(1);
	}
}
