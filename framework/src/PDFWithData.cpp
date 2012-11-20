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
PDFWithData::PDFWithData( IPDF * InputPDF, PhaseSpaceBoundary * InputBoundary, vector< DataSetConfiguration* > DataConfig ) :
	fitPDF(ClassLookUp::CopyPDF(InputPDF)), inputBoundary(new PhaseSpaceBoundary(*InputBoundary)),  parametersAreSet(false),
	dataSetMakers(), cached_data(), useCache(false), debug( new DebugClass(false) ), WeightName(""), useWeights(false)
{
	if( DataConfig.size() < 1 || DataConfig.empty() )
	{
		cerr << "No data sets configured" << endl;
		exit(1);
	}
	for( vector<DataSetConfiguration*>::const_iterator config_i = DataConfig.begin(); config_i != DataConfig.end(); ++config_i )
	{
		dataSetMakers.push_back( new DataSetConfiguration( *(*config_i) ) );
	}
}

//Destructor
PDFWithData::~PDFWithData()
{
	fitPDF->Can_Remove_Cache( true );
	if( fitPDF != NULL ) delete fitPDF;
	if( inputBoundary != NULL ) delete inputBoundary;
	while( !dataSetMakers.empty() )
	{
		if( dataSetMakers.back() != NULL ) delete dataSetMakers.back();
		dataSetMakers.pop_back();
	}
	this->ClearCache();
	if( debug != NULL ) delete debug;
}

//Return the PDF
IPDF * PDFWithData::GetPDF() const
{
	return fitPDF;
}

void PDFWithData::AddCachedData( vector<IDataSet*> input_cache )
{
	for( unsigned short int element=0; element < input_cache.size(); ++element )
	{
		cached_data.push_back( input_cache[element] );
	}
}

void PDFWithData::AddCachedData( IDataSet* input_cache )
{
	cached_data.push_back( input_cache );
}


DataSetConfiguration* PDFWithData::GetDataSetConfig( int i )
{
	return dataSetMakers[(unsigned)i];
}

vector<DataSetConfiguration*> PDFWithData::GetAllDataSetConfigs()
{
	return dataSetMakers;
}

void PDFWithData::ClearCache()
{
	//cout << "CacheSize: " << cached_data.size() << endl;

	//	Problems Happen so protect against it
	sort( cached_data.begin(), cached_data.end() );
	cached_data.erase( unique( cached_data.begin(), cached_data.end() ), cached_data.end() );

	while( !cached_data.empty() )
	{
		//cout << "Removing DataSet At: " << cached_data.back() << endl;
		if( cached_data.back() != NULL ) delete cached_data.back();
		cached_data.pop_back();
	}
}

void PDFWithData::SetUseCache( bool input )
{
	useCache = input;
}

bool PDFWithData::GetUseCache() const
{
	return useCache;
}

vector<IDataSet*> PDFWithData::GetCacheList()
{
	return cached_data;
}

IDataSet* PDFWithData::GetFromCache( int request )
{
	if( (unsigned)request >= cached_data.size() ) return NULL;
	else if( request < 0 ) return NULL;
	else return cached_data[(unsigned)request];
}

void PDFWithData::RemoveFromCache( int request )
{
	if( request > 0 && (unsigned)request < cached_data.size() )
	{
		int i=0;
		vector<IDataSet*>::iterator table_iter = cached_data.begin();
		while( i != request )
		{
			++table_iter;
			++i;
		}
		cached_data.erase( table_iter );
	}
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
	if( dataSetMakers[0]->GetSource() == "File" ) useCache = true;

	if( cached_data.empty() || !useCache )
	{
		newDataSet = dataSetMakers[0]->MakeDataSet( inputBoundary, fitPDF );
		for( unsigned int sourceIndex = 1; sourceIndex < dataSetMakers.size(); ++sourceIndex )
		{
			IDataSet * extraData = dataSetMakers[sourceIndex]->MakeDataSet( inputBoundary, fitPDF );
			for( int dataIndex = 0; dataIndex < extraData->GetDataNumber(); ++dataIndex )
			{
				newDataSet->AddDataPoint( extraData->GetDataPoint(dataIndex) );
			}
			delete extraData;
		}

		cached_data.push_back( newDataSet );
	}
	else
	{
		newDataSet = cached_data.back();
	}

	if( useWeights ) newDataSet->UseEventWeights( WeightName );

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
	for (unsigned int dataIndex = 0; dataIndex < dataSetMakers.size(); ++dataIndex )
	{
		dataSetMakers[dataIndex]->SetPhysicsParameters( NewParameters );
	}
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

	for( vector<DataSetConfiguration*>::const_iterator set_i = dataSetMakers.begin(); set_i != dataSetMakers.end(); ++set_i )
	{
		xml << endl;
		xml << (*set_i)->XML() << endl;
		xml << endl;
	}

	xml << "</ToFit>" << endl;

	return xml.str();
}

void PDFWithData::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
}

void PDFWithData::UseEventWeights( const string Name )
{
	WeightName = Name;
	useWeights = true;
}

bool PDFWithData::GetWeightsWereUsed() const
{
	return useWeights;
}

string PDFWithData::GetWeightName() const
{
	return WeightName;
}

