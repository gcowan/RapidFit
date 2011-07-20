//
/**
        @class MemoryDataSet

        A data set which simply stores a vector of pointers to datapoint objects

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "MemoryDataSet.h"
//	System Headers
#include <iostream>
#include <vector>
#include <algorithm>

//Default constructor
MemoryDataSet::MemoryDataSet() : allData(), dataBoundary()
{
}

//Constructor with correct argument
MemoryDataSet::MemoryDataSet( PhaseSpaceBoundary * NewBoundary ) : allData(), dataBoundary(NewBoundary)
{
}

//Destructor
MemoryDataSet::~MemoryDataSet()
{
//	delete dataBoundary;
}

void MemoryDataSet::ReserveDataSpace( int numberOfPoints )
{
	allData.reserve( unsigned(numberOfPoints) );
}

//Add a data point to the set
bool MemoryDataSet::AddDataPoint( DataPoint * NewDataPoint )
{
	if ( dataBoundary->IsPointInBoundary(NewDataPoint) )
	{
		allData.push_back( *NewDataPoint );
		return true;
	}
	else
	{
		//cerr << "Data point is not within data set boundary" << endl;
		return false;
	}
}

//Retrieve the data point with the given index
DataPoint * MemoryDataSet::GetDataPoint(  int Index )
{
	if ( Index < int(allData.size()) )
	{
		return &allData[unsigned(Index)];
	}
	else
	{
		cerr << "Index (" << Index << ") out of range in DataSet" << endl;
	}
	return NULL;
}

//Get the number of data points in the set
int MemoryDataSet::GetDataNumber()
{
	return int(allData.size());
}

//Get the data bound
PhaseSpaceBoundary * MemoryDataSet::GetBoundary()
{
	return dataBoundary;
}

//Empty the data set
void MemoryDataSet::Clear()
{
	allData.clear();
}

void MemoryDataSet::SortBy( string parameter )
{
	cout << "Sorting" << endl;
	if( allData.size() > 0 )
	{
		vector<pair<DataPoint,pair<string,int> > > allData_sort;

		for( vector<DataPoint>::iterator data_i = allData.begin(); data_i != allData.end(); ++data_i )
		{
			allData_sort.push_back( make_pair( *data_i, make_pair( parameter, -1 ) ) );
		}

		cout << "hello" << endl;
		sort( allData_sort.begin(), allData_sort.end(), DataPoint() );
		cout << "sorted" << endl;
		//	Sort the data in memory

		while( !allData.empty() ) allData.pop_back();		

		for( vector<pair<DataPoint,pair<string,int> > >::iterator sort_i = allData_sort.begin(); sort_i != allData_sort.end(); ++sort_i )
		{
			allData.push_back( sort_i->first );
		}

		cout << allData.size() << endl;
	}
	cout << "Sorted" << endl;
}

