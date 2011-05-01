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
	delete dataBoundary;
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
