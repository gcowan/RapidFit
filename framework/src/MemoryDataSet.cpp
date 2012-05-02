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
#include <math.h>

#define DOUBLE_TOLERANCE_DATA 1E-6

using namespace::std;

//Constructor with correct argument
MemoryDataSet::MemoryDataSet( PhaseSpaceBoundary * NewBoundary ) : allData(), dataBoundary( new PhaseSpaceBoundary(*NewBoundary) )
{
}

//Destructor
MemoryDataSet::~MemoryDataSet()
{
	delete dataBoundary;
	while( !allData.empty() )
	{
		if( allData.back() != NULL ) delete allData.back();
		allData.pop_back();
	}
}

void MemoryDataSet::ReserveDataSpace( int numberOfPoints )
{
	allData.reserve( unsigned(numberOfPoints) );
}

//Add a data point to the set
bool MemoryDataSet::AddDataPoint( DataPoint* NewDataPoint )
{
	if ( dataBoundary->IsPointInBoundary(NewDataPoint) )
	{
		allData.push_back( NewDataPoint );
		return true;
	}
	else
	{
		delete NewDataPoint;
		//cerr << "Data point is not within data set boundary" << endl;
		return false;
	}
}

//Retrieve the data point with the given index
DataPoint * MemoryDataSet::GetDataPoint( int Index ) const
{
	if ( Index < int(allData.size()) )
	{
		return allData[unsigned(Index)];
	}
	else
	{
		cerr << "Index (" << Index << ") out of range in DataSet" << endl;
	}
	return NULL;
}

//Get the number of data points in the set
int MemoryDataSet::GetDataNumber() const
{
	return int(allData.size());
}

//Get the data bound
PhaseSpaceBoundary * MemoryDataSet::GetBoundary() const
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
		vector<pair<DataPoint*,pair<string,int> > > allData_sort;

		for( vector<DataPoint*>::iterator data_i = allData.begin(); data_i != allData.end(); ++data_i )
		{
			allData_sort.push_back( make_pair( *data_i, make_pair( parameter, -1 ) ) );
		}

		//cout << "hello" << endl;
		sort( allData_sort.begin(), allData_sort.end(), DataPoint() );
		cout << "sorted" << endl;
		//	Sort the data in memory

		while( !allData.empty() ) allData.pop_back();

		for( vector<pair<DataPoint*,pair<string,int> > >::iterator sort_i = allData_sort.begin(); sort_i != allData_sort.end(); ++sort_i )
		{
			allData.push_back( sort_i->first );
		}

		cout << allData.size() << endl;
	}
	cout << "Sorted" << endl;
}

vector<DataPoint*> MemoryDataSet::GetDiscreteSubSet( vector<string> discreteParam, vector<double> discreteVal )
{
	vector<DataPoint*> returnable_subset;
	if( discreteParam.size() != discreteVal.size() )
	{
		cerr << "\n\n\t\tBadly Defined definition of a subset, returning 0 events!\n\n" << endl;
		return returnable_subset;
	}

	for( unsigned int i=0; i< allData.size(); ++i )
	{
		DataPoint* data_i = allData[i];
		vector<ObservableRef*> temp_ref;
		for( unsigned int j=0; j< discreteParam.size(); ++j )
		{
			temp_ref.push_back( new ObservableRef( discreteParam[j] ) );
		}

		bool decision = true;
		for( unsigned int j=0; j< discreteParam.size(); ++j )
		{
			if( !( fabs( data_i->GetObservable( *(temp_ref[j]) )->GetValue() - discreteVal[j] ) < DOUBLE_TOLERANCE_DATA ) )
			{
				decision = false;
			}
		}

		if( decision ) returnable_subset.push_back( data_i );

		while( !temp_ref.empty() )
		{
			if( temp_ref.back() != NULL ) delete temp_ref.back();
			temp_ref.pop_back();
		}
	}

	return returnable_subset;
}

void MemoryDataSet::Print()
{
	for( vector<DataPoint*>::iterator data_i = allData.begin(); data_i != allData.end(); ++data_i )
	{
		(*data_i)->Print();
	}
}

int MemoryDataSet::Yield()
{
	return this->GetDataNumber();
}

