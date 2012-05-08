/**
  @class MemoryDataSet

  A data set which simply stores a vector of pointers to datapoint objects

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

#pragma once
#ifndef MEMORY_DATA_SET_H
#define MEMORY_DATA_SET_H

//	RapidFit Headers
#include "IDataSet.h"
#include "DataPoint.h"
//	System Headers
#include <vector>

using namespace::std;

class MemoryDataSet : public IDataSet
{
	public:
		MemoryDataSet( PhaseSpaceBoundary* );
		~MemoryDataSet();

		//Interface functions
		virtual DataPoint * GetDataPoint(int) const;
		virtual void ReserveDataSpace( int numberOfPoints );
		virtual bool AddDataPoint( DataPoint* );
		virtual int GetDataNumber( DataPoint* templateDataPoint =NULL ) const;
		virtual PhaseSpaceBoundary * GetBoundary() const;

		virtual void SortBy( string );

		virtual vector<DataPoint*> GetDiscreteSubSet( vector<string> discreteParam, vector<double> discreteVal );

		void Clear();

		virtual void Print();

		/*!                      
		 * @brief Returns an estimate of the total Yield
		 */                      
		virtual int Yield();

	private:
		//	Uncopyable!
		MemoryDataSet ( const MemoryDataSet& );
		MemoryDataSet& operator = ( const MemoryDataSet& );
		vector<DataPoint*> allData;
		PhaseSpaceBoundary * dataBoundary;
		mutable vector<int> allSubSets;
};

bool compare_datapoints ( pair<DataPoint,pair<string,int> > first, pair<DataPoint,pair<string,int> > second );

#endif

