/**
        @class MemoryDataSet

        A data set which simply stores a vector of pointers to datapoint objects

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef MEMORY_DATA_SET_H
#define MEMORY_DATA_SET_H

//	RapidFit Headers
#include "IDataSet.h"
#include "DataPoint.h"
//	System Headers
#include <vector>

using namespace std;

class MemoryDataSet : public IDataSet
{
	public:
		MemoryDataSet();
		MemoryDataSet( PhaseSpaceBoundary* );
		~MemoryDataSet();

		//Interface functions
		virtual DataPoint * GetDataPoint(int);
		virtual void ReserveDataSpace( int numberOfPoints );
		virtual bool AddDataPoint( DataPoint* );
		virtual int GetDataNumber();
		virtual PhaseSpaceBoundary * GetBoundary();

		virtual void SortBy( string );

		void Clear();

	private:
		//	Uncopyable!
		//MemoryDataSet ( const MemoryDataSet& );
		//MemoryDataSet& operator = ( const MemoryDataSet& );
		vector<DataPoint> allData;
		PhaseSpaceBoundary * dataBoundary;

};

bool compare_datapoints ( pair<DataPoint,pair<string,int> > first, pair<DataPoint,pair<string,int> > second );

#endif
