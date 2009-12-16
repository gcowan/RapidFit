/**
        @class MemoryDataSet

        A data set which simply stores a vector of pointers to datapoint objects

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef MEMORY_DATA_SET_H
#define MEMORY_DATA_SET_H

#include "IDataSet.h"
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
		virtual bool AddDataPoint( DataPoint* );
		virtual int GetDataNumber();
		virtual PhaseSpaceBoundary * GetBoundary();

		void Clear();

	private:
		vector<DataPoint> allData;
		PhaseSpaceBoundary * dataBoundary;
};

#endif
