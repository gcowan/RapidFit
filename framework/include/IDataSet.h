/**
        @interface IDataSet

        Interface for all collections of data points

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef I_DATA_SET_H
#define I_DATA_SET_H

#include "DataPoint.h"
#include "PhaseSpaceBoundary.h"

class IDataSet
{
	public:
		virtual DataPoint * GetDataPoint(int) = 0;
		virtual bool AddDataPoint( DataPoint* ) = 0;
		virtual int GetDataNumber() = 0;
		virtual PhaseSpaceBoundary * GetBoundary() = 0;
};

#endif
