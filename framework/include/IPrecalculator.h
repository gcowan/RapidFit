/**
  @Interface IPrecalculator

  An interface for classes that perform some operation on an input IDataSet before it is used in the fit

  @author Benjamin Wynne bwynne@cern.ch
  @date 2009-12-14
  */

#ifndef I_PRECALCULATOR_H
#define I_PRECALCULATOR_H

//	RapidFit Headers
#include "IDataSet.h"

class IPrecalculator
{
	public:
		virtual IDataSet * ProcessDataSet( IDataSet* ) = 0;
};

#endif
