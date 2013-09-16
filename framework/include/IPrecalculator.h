/*!
 * @Interface IPrecalculator
 *
 * An interface for classes that perform some operation on an input IDataSet before it is used in the fit
 *
 * @author Benjamin Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef I_PRECALCULATOR_H
#define I_PRECALCULATOR_H

///	RapidFit Headers
#include "IDataSet.h"
#include "IPDF.h"

using namespace::std;

class IPrecalculator
{
	public:
		/*!
		 * @brief Interface Function to Process the given DataSet and add a set of Weights
		 *
		 * @param Input   This is a pointer to the DataSet wished to be processed
		 */
		virtual IDataSet * ProcessDataSet( IDataSet* Input, IPDF* InputPDF ) = 0;

		virtual void SetApplyAlphaCorrection( bool ) = 0;

		/*!
		 * @brief Virtual Destructor
		 */
		virtual ~IPrecalculator() {};

	protected:
		IPrecalculator() {};
};

#endif

