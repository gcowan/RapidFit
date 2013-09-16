/*!
 * @ingroup Generators
 * @interface IDataGenerator
 *
 * @brief Interface for any class used for generating toy data from a PDF
 *
 * @author Benjamin M Wynne bwynne@cern.ch
*/

#pragma once
#ifndef I_DATA_GENERATOR_H
#define I_DATA_GENERATOR_H

//	RapidFit Headers
#include "IDataSet.h"

using namespace::std;

class IDataGenerator
{
	public:
		/*!
		 * Interface Function:
		 * Request a DataSet of a given size be Constructed
		 */
		virtual int GenerateData(int) = 0;

		/*!
		 * Interface Function:
		 * Request a Pointer to the newly created DataSet
		 */
		virtual IDataSet * GetDataSet() const = 0;

		/*!
		 * Virtual Destructor
		 */
		virtual ~IDataGenerator() {};

	protected:
		IDataGenerator() {};
};

#endif


