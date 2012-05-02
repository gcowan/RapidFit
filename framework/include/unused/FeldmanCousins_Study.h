/*!
 * @class FeldmanCousins_Study
 *
 * Initial implementation of FelmanCousins Study
 * Linear code, difficult to understand, useful for reference but no longer called
 *
 */

#ifndef FeldmanCousins_Study_H
#define FeldmanCousins_Study_H

///	RapidFit Headers
#include "FitResultVector.h"
#include "OutputConfiguration.h"
#include "MinimiserConfiguration.h"
#include "XMLConfigReader.h"
#include "PDFWithData.h"
///	System Headers
#include <vector>

using namespace::std;

class FeldmanCousins_Study
{

	public:
		/*!
		 * Std Feldman-Cousins Code, relying on many already defined objects even though it pulls in the XMLConfigReader object
		 */
		static FitResultVector* FeldmanCousins( FitResultVector*, FitResultVector*, const vector<unsigned int>, const unsigned int, const bool, OutputConfiguration*,
				MinimiserConfiguration*, FitFunctionConfiguration*, XMLConfigReader*, const vector< PDFWithData* >, const int=-999);
};

#endif

