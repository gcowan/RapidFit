/*!
 * @class InputParsing
 *
 * @brief Collection of static functions for creating RapidFit data objects from string templates
 *
 * @author Benjamin M Wynne bwynne@cern.ch
*/

#pragma once
#ifndef INPUT_PARSING_H
#define INPUT_PARSING_H

//	RapidFit Headers
#include "ParameterSet.h"
#include "PhaseSpaceBoundary.h"
#include "PDFWithData.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class InputParsing
{
	public:
		/*!
		 * Make a ParameterSet from the input string
		 */
		static ParameterSet * MakeParameterSet( string );

		/*!
		 * Make a ParameterSet from a vector of strings
		 */
		static ParameterSet * MakeParameterSet( vector<string> );

		/*!
		 * Make a PhaseSpaceBoundary from a string
		 */
		static PhaseSpaceBoundary * MakePhaseSpaceBoundary( string );

		/*!
		 * Make a PDFWithData from a string inputs
		 */
		static PDFWithData * MakePDFWithData( string, string, string );
};

#endif

