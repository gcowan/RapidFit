/*!
 * @class Blinder
 *
 * @brief Code to carry out blinding
 *
 * @author Pete Clarke 
*/

#pragma once
#ifndef BLINDER_RESULT_H
#define BLINDER_RESULT_H

///	ROOT Headers
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooUnblindUniform.h"
///	System Headers
#include <string>

using namespace::std;

class Blinder
{
	public:
		/*!
		 * @brief Function to unblind a blinded value
		 * 
		 * @param blindValue   Blinded Value
		 *
		 * @param blindString  Blinding String as defined in the XML
		 *
		 * @param scale        Blinding Scale as defined in the XML
		 */
		static double unBlindParameter( double blindValue, const char * blindString, double scale );

		/*!
		 * @brief Function to get the blinding offset defined with the string and scale
		 *
		 * @param blindString  Blinding String as defined in the XML
		 *
		 * @param scale        Blinding Scale as defined in the XML
		 */
		static double getBlindOffset( const char * blindString, double scale );	
};

#endif

