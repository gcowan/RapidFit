/**
        @class InputParsing

        Collection of static functions for creating RapidFit data objects from string templates

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef INPUT_PARSING_H
#define INPUT_PARSING_H

//	RapidFit Headers
#include "ParameterSet.h"
#include "PhaseSpaceBoundary.h"
#include "PDFWithData.h"
//	System Headers
#include <string>
#include <vector>

class InputParsing
{
	public:
		static ParameterSet * MakeParameterSet( string );
		static ParameterSet * MakeParameterSet( vector<string> );
		static PhaseSpaceBoundary * MakePhaseSpaceBoundary( string );
		static PDFWithData * MakePDFWithData( string, string, string );
};

#endif
