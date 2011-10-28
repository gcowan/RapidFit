/**
        @class ClassLookUp

        Central place to hold the methods for returning an instance of a class with a given name.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef CLASS_LOOK_UP_H
#define CLASS_LOOK_UP_H

//	RapidFit Headers
#include "IPDF.h"
#include "BasePDF.h"
#include "PDFConfigurator.h"
#include "FitFunction.h"
#include "IMinimiser.h"
#include "IDataGenerator.h"
#include "IPrecalculator.h"
//	System Headers
#include <iostream>
#include <vector>

using namespace std;

class ClassLookUp
{
	public:
		//	Function to return a named PDF by looking for it's symbol in the compiled object
		static IPDF* LookUpPDFName( string, vector<string>, vector<string>, PDFConfigurator* );
		//	Copy the input PDF using the correct copy constructor object
		static IPDF* CopyPDF( IPDF* );

		static FitFunction * LookUpFitFunctionName( string );
		static IMinimiser * LookUpMinimiserName( string, int );
		static IDataGenerator * LookUpDataGenerator( string, PhaseSpaceBoundary*, IPDF* );
		static IPrecalculator * LookUpPrecalculator( string, IPDF*, IPDF*, vector<ParameterSet*>, string );

		//	Return the path to the running executable
		static char* getSelfPath();
		//	Resolve the corresponding object within the library (or PIE executable)
		static void* getObject( string );
};

#endif
