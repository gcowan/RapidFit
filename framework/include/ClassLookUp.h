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
		static IPDF * LookUpPDFName( string, vector<string>, vector<string>, PDFConfigurator );
		static FitFunction * LookUpFitFunctionName( string );
		static IMinimiser * LookUpMinimiserName( string, int );
		static IDataGenerator * LookUpDataGenerator( string, PhaseSpaceBoundary*, IPDF* );
		static IPrecalculator * LookUpPrecalculator( string, IPDF*, IPDF*, vector<ParameterSet*>, string );
};

#endif
