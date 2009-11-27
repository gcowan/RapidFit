/**
        @class ClassLookUp

        Central place to hold the methods for returning an instance of a class with a given name.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef CLASS_LOOK_UP_H
#define CLASS_LOOK_UP_H

#include "IPDF.h"
#include "FitFunction.h"
#include "IMinimiser.h"
#include "IDataGenerator.h"
#include <iostream>

using namespace std;

class ClassLookUp
{
	public:
		static IPDF * LookUpPDFName( string, vector<string>, vector<string> );
		static FitFunction * LookUpFitFunctionName( string, string );
		static IMinimiser * LookUpMinimiserName( string, int );
		static IDataGenerator * LookUpDataGenerator( string, PhaseSpaceBoundary*, IPDF* );
};

#endif
