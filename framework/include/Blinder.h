/**
        @class Blinder

        Code to carry out blinding

        @author Pete Clarke 
	@date 2010-11-25
*/


#ifndef BLINDER_RESULT_H
#define BLINDER_RESULT_H

//	ROOT Headers
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooUnblindUniform.h"
//	System Headers
#include <string>

using namespace std;

class Blinder
{
	public:
	static double unBlindParameter( double blindValue, const char * blindString, double scale );
	static double getBlindOffset( const char * blindString, double scale );	
};

#endif
