/**
  @class Blinder

  Functions to blind a parameter 
 
  @author Pete Clarke
  @date 2010-11-25
 */

//	RapidFit Headers
#include "Blinder.h"

using namespace::std;

//Unblinds a given value -probably not used.
double Blinder::unBlindParameter( double blindValue, const char * blindString, double scale )
{
	RooRealVar bv("b", "b", blindValue ) ;
	RooUnblindUniform tv("t", "t", blindString , scale, bv ) ;
	RooAbsReal& tv1 = tv ;
	return tv1.getVal() ;	
}


// Returns a blind offset relative to zero input as the blinded parameter.
double Blinder::getBlindOffset( const char * blindString, double scale )
{
	double blindValue = 0;
	RooRealVar bv("b", "b", blindValue ) ;
	RooUnblindUniform tv("t", "t", blindString , scale, bv ) ;
	RooAbsReal& tv1 = tv ;
	return tv1.getVal() ;	
}

