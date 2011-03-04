// $Id: Bs2JpsiPhi_SignalAlt_BaseClass.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_BaseClass Bs2JpsiPhi_SignalAlt_BaseClass.cpp
 *
 *  Base class for Bs2JpsiPhi_SignalAlt..... PDFs
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-12
 */

#include "Bs2JpsiPhi_SignalAlt_BaseClass.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"
#include "Mathematics.h"

#define DEBUGFLAG true

//......................................
//Constructor

Bs2JpsiPhi_SignalAlt_BaseClass::Bs2JpsiPhi_SignalAlt_BaseClass() : 
	// Physics parameters
	  gammaName     ( "gamma" )
	, deltaGammaName( "deltaGamma" )
	, deltaMName    ( "deltaM")
	, Phi_sName     ( "Phi_s")
	, Azero_sqName  ( "Azero_sq" )
	, Aperp_sqName  ( "Aperp_sq" )
	, As_sqName		( "As_sq" )
	, delta_zeroName( "delta_zero" )
	, delta_paraName( "delta_para" )
	, delta_perpName( "delta_perp" )
	, delta_sName( "delta_s" )
	// Detector parameters
	, mistagName	( "mistag" )
	, res1Name	( "timeResolution1" )
	, res2Name	( "timeResolution2" )
	, res1FractionName	( "timeResolution1Fraction" )
	, timeOffsetName	( "timeOffset" )
	// Angular acceptance factors
	, angAccI1Name ( "angAccI1" )
	, angAccI2Name ( "angAccI2" )
	, angAccI3Name ( "angAccI3" )
	, angAccI4Name ( "angAccI4" )
	, angAccI5Name ( "angAccI5" )
	, angAccI6Name ( "angAccI6" )
	, angAccI7Name ( "angAccI7" )
	, angAccI8Name ( "angAccI8" )
	, angAccI9Name ( "angAccI9" )
	, angAccI10Name ( "angAccI10" )
	// Observables
	, timeName	    ( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	    ( "phi" )
	, cosPsiName	( "cosPsi" )
	, tagName	    ( "tag" )
	, timeAcceptanceCategoryName ( "timeAcceptanceCategory" )
	// Other things
{
	if( ! USE_LOWER_TIME_ACCEPTANCE ) {
		cout << "=====>WARNING " << endl ;
		cout << "=====>WARNING YOU APPEAR TO **NOT** BE USING THE LOWER TIME ACCEPTANCE" << endl ;
		cout << "       - The define flag is turned off " << endl << endl ;
	}

	if( USE_UPPER_TIME_ACCEPTANCE ) {
		cout << "=====>WARNING " << endl ;
		cout << "=====>WARNING YOU APPEAR TO BE USING THE UPPER TIME ACCEPTANCE. BE WARNED THAT: " << endl ;
		cout << "       - It doesnt work for Tagged fits as we havnt looked up how to combine with Exp*Cos in the normalisation " << endl ;
		cout << "       - It doesnt work with resolution as we dont know how to do the integral for normalisation " << endl << endl ;
	}
		
}


//........................................................
//Destructor
Bs2JpsiPhi_SignalAlt_BaseClass::~Bs2JpsiPhi_SignalAlt_BaseClass()
{
}

//....................................
//Internal helper functions

double Bs2JpsiPhi_SignalAlt_BaseClass::AT() const { 
	if( Aperp_sq <= 0. ) return 0. ;
	else return sqrt(Aperp_sq) ; 
};
double Bs2JpsiPhi_SignalAlt_BaseClass::AP() const { 
	if( Apara_sq <= 0. ) return 0. ;
	else return sqrt(Apara_sq) ; 
};
double Bs2JpsiPhi_SignalAlt_BaseClass::A0() const { 
	if( Azero_sq <= 0. ) return 0. ;
	else return sqrt(Azero_sq) ; 
};
double Bs2JpsiPhi_SignalAlt_BaseClass::AS() const { 
	if( As_sq <= 0. ) return 0. ;
	else return sqrt(As_sq) ; 
};

double Bs2JpsiPhi_SignalAlt_BaseClass::ctrsq() const { return (ctheta_tr*ctheta_tr) ; }
double Bs2JpsiPhi_SignalAlt_BaseClass::strsq() const { return (1.0 - ctheta_tr*ctheta_tr) ; }
double Bs2JpsiPhi_SignalAlt_BaseClass::ct1sq() const { return (ctheta_1*ctheta_1) ; }
double Bs2JpsiPhi_SignalAlt_BaseClass::st1sq() const { return (1.0 - ctheta_1*ctheta_1) ; }
double Bs2JpsiPhi_SignalAlt_BaseClass::cphsq() const { return (cos(phi_tr)*cos(phi_tr)) ; }
double Bs2JpsiPhi_SignalAlt_BaseClass::sphsq() const { return (sin(phi_tr)*sin(phi_tr)) ; }

double Bs2JpsiPhi_SignalAlt_BaseClass::gamma_l() const { 
	double gl = gamma() + ( dgam / 2.0 ) ;
	if( gl <= 0. ) {
		cout << " In Bs2JpsiPhi_SignalAlt_BaseClass : gamma_l() < 0 so setting it to 0.0000001 " << endl ;
		return 0.0000001 ;
	}
	else
		return gl ; 
}

double Bs2JpsiPhi_SignalAlt_BaseClass::gamma_h() const { 
	double gh = gamma() - ( dgam / 2.0 ) ;
	if( gh <= 0. ) {
		cout << " In Bs2JpsiPhi_SignalAlt_BaseClass : gamma_h() < 0 so setting it to 0.0000001 " << endl ;
		return 0.0000001 ;
	}
	else
		return gh ;   
}

double Bs2JpsiPhi_SignalAlt_BaseClass::gamma() const { return gamma_in ; }

double Bs2JpsiPhi_SignalAlt_BaseClass::q() const { return tag ;}

bool Bs2JpsiPhi_SignalAlt_BaseClass::useLowerTimeAcceptance() const { return (USE_LOWER_TIME_ACCEPTANCE && (timeAcceptanceCategory > 0)) ; }

bool Bs2JpsiPhi_SignalAlt_BaseClass::useUpperTimeAcceptance() const { return (USE_UPPER_TIME_ACCEPTANCE && ( fabs(UPPER_TIME_ACCEPTANCE_FACTOR - 0) < DOUBLE_TOLERANCE) ) ; }

double Bs2JpsiPhi_SignalAlt_BaseClass::upperTimeAcceptanceBeta() const { return UPPER_TIME_ACCEPTANCE_FACTOR ; }



//--------------------------------------------------------------------------
// Time primitives including single gaussian resolution
// These now interface to an external helper library
//

//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_SignalAlt_BaseClass::preCalculateTimeFactors() const
{
	if( useUpperTimeAcceptance() ) expL_stored = Mathematics::Exp_betaAcceptance( t, gamma_l(), resolution, upperTimeAcceptanceBeta() ) ;
	else expL_stored = Mathematics::Exp( t, gamma_l(), resolution ) ;
	if( useUpperTimeAcceptance() ) expH_stored = Mathematics::Exp_betaAcceptance( t, gamma_h(), resolution, upperTimeAcceptanceBeta() ) ;
	else expH_stored = Mathematics::Exp( t, gamma_h(), resolution ) ;
	expSin_stored = Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
	expCos_stored = Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
	return ;	
}


//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_SignalAlt_BaseClass::preCalculateTimeIntegrals() const
{
	if( useUpperTimeAcceptance() ) intExpL_stored = Mathematics::ExpInt_betaAcceptance( tlo, thi, gamma_l(), resolution, upperTimeAcceptanceBeta() )  ;
	else intExpL_stored = Mathematics::ExpInt( tlo, thi, gamma_l(), resolution )  ;
	if( useUpperTimeAcceptance() ) intExpH_stored = Mathematics::ExpInt_betaAcceptance( tlo, thi, gamma_h(), resolution, upperTimeAcceptanceBeta() )  ;
	else intExpH_stored = Mathematics::ExpInt( tlo, thi, gamma_h(), resolution )  ;
	intExpSin_stored = Mathematics::ExpSinInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
	intExpCos_stored = Mathematics::ExpCosInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
	return ;
}

//...................................................
//Exponentials 

double Bs2JpsiPhi_SignalAlt_BaseClass::expL() const 
{
	return expL_stored ; //Mathematics::Exp( t, gamma_l(), resolution ) ;
}

double Bs2JpsiPhi_SignalAlt_BaseClass::expH() const 
{
	return expH_stored ; //Mathematics::Exp( t, gamma_h(), resolution ) ;
}

double Bs2JpsiPhi_SignalAlt_BaseClass::intExpL( ) const {
	return intExpL_stored ;
	//if( useUpperTimeAcceptance() ) return Mathematics::ExpInt_betaAcceptance( tlo, thi, gamma_l(), resolution, upperTimeAcceptanceBeta() )  ;
	//else return Mathematics::ExpInt( tlo, thi, gamma_l(), resolution )  ;
}

double Bs2JpsiPhi_SignalAlt_BaseClass::intExpH( ) const {
	return intExpH_stored ;
	//if( useUpperTimeAcceptance() )return Mathematics::ExpInt_betaAcceptance( tlo, thi, gamma_h(), resolution, upperTimeAcceptanceBeta() )  ;
	//else return Mathematics::ExpInt( tlo, thi, gamma_h(), resolution )  ;
}


//......................................................
// Exponential x sine  and cosine

double Bs2JpsiPhi_SignalAlt_BaseClass::expSin() const  
{
    return expSin_stored ; // Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
}

double Bs2JpsiPhi_SignalAlt_BaseClass::expCos() const 
{
    return expCos_stored ; // Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
}

double Bs2JpsiPhi_SignalAlt_BaseClass::intExpSin( ) const 
{
	return intExpSin_stored ;   // Mathematics::ExpSinInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
}

// Integral of exp( - G * t ) * cos( dm * t )  
double Bs2JpsiPhi_SignalAlt_BaseClass::intExpCos( ) const 
{
	return intExpCos_stored ;   // Mathematics::ExpCosInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
}



//------------------------------------------------------------------------------
// These are the time factors and their analytic integrals for the one angle PDF

//..................................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorEven(  )  const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = 
	( 1.0 + cos(phi_s) ) * expL( ) 
	+ ( 1.0 - cos(phi_s) ) * expH( ) 
	+ q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) ;
	
	//DEBUG
	if( DEBUGFLAG && (result < 0) ) {
		cout << " Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorEven() : result < 0 " << endl ;
		cout << " ->term1 " << ( 1.0 + cos(phi_s) ) * expL( ) << endl ;
		cout << " ->term2 " << ( 1.0 - cos(phi_s) ) * expH( ) << endl ;
		cout << " ->term3 " << q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) << endl ;
		cout << "   -->sin(phis) "  << sin(phi_s) << endl ;
		cout << "   -->expSin    "  << expSin() << endl ;
		cout << "   -->tagFrac   "  << tagFraction << endl ;
		cout << "   -->delta_ms  "  << delta_ms << endl ;
	}
	return result ;
};

double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorEvenInt(  )  const
{

	double result = 
	( 1.0 + cos(phi_s) )  * intExpL()     
	+ ( 1.0 - cos(phi_s) )  * intExpH()          
	+ q() * ( 2.0 * sin(phi_s)   ) * intExpSin( ) * (1.0 - 2.0*tagFraction) ;
	return result ;
};


//..................................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorOdd(  )   const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = 
	( 1.0 - cos(phi_s) ) * expL( ) 
	+ ( 1.0 + cos(phi_s) ) * expH( ) 
	- q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) ;
	return result ;
};

double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorOddInt(  )  const
{
	double result = 
	( 1.0 - cos(phi_s) ) * intExpL()
	+ ( 1.0 + cos(phi_s) ) * intExpH() 
	- q() * ( 2.0 * sin(phi_s)   ) * intExpSin( ) * (1.0 - 2.0*tagFraction) ;
	return result ;
};


//----------------------------------------------------------
// These are the time factors and their analytic integrals for the three angle PDF

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorA0A0( )    const { return timeFactorEven( ) ; } ;      
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorA0A0Int( ) const { return timeFactorEvenInt( ) ; } ;

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorAPAP( )    const { return timeFactorEven( ) ; } ;
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorAPAPInt( ) const { return timeFactorEvenInt( ) ; } ;

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorATAT( )    const { return timeFactorOdd( ) ; } ;
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorATATInt( ) const { return timeFactorOddInt( ) ; } ;

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorImAPAT( ) const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = 
	q() * 2.0  * ( sin(delta1)*expCos( ) - cos(delta1)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( expH( ) - expL( ) ) * cos(delta1) * sin(phi_s)  ;
	
	return result ;
} ;

double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorImAPATInt( ) const
{
	double _tlo = tlo ;
	if(_tlo < 0.) _tlo = 0. ;
	
	double result = 
	q() * 2.0  * ( sin(delta1)*intExpCos() - cos(delta1)*cos(phi_s)*intExpSin() ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( intExpH() - intExpL() ) * cos(delta1) * sin(phi_s) ;	    
	return result ;
} ;


//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorReA0AP( )  const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = cos(delta2-delta1) * this->timeFactorEven(  ) ;
	return result ;
} ;

double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorReA0APInt( ) const
{
	double result = cos(delta2-delta1) * this->timeFactorEvenInt( ) ;
	return result ;
} ;


//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorImA0AT(  ) const
{
	//if( t < 0.0 ) return 0.0 ;
	double result =
	q() * 2.0  * ( sin(delta2)*expCos( ) - cos(delta2)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)	
	-1.0 * ( expH( ) - expL( ) ) * cos(delta2) * sin(phi_s) ;
	return result ;
} ;

double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorImA0ATInt( ) const
{
	double _tlo = tlo ;
	if(_tlo < 0.) _tlo = 0. ;
	
	double result = 
	q() * 2.0  * ( sin(delta2)*intExpCos() - cos(delta2)*cos(phi_s)*intExpSin()  ) * (1.0 - 2.0*tagFraction)
	-1.0 * ( intExpH() - intExpL()  ) * cos(delta2) * sin(phi_s) ;
	return result ;
} ;

//.... S wave additions.......

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorASAS( )    const { return timeFactorOdd( ) ; } ;
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorASASInt( ) const { return timeFactorOddInt( ) ; } ;


//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorReASAP( ) const
{
	//if( t < 0.0 ) return 0.0 ;
	
	double delta = delta_para - delta_s ;
	double result = 
	q() * 2.0  * ( cos(delta)*expCos( ) - sin(delta)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( expH( ) - expL( ) ) * sin(delta) * sin(phi_s)  ;
	
	return result ;
} ;

double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorReASAPInt( ) const
{
	double _tlo = tlo ;
	if(_tlo < 0.) _tlo = 0. ;

	double delta = delta_para - delta_s ;

	double result = 
	q() * 2.0  * ( cos(delta)*intExpCos() - sin(delta)*cos(phi_s)*intExpSin() ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( intExpH() - intExpL() ) * sin(delta) * sin(phi_s) ;	    
	return result ;
} ;


//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorImASAT( )  const
{
	//if( t < 0.0 ) return 0.0 ;
	double result = sin(delta_perp-delta_s) * this->timeFactorOdd(  ) ;
	return result ;
} ;

double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorImASATInt( ) const
{
	double result = sin(delta_perp-delta_s) * this->timeFactorOddInt( ) ;
	return result ;
} ;


//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorReASA0( ) const
{
	//if( t < 0.0 ) return 0.0 ;
	
	double delta = delta_zero - delta_s ;
	double result = 
	q() * 2.0  * ( cos(delta)*expCos( ) - sin(delta)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( expH( ) - expL( ) ) * sin(delta) * sin(phi_s)  ;
	
	return result ;
} ;

double Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorReASA0Int( ) const
{
	double _tlo = tlo ;
	if(_tlo < 0.) _tlo = 0. ;
	
	double delta = delta_zero - delta_s ;
	
	double result = 
	q() * 2.0  * ( cos(delta)*intExpCos() - sin(delta)*cos(phi_s)*intExpSin() ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( intExpH() - intExpL() ) * sin(delta) * sin(phi_s) ;	    
	return result ;
} ;


//------------------------------------------------------
// Angle factors for three angle PDFs

//........ P Wave ..........

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::angleFactorA0A0(  ) const
{
	// Normalised to  1	
	double result = 2.0 * ct1sq() * (1.0 - strsq()*cphsq() ) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::angleFactorAPAP(  ) const
{
	// Normalised to  1
	double result =  st1sq() * (1.0 - strsq()*sphsq() ) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::angleFactorATAT(  ) const
{
	// Normalised to  1
	double result = st1sq() * strsq() * (9.0/32.0/TMath::Pi());
	return result ;
	
};

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::angleFactorImAPAT(  ) const
{
	// Normalised to  0
	double theta_tr = acos(ctheta_tr) ;		
	double result =   -1.0 *  st1sq() * sin(2.0*theta_tr) * sin(phi_tr) * (9.0/32.0/TMath::Pi()) ;
	return result ;	
};

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::angleFactorReA0AP( ) const
{
	// Normalised to  0
	double theta_1 = acos(ctheta_1) ;	
	double result =    sin(2.0*theta_1) * strsq() * sin(2.0*phi_tr) / sqrt(2.0) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::angleFactorImA0AT(  ) const
{
	// Normalised to  0
	double theta_tr = acos(ctheta_tr) ;		
	double theta_1 = acos(ctheta_1) ;		
	double result =  +1.0*   sin(2.0*theta_1) * sin(2.0*theta_tr) * cos(phi_tr) / sqrt(2.0) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//......  S wave additions ....

//.............................
double Bs2JpsiPhi_SignalAlt_BaseClass::angleFactorASAS(  ) const
{
	double result =  (1.0 - strsq()*cphsq() ) * (2./3.) * (9.0/32.0/TMath::Pi());
	return result ;	
};

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::angleFactorReASAP(  ) const
{
	double stheta_1 =  sqrt(st1sq());		
	double result =   strsq() * stheta_1 * sin(2.0*phi_tr) * (sqrt(6.)/3.) * (9.0/32.0/TMath::Pi()) ;
	return result ;	
};

//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::angleFactorImASAT(  ) const
{
	double theta_tr = acos(ctheta_tr) ;		
	double stheta_1 =  sqrt(st1sq());		
	double result = -1.0 *  sin(2.0*theta_tr) * stheta_1 * cos(phi_tr) * (sqrt(6.)/3.) * (9.0/32.0/TMath::Pi()) ;
	return result ;
};


//...........................
double Bs2JpsiPhi_SignalAlt_BaseClass::angleFactorReASA0(  ) const
{
	double result = -1.0 *  ( 1.0 -  strsq()* cphsq() ) * ctheta_1 *  (4.0*sqrt(3.)/3.) * (9.0/32.0/TMath::Pi()) ;
	return result ;	
};
	

