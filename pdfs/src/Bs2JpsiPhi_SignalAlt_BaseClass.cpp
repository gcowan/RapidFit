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
	, As_sqName	( "As_sq" )
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

