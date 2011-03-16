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
	  gammaName     ( make_pair("gamma",-1) )
	, deltaGammaName( make_pair("deltaGamma",-1) )
	, deltaMName    ( make_pair("deltaM",-1))
	, Phi_sName     ( make_pair("Phi_s",-1))
	, Azero_sqName  ( make_pair("Azero_sq",-1) )
	, Aperp_sqName  ( make_pair("Aperp_sq",-1) )
	, As_sqName	( make_pair("As_sq",-1) )
	, delta_zeroName( make_pair("delta_zero",-1) )
	, delta_paraName( make_pair("delta_para",-1) )
	, delta_perpName( make_pair("delta_perp",-1) )
	, delta_sName( make_pair("delta_s",-1) )
	// Detector parameters
	, mistagName	( make_pair("mistag",-1) )
	, res1Name	( make_pair("timeResolution1",-1) )
	, res2Name	( make_pair("timeResolution2",-1) )
	, res1FractionName	( make_pair("timeResolution1Fraction",-1) )
	, timeOffsetName	( make_pair("timeOffset",-1) )
	// Angular acceptance factors
	, angAccI1Name ( make_pair("angAccI1",-1) )
	, angAccI2Name ( make_pair("angAccI2",-1) )
	, angAccI3Name ( make_pair("angAccI3",-1) )
	, angAccI4Name ( make_pair("angAccI4",-1) )
	, angAccI5Name ( make_pair("angAccI5",-1) )
	, angAccI6Name ( make_pair("angAccI6",-1) )
	, angAccI7Name ( make_pair("angAccI7",-1) )
	, angAccI8Name ( make_pair("angAccI8",-1) )
	, angAccI9Name ( make_pair("angAccI9",-1) )
	, angAccI10Name ( make_pair("angAccI10",-1) )
	// Observables
	, timeName	    ( make_pair("time",-1) )
	, cosThetaName	( make_pair("cosTheta",-1) )
	, phiName	    ( make_pair("phi",-1) )
	, cosPsiName	( make_pair("cosPsi",-1) )
	, tagName	    ( make_pair("tag",-1) )
	, timeAcceptanceCategoryName ( make_pair("timeAcceptanceCategory",-1) )
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

