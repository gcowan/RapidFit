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
//Constructor(s)

//...............................
// Old default constructor
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
	// PELC NEW additions for v2
	, cosphisName( make_pair("cosphis",-1) )
	, sinphisName( make_pair("sinphis",-1) )
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
		
	//PELC  - debug to plot the distribution of PDF values for each event 
	//histOfPdfValues = new TH1D( "HistOfPdfValue" ,  "HistOfPdfValue" , 110, -0.00001, 0.00001 ) ;
	//c0  = new TCanvas;
	//histCounter = 0;
	//~PELC
}


//.....................................
// New Constructor which takes configuration object
Bs2JpsiPhi_SignalAlt_BaseClass::Bs2JpsiPhi_SignalAlt_BaseClass(PDFConfigurator configurator ) : 
// Physics parameters
	  gammaName     ( make_pair(configurator.getName("gamma"),-1) )
	, deltaGammaName( make_pair(configurator.getName("deltaGamma"),-1) )
	, deltaMName    ( make_pair(configurator.getName("deltaM"),-1))
	, Phi_sName     ( make_pair(configurator.getName("Phi_s"),-1))
	, Azero_sqName  ( make_pair(configurator.getName("Azero_sq"),-1) )
	, Aperp_sqName  ( make_pair(configurator.getName("Aperp_sq"),-1) )
	, As_sqName		( make_pair(configurator.getName("As_sq"),-1) )
	, delta_zeroName( make_pair(configurator.getName("delta_zero"),-1) )
	, delta_paraName( make_pair(configurator.getName("delta_para"),-1) )
	, delta_perpName( make_pair(configurator.getName("delta_perp"),-1) )
	, delta_sName	( make_pair(configurator.getName("delta_s"),-1) )
	// PELC NEW additions for v2
	, cosphisName( make_pair(configurator.getName("cosphis"),-1) )
	, sinphisName( make_pair(configurator.getName("sinphis"),-1) )
	// Detector parameters
	, mistagName		( make_pair(configurator.getName("mistag"),-1) )
	, res1Name			( make_pair(configurator.getName("timeResolution1"),-1) )
	, res2Name			( make_pair(configurator.getName("timeResolution2"),-1) )
	, res1FractionName	( make_pair(configurator.getName("timeResolution1Fraction"),-1) )
	, timeOffsetName	( make_pair(configurator.getName("timeOffset"),-1) )
	// Angular acceptance factors
	, angAccI1Name ( make_pair(configurator.getName("angAccI1"),-1) )
	, angAccI2Name ( make_pair(configurator.getName("angAccI2"),-1) )
	, angAccI3Name ( make_pair(configurator.getName("angAccI3"),-1) )
	, angAccI4Name ( make_pair(configurator.getName("angAccI4"),-1) )
	, angAccI5Name ( make_pair(configurator.getName("angAccI5"),-1) )
	, angAccI6Name ( make_pair(configurator.getName("angAccI6"),-1) )
	, angAccI7Name ( make_pair(configurator.getName("angAccI7"),-1) )
	, angAccI8Name ( make_pair(configurator.getName("angAccI8"),-1) )
	, angAccI9Name ( make_pair(configurator.getName("angAccI9"),-1) )
	, angAccI10Name( make_pair(configurator.getName("angAccI10"),-1) )
	// Observables
	, timeName      ( make_pair(configurator.getName("time"),-1) )
	, cosThetaName	( make_pair(configurator.getName("cosTheta"),-1) )
	, phiName	    ( make_pair(configurator.getName("phi"),-1) )
	, cosPsiName	( make_pair(configurator.getName("cosPsi"),-1) )
	, tagName	    ( make_pair(configurator.getName("tag"),-1) )
	, timeAcceptanceCategoryName ( make_pair(configurator.getName("timeAcceptanceCategory"),-1) )
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
	
	//PELC  - debug to plot the distribution of PDF values for each event 
	//histOfPdfValues = new TH1D( "HistOfPdfValue" ,  "HistOfPdfValue" , 110, -0.00001, 0.00001 ) ;
	//c0  = new TCanvas;
	//histCounter = 0;
	//~PELC
}


//........................................................
//Destructor
Bs2JpsiPhi_SignalAlt_BaseClass::~Bs2JpsiPhi_SignalAlt_BaseClass()
{
	cout << "hello from BaseClass PDF" << endl;
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


//-------------------------------------------------------------------------
// Differential crtoss sections and normalisations
//

//...................................
// Main Diff cross section

double Bs2JpsiPhi_SignalAlt_BaseClass::diffXsec(  )  const
{   
	preCalculateTimeFactors() ;
	
	double xsec = 
	
	0.5 * A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) +
	0.5 * AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) +
	0.5 * AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) +
	
	0.5 * AP()*AT() * timeFactorImAPAT(  ) * angleFactorImAPAT( ) +
	0.5 * A0()*AP() * timeFactorReA0AP(  ) * angleFactorReA0AP( ) +
	0.5 * A0()*AT() * timeFactorImA0AT(  ) * angleFactorImA0AT( ) +
	
	0.5 * AS()*AS() * timeFactorASAS(  ) * angleFactorASAS( ) +
	
	0.5 * AS()*AP() * timeFactorReASAP(  ) * angleFactorReASAP( ) +
	0.5 * AS()*AT() * timeFactorImASAT(  ) * angleFactorImASAT( ) +
	0.5 * AS()*A0() * timeFactorReASA0(  ) * angleFactorReASA0( ) ;
	
	//PELC DEBUG 
	if( DEBUGFLAG && (xsec < 0) ) {
		cout << " Bs2JpsiPhi_SignalAlt_MO_v1::diffXsec( ) : return value < 0 = " << xsec << endl ;
		cout << "   A0()*A0() term: " <<  A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) << endl ;
		cout << "   AP()*AP() term: " <<AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) << endl ;
		cout << "   AT()*AT() term: " <<AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) << endl << endl ;
		
		cout << "   AP()*AT() term: " <<AP()*AT() * timeFactorImAPAT(  ) * angleFactorImAPAT( )<< endl ;
		cout << "                 : " <<AP()*AT() <<" / "<<  timeFactorImAPAT( )  <<" / "<<  angleFactorImAPAT( )<< endl ;
		cout << "   A0()*AP() term: " <<A0()*AP() * timeFactorReA0AP(  ) * angleFactorReA0AP( )<< endl ;
		cout << "                 : " <<A0()*AP() <<" / "<<  timeFactorReA0AP(  ) <<" / "<<  angleFactorReA0AP( )<< endl ;
		cout << "   A0()*AT() term: " <<A0()*AT() * timeFactorImA0AT(  ) * angleFactorImA0AT( )<< endl << endl;
		cout << "                 : " <<A0()*AT() <<" / "<<  timeFactorImA0AT(  ) <<" / "<<  angleFactorImA0AT( )<< endl << endl;
		
		cout << "   AS()*AS() term: " <<AS()*AS() * timeFactorASAS(  ) * angleFactorASAS( ) << endl << endl ;

		cout << "   AS()*AP() term: " <<AS()*AP() * timeFactorReASAP(  ) * angleFactorReASAP( )<< endl ;
		cout << "                 : " <<AS()*AP() <<" / "<<   timeFactorReASAP(  ) <<" / "<<  angleFactorReASAP( )<< endl ;
		cout << "   AS()*AT() term: " <<AS()*AT() * timeFactorImASAT(  ) * angleFactorImASAT( )<< endl ;
		cout << "                 : " <<AS()*AT() <<" / "<<   timeFactorImASAT(  ) <<" / "<<   angleFactorImASAT( )<< endl ;
		cout << "   AS()*A0() term: " <<AS()*A0() * timeFactorReASA0(  ) * angleFactorReASA0( )<< endl ;
		cout << "                 : " <<AS()*A0() <<" / "<<   timeFactorReASA0(  ) <<" / "<<  angleFactorReASA0( )<< endl << endl ;

		double PwaveTot =
		0.5 * A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) +
		0.5 * AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) +
		0.5 * AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) +		
		0.5 * AP()*AT() * timeFactorImAPAT(  ) * angleFactorImAPAT( ) +
		0.5 * A0()*AP() * timeFactorReA0AP(  ) * angleFactorReA0AP( ) +
		0.5 * A0()*AT() * timeFactorImA0AT(  ) * angleFactorImA0AT( ) ;
		
		double SwaveAdditions =
		0.5 * AS()*AS() * timeFactorASAS(  ) * angleFactorASAS( ) +
		0.5 * AS()*AP() * timeFactorReASAP(  ) * angleFactorReASAP( ) +
		0.5 * AS()*AT() * timeFactorImASAT(  ) * angleFactorImASAT( ) +
		0.5 * AS()*A0() * timeFactorReASA0(  ) * angleFactorReASA0( ) ;

		cout << "   Pwave Only : " << PwaveTot << endl ;
		cout << "   Swave add : " <<  SwaveAdditions << endl ;

	}
	
	//PELC - This turned out to be an important debugging tool 
	//switch it on to see the values of PDF being returend.  If ANY go negative, it means there is a sign wrong in one or more of the terms
	//You need to enable in the .h file as well
	//histOfPdfValues->Fill(xsec) ;	
	//histCounter++ ;
	//if( histCounter > 10000 ) {
	//	histOfPdfValues->Draw() ;
	//	c0->Update() ;	
	//	c0->SaveAs( "histOfPdfValues-from-Evaluate.eps" ) ;
	//	histCounter = 0 ;
	//}
	 
		
	
	return xsec ;
}

//...................................
// Integral over all variables: t + angles

double Bs2JpsiPhi_SignalAlt_BaseClass::diffXsecNorm1(  ) const
{ 
	preCalculateTimeIntegrals() ;
	
	double norm =
	
	0.5 * A0()*A0() * timeFactorA0A0Int(  ) * angAccI1   +  
	0.5 * AP()*AP() * timeFactorAPAPInt(  ) * angAccI2   +  
	0.5 * AT()*AT() * timeFactorATATInt(  ) * angAccI3   +  
	
	0.5 * AP()*AT() * timeFactorImAPATInt(  ) * angAccI4 +  
	0.5 * A0()*AP() * timeFactorReA0APInt(  ) * angAccI5 +  
	0.5 * A0()*AT() * timeFactorImA0ATInt(  ) * angAccI6 +  
	
	0.5 * AS()*AS() * timeFactorASASInt(  ) * angAccI7   +  
	
	0.5 * AS()*AP() * timeFactorReASAPInt(  ) * angAccI8 +  
	0.5 * AS()*AT() * timeFactorImASATInt(  ) * angAccI9 +  
	0.5 * AS()*A0() * timeFactorReASA0Int(  ) * angAccI10 ;  
	
	//PELC DEBUG 	
	 if( DEBUGFLAG && (norm < 0) ) {
	 cout << endl ;
	 cout <<  A0()*A0() * timeFactorA0A0Int(  )* angAccI1  << endl ;
	 cout <<  AP()*AP() * timeFactorAPAPInt(  )* angAccI2 << endl ;
	 cout <<  AT()*AT() * timeFactorATATInt(  )* angAccI3 << endl << endl ;
	 cout <<  AP()*AT() * timeFactorImAPATInt(  ) * angAccI4<< endl ;
	 cout <<  A0()*AP() * timeFactorReA0APInt(  ) * angAccI5<< endl ;
	 cout <<  A0()*AT() * timeFactorImA0ATInt(  ) * angAccI6<< endl << endl;
	 cout <<  AS()*AS() * timeFactorASASInt(  ) * angAccI7 << endl ;
	 cout <<  AS()*AP() * timeFactorReASAPInt(  ) * angAccI8<< endl ;
	 cout <<  AS()*AT() * timeFactorImASATInt(  ) * angAccI9<< endl ;
	 cout <<  AS()*A0() * timeFactorReASA0Int(  ) * angAccI10<< endl ;
	 }
	 
	return norm ;
}



//...................................
// Integral over angles only for a fixed time.

double Bs2JpsiPhi_SignalAlt_BaseClass::diffXsecNorm2(  ) const
{          
	preCalculateTimeFactors() ;
	
	double norm = 
	
	0.5 * A0()*A0() * timeFactorA0A0(  ) * angAccI1 +
	0.5 * AP()*AP() * timeFactorAPAP(  ) * angAccI2 +
	0.5 * AT()*AT() * timeFactorATAT(  ) * angAccI3 +
	
	0.5 * AP()*AT() * timeFactorImAPAT(  ) * angAccI4 +
	0.5 * A0()*AP() * timeFactorReA0AP(  ) * angAccI5 +
	0.5 * A0()*AT() * timeFactorImA0AT(  ) * angAccI6 +
	
	0.5 * AS()*AS() * timeFactorASAS(  ) * angAccI7 +
	
	0.5 * AS()*AP() * timeFactorReASAP(  ) * angAccI8 +
	0.5 * AS()*AT() * timeFactorImASAT(  ) * angAccI9 +
	0.5 * AS()*A0() * timeFactorReASA0(  ) * angAccI10 ;
	
	return norm ;
}


//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,

double Bs2JpsiPhi_SignalAlt_BaseClass::diffXsecCompositeNorm1(  )  
{   
	if( useLowerTimeAcceptance() ) 
	{
		double norm = 0 ;
		double tlo_remember = tlo ;
		
		timeAcceptance.configure( tlo ) ;
		
		bool firstBin = true ;
		for( int islice = timeAcceptance.firstSlice( ); islice <= timeAcceptance.lastSlice( ); ++islice ) 
		{
			if( firstBin )firstBin = false ;
			else tlo = timeAcceptance.sliceStart( islice ) ;
			norm += this->diffXsecNorm1(  ) * timeAcceptance.fraction( islice ) ;
		}
		
		tlo =  tlo_remember ;
		return norm ;	
	}
	
	else 
	{
		return this->diffXsecNorm1() ;
	}
	
}



