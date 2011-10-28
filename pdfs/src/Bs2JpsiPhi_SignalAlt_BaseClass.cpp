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
Bs2JpsiPhi_SignalAlt_BaseClass::Bs2JpsiPhi_SignalAlt_BaseClass() : useCosAndSin(), allowNegativeAsSq(),
	// Physics parameters
	  gammaName     		( "gamma" )
	, deltaGammaName		( "deltaGamma" )
	, deltaMName    		( "deltaM" )
	, Phi_sName     		( "Phi_s" )
	, Azero_sqName  		( "Azero_sq" )
	, Apara_sqName  		( "Apara_sq" )
	, Aperp_sqName  		( "Aperp_sq" )
	, As_sqName			( "As_sq" )
	, delta_zeroName		( "delta_zero" )
	, delta_paraName		( "delta_para" )
	, delta_perpName		( "delta_perp" )
	, delta_sName			( "delta_s" )
	// PELC NEW additions for v2
	, cosphisName			( "cosphis" )
	, sinphisName			( "sinphis" )
	// Mistag parameters
	, mistagName			( "mistag" )
	, mistagP1Name			( "mistagP1" )
	, mistagP0Name			( "mistagP0" )
	, mistagSetPointName		( "mistagSetPoint" )
	// Detector Parameters
	, res1Name			( "timeResolution1" )
	, res2Name			( "timeResolution2" )
	, res1FractionName		( "timeResolution1Fraction" )
	, timeOffsetName		( "timeOffset" )
	// Angular acceptance factors
	, angAccI1Name			( "angAccI1" )
	, angAccI2Name			( "angAccI2" )
	, angAccI3Name			( "angAccI3" )
	, angAccI4Name			( "angAccI4" )
	, angAccI5Name			( "angAccI5" )
	, angAccI6Name			( "angAccI6" )
	, angAccI7Name			( "angAccI7" )
	, angAccI8Name			( "angAccI8" )
	, angAccI9Name			( "angAccI9" )
	, angAccI10Name 		( "angAccI10" )
	// Observables
	, timeName			( "time" )
	, cosThetaName			( "cosTheta" )
	, phiName			( "phi" )
	, cosPsiName			( "cosPsi" )
	, tagName			( "tag" )
	, timeAcceptanceCategoryName	( "timeAcceptanceCategory" )
	, timeConstraintName		( "time" )
	// Other things
	//objects
	,t(), ctheta_tr(), phi_tr(), ctheta_1(), tag(), timeAcceptanceCategory(), _gamma(), dgam(), Aperp_sq(), Apara_sq(), Azero_sq(), As_sq(), delta_para(),
	delta_perp(), delta_zero(), delta_s(), delta1(), delta2(), delta_ms(), phi_s(), _cosphis(), _sinphis(), _mistag(), _mistagP1(), _mistagP0(), _mistagSetPoint(),
	resolution(), resolution1(), resolution2(), resolution1Fraction(), timeOffset(), angAccI1(), angAccI2(), angAccI3(), angAccI4(), angAccI5(), angAccI6(),
	angAccI7(), angAccI8(), angAccI9(), angAccI10(), tlo(), thi(), timeAcceptance(), expL_stored(), expH_stored(), expSin_stored(), expCos_stored(),
	intExpL_stored(), intExpH_stored(), intExpSin_stored(), intExpCos_stored()
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

//	Copy Constructor
Bs2JpsiPhi_SignalAlt_BaseClass::Bs2JpsiPhi_SignalAlt_BaseClass( const Bs2JpsiPhi_SignalAlt_BaseClass& input ) :
	useCosAndSin(input.useCosAndSin), allowNegativeAsSq(input.allowNegativeAsSq), gammaName(input.gammaName), deltaGammaName(input.deltaGammaName), deltaMName(input.deltaMName),
	Phi_sName(input.Phi_sName), Azero_sqName(input.Azero_sqName), Apara_sqName(input.Apara_sqName), Aperp_sqName(input.Aperp_sqName), As_sqName(input.As_sqName),
	delta_zeroName(input.delta_zeroName), delta_paraName(input.delta_paraName), delta_perpName(input.delta_perpName), delta_sName(input.delta_sName), cosphisName(input.cosphisName),
	sinphisName(input.sinphisName), mistagName(input.mistagName), mistagP1Name(input.mistagP1Name), mistagP0Name(input.mistagP0Name), mistagSetPointName(input.mistagSetPointName),
	res1Name(input.res1Name), res2Name(input.res2Name), res1FractionName(input.res1FractionName), timeOffsetName(input.timeOffsetName), angAccI1Name(input.angAccI1Name),
	angAccI2Name(input.angAccI2Name), angAccI3Name(input.angAccI3Name), angAccI4Name(input.angAccI4Name), angAccI5Name(input.angAccI5Name), angAccI6Name(input.angAccI6Name),
	angAccI7Name(input.angAccI7Name), angAccI8Name(input.angAccI8Name), angAccI9Name(input.angAccI9Name), angAccI10Name(input.angAccI10Name), timeName(input.timeName),
	cosThetaName(input.cosThetaName), phiName(input.phiName), cosPsiName(input.cosPsiName), tagName(input.tagName), timeAcceptanceCategoryName(input.timeAcceptanceCategoryName),
	timeConstraintName(input.timeConstraintName), t(input.t), ctheta_tr(input.ctheta_tr), phi_tr(input.phi_tr), ctheta_1(input.ctheta_1), tag(input.tag),
	timeAcceptanceCategory(input.timeAcceptanceCategory), _gamma(input._gamma), dgam(input.dgam), Aperp_sq(input.Aperp_sq), Apara_sq(input.Apara_sq), Azero_sq(input.Azero_sq),
	As_sq(input.As_sq), delta_para(input.delta_para), delta_perp(input.delta_perp), delta_zero(input.delta_zero), delta_s(input.delta_s), delta1(input.delta1), delta2(input.delta2),
	delta_ms(input.delta_ms), phi_s(input.phi_s), _cosphis(input._cosphis), _sinphis(input._sinphis), _mistag(input._mistag), _mistagP1(input._mistagP1), _mistagP0(input._mistagP0),
	_mistagSetPoint(input._mistagSetPoint), resolution(input.resolution), resolution1(input.resolution1), resolution2(input.resolution2), resolution1Fraction(input.resolution1Fraction),
	timeOffset(input.timeOffset), angAccI1(input.angAccI1), angAccI2(input.angAccI2), angAccI3(input.angAccI3), angAccI4(input.angAccI4), angAccI5(input.angAccI5), angAccI6(input.angAccI6),
	angAccI7(input.angAccI7), angAccI8(input.angAccI8), angAccI9(input.angAccI9), angAccI10(input.angAccI10), tlo(input.tlo), thi(input.thi), timeAcceptance( *( new TimeAcceptanceFunction() )),
	expL_stored(input.expL_stored), expH_stored(input.expH_stored), expSin_stored(input.expSin_stored), expCos_stored(input.expCos_stored), intExpL_stored(input.intExpL_stored),
	intExpH_stored(input.intExpH_stored), intExpSin_stored(input.intExpSin_stored), intExpCos_stored(input.intExpCos_stored)
{
}

//.....................................
// New Constructor which takes configuration object
Bs2JpsiPhi_SignalAlt_BaseClass::Bs2JpsiPhi_SignalAlt_BaseClass(PDFConfigurator* configurator ) : useCosAndSin(), allowNegativeAsSq(),
// Physics parameters
	  gammaName     		( configurator->getName("gamma") )
	, deltaGammaName		( configurator->getName("deltaGamma") )
	, deltaMName			( configurator->getName("deltaM") )
	, Phi_sName			( configurator->getName("Phi_s") )
	, Azero_sqName			( configurator->getName("Azero_sq") )
	, Apara_sqName			( configurator->getName("Apara_sq") )
	, Aperp_sqName			( configurator->getName("Aperp_sq") )
	, As_sqName			( configurator->getName("As_sq") )
	, delta_zeroName		( configurator->getName("delta_zero") )
	, delta_paraName		( configurator->getName("delta_para") )
	, delta_perpName		( configurator->getName("delta_perp") )
	, delta_sName			( configurator->getName("delta_s") )
	// PELC NEW additions for v2
	, cosphisName			( configurator->getName("cosphis") )
	, sinphisName			( configurator->getName("sinphis") )
	// Detector parameters
	, mistagName			( configurator->getName("mistag") )
	, mistagP1Name			( configurator->getName("mistagP1") )
	, mistagP0Name			( configurator->getName("mistagP0") )
	, mistagSetPointName		( configurator->getName("mistagSetPoint") )
	// Detector parameters
	, res1Name			( configurator->getName("timeResolution1") )
	, res2Name			( configurator->getName("timeResolution2") )
	, res1FractionName		( configurator->getName("timeResolution1Fraction") )
	, timeOffsetName		( configurator->getName("timeOffset") )
	// Angular acceptance factors
	, angAccI1Name			( configurator->getName("angAccI1") )
	, angAccI2Name			( configurator->getName("angAccI2") )
	, angAccI3Name			( configurator->getName("angAccI3") )
	, angAccI4Name			( configurator->getName("angAccI4") )
	, angAccI5Name			( configurator->getName("angAccI5") )
	, angAccI6Name			( configurator->getName("angAccI6") )
	, angAccI7Name			( configurator->getName("angAccI7") )
	, angAccI8Name			( configurator->getName("angAccI8") )
	, angAccI9Name			( configurator->getName("angAccI9") )
	, angAccI10Name			( configurator->getName("angAccI10") )
	// Observables
	, timeName			( configurator->getName("time") )
	, cosThetaName			( configurator->getName("cosTheta") )
	, phiName			( configurator->getName("phi") )
	, cosPsiName			( configurator->getName("cosPsi") )
	, tagName			( configurator->getName("tag") )
	, timeAcceptanceCategoryName	( configurator->getName("timeAcceptanceCategory") )
	, timeConstraintName		( configurator->getName("time") )
	// Other things
	//objects
	,t(), ctheta_tr(), phi_tr(), ctheta_1(), tag(), timeAcceptanceCategory(), _gamma(), dgam(), Aperp_sq(), Apara_sq(), Azero_sq(), As_sq(), delta_para(),
	delta_perp(), delta_zero(), delta_s(), delta1(), delta2(), delta_ms(), phi_s(), _cosphis(), _sinphis(), _mistag(), _mistagP1(), _mistagP0(), _mistagSetPoint(),
	resolution(), resolution1(), resolution2(), resolution1Fraction(), timeOffset(), angAccI1(), angAccI2(), angAccI3(), angAccI4(), angAccI5(), angAccI6(),
	angAccI7(), angAccI8(), angAccI9(), angAccI10(), tlo(), thi(), timeAcceptance(), expL_stored(), expH_stored(), expSin_stored(), expCos_stored(),
	intExpL_stored(), intExpH_stored(), intExpSin_stored(), intExpCos_stored()
{

	useCosAndSin = configurator->isTrue( "UseCosAndSin" ) ;
	allowNegativeAsSq = configurator->isTrue( "AllowNegativeAsSq" ) ;
		
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
	//cout << "hello from BaseClass PDF" << endl;
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
		cout << " Bs2JpsiPhi_SignalAlt_BaseClass_v1::diffXsec( ) : return value < 0 = " << xsec << endl ;
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



