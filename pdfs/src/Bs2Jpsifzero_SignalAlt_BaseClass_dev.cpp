// $Id: Bs2Jpsifzero_SignalAlt_BaseClass_dev.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2Jpsifzero_SignalAlt_BaseClass_dev Bs2Jpsifzero_SignalAlt_BaseClass_dev.cpp
 *
 *  Base class for Bs2Jpsifzero_SignalAlt..... PDFs
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-12
 */

//	ROOT Headers
#include "TMath.h"
//	RapidFit Headers
#include "Bs2Jpsifzero_SignalAlt_BaseClass_dev.h"
#include "SlicedAcceptance.h"
//	System Headers
#include <iostream>
#include <math.h>

PDF_CREATOR( Bs2Jpsifzero_SignalAlt_BaseClass_dev );

//......................................
//Constructor(s)

//.....................................
// New Constructor which takes configuration object
Bs2Jpsifzero_SignalAlt_BaseClass_dev::Bs2Jpsifzero_SignalAlt_BaseClass_dev(PDFConfigurator* configurator ) : 
// Physics parameters
	  gammaName     		( configurator->getName("gamma") )
	, deltaGammaName		( configurator->getName("deltaGamma") )
	, deltaMName			( configurator->getName("deltaM") )
	, Phi_sName				( configurator->getName("Phi_s") )
	, Azero_sqName			( configurator->getName("Azero_sq") )
	, Apara_sqName			( configurator->getName("Apara_sq") )
	, Aperp_sqName			( configurator->getName("Aperp_sq") )
	, As_sqName				( configurator->getName("As_sq") )
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
	, mistagSetPointName	( configurator->getName("mistagSetPoint") )
	// Detector parameters
	, resScaleName			( configurator->getName("timeResolutionScale") )
	, res1Name				( configurator->getName("timeResolution1") )
	, res2Name				( configurator->getName("timeResolution2") )
	, res3Name				( configurator->getName("timeResolution3") )
	, res2FractionName		( configurator->getName("timeResolution2Fraction") )
	, res3FractionName		( configurator->getName("timeResolution3Fraction") )
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
	, timeName				( configurator->getName("time") )
	, cosThetaName			( configurator->getName("cosTheta") )
	, phiName				( configurator->getName("phi") )
	, cosPsiName			( configurator->getName("cosPsi") )
	, tagName				( configurator->getName("tag") )
	//X, timeAcceptanceCategoryName	( configurator->getName("timeAcceptanceCategory") )
	, timeConstraintName	( configurator->getName("time") )
	// Other things
	, _useTimeAcceptance(false)
	, _numericIntegralForce(false)
	, _numericIntegralTimeOnly(false)
	, _useCosAndSin(false) 
	, allowNegativeAsSq(false)
	//objects
	,t(), ctheta_tr(), phi_tr(), ctheta_1(), tag(), /*timeAcceptanceCategory(),*/ _gamma(), dgam(), Aperp_sq(), Apara_sq(), Azero_sq(), As_sq(), delta_para(),
	delta_perp(), delta_zero(), delta_s(), delta1(), delta2(), delta_ms(), phi_s(), _cosphis(), _sinphis(), _mistag(), _mistagP1(), _mistagP0(), _mistagSetPoint(),
	resolutionScale(), resolution1(), resolution2(), resolution3(), resolution2Fraction(), resolution3Fraction(), timeOffset(), 
	angAccI1(), angAccI2(), angAccI3(), angAccI4(), angAccI5(), angAccI6(), angAccI7(), angAccI8(), angAccI9(), angAccI10(), 
	tlo(), thi(), expL_stored(), expH_stored(), expSin_stored(), expCos_stored(),
	intExpL_stored(), intExpH_stored(), intExpSin_stored(), intExpCos_stored(), timeAcc(NULL), normalisationCacheValid(false), resolution()
{

	//...........................................
	// Configure to use time acceptance machinery 
	_useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;
	
	if( useTimeAcceptance() ) {
		if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) ) {
			timeAcc = new SlicedAcceptance( 0., 14.0, 0.0157 ) ;
			cout << "Bs2Jpsifzero_SignalAlt_BaseClass_dev:: Constructing timeAcc: Upper time acceptance beta=0.0157 [0 < t < 14] " << endl ;
		}
		else if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Lower2010" ) ) {
			timeAcc = new SlicedAcceptance( "Lower2010" ) ;
			cout << "Bs2Jpsifzero_SignalAlt_BaseClass_dev:: Constructing timeAcc: Lower time acceptance 2010  [0 < t < 14] " << endl ;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" ) {
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ) ) ;
			cout << "Bs2Jpsifzero_SignalAlt_BaseClass_dev:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl ;
		}
	}
	else {
			timeAcc = new SlicedAcceptance( 0., 14. ) ;
			cout << "Bs2Jpsifzero_SignalAlt_BaseClass_dev:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl ;
	}
	
	
	//...........................................
	// Configure numerical integration options 
	_numericIntegralForce    = configurator->isTrue( "NumericIntegralForce") ;
	_numericIntegralTimeOnly = configurator->isTrue( "NumericIntegralTimeOnly" ) ;
	
	//...........................................
	// Configure fit parameter options 
	_useCosAndSin = configurator->isTrue( "UseCosAndSin" ) ;
	allowNegativeAsSq = configurator->isTrue( "AllowNegativeAsSq" ) ;
	

	
	//PELC  - debug to plot the distribution of PDF values for each event 
	//histOfPdfValues = new TH1D( "HistOfPdfValue" ,  "HistOfPdfValue" , 110, -0.00001, 0.00001 ) ;
	//c0  = new TCanvas;
	//histCounter = 0;
	//~PELC
}


Bs2Jpsifzero_SignalAlt_BaseClass_dev::Bs2Jpsifzero_SignalAlt_BaseClass_dev( const Bs2Jpsifzero_SignalAlt_BaseClass_dev& input )	:
	gammaName(input.gammaName), deltaGammaName(input.deltaGammaName), deltaMName(input.deltaMName), Phi_sName( input.Phi_sName ),
	Azero_sqName(input.Azero_sqName), Apara_sqName(input.Apara_sqName), Aperp_sqName(input.Aperp_sqName), As_sqName(input.As_sqName),
	delta_zeroName(input.delta_zeroName), delta_paraName(input.delta_paraName), delta_perpName(input.delta_perpName),
	delta_sName(input.delta_sName), cosphisName(input.cosphisName), sinphisName(input.sinphisName), mistagName(input.mistagName),
	mistagP1Name(input.mistagP1Name), mistagP0Name(input.mistagP0Name), mistagSetPointName(input.mistagSetPointName),
	resScaleName(input.resScaleName), res1Name(input.res1Name), res2Name(input.res2Name), res3Name(input.res3Name),
	res2FractionName(input.res2FractionName), res3FractionName(input.res3FractionName), timeOffsetName(input.timeOffsetName),
	angAccI1Name(input.angAccI1Name), angAccI2Name(input.angAccI2Name), angAccI3Name(input.angAccI3Name), angAccI4Name(input.angAccI4Name),
	angAccI5Name(input.angAccI5Name), angAccI6Name(input.angAccI6Name), angAccI7Name(input.angAccI7Name), angAccI8Name(input.angAccI8Name),
	angAccI9Name(input.angAccI9Name), angAccI10Name(input.angAccI10Name), timeName(input.timeName), cosThetaName(input.cosThetaName),
	phiName(input.phiName), cosPsiName(input.cosPsiName), tagName(input.tagName), timeConstraintName(input.timeConstraintName), t(input.t),
	ctheta_tr(input.ctheta_tr), phi_tr(input.phi_tr), ctheta_1(input.ctheta_1), tag(input.tag), _gamma(input._gamma), dgam(input.dgam),
	Aperp_sq(input.Aperp_sq), Apara_sq(input.Apara_sq), Azero_sq(input.Azero_sq), As_sq(input.As_sq), delta_para(input.delta_para),
	delta_perp(input.delta_perp), delta_zero(input.delta_zero), delta_s(input.delta_s), delta1(input.delta1), delta2(input.delta2),
	delta_ms(input.delta_ms), phi_s(input.phi_s), _cosphis(input._cosphis), _sinphis(input._sinphis), _mistag(input._mistag),
	_mistagP1(input._mistagP1), _mistagP0(input._mistagP0), _mistagSetPoint(input._mistagSetPoint), resolution(input.resolution),
	resolutionScale(input.resolutionScale), resolution1(input.resolution1), resolution2(input.resolution2), resolution3(input.resolution3),
	resolution2Fraction(input.resolution2Fraction), resolution3Fraction(input.resolution3Fraction), timeOffset(input.timeOffset),
	angAccI1(input.angAccI1), angAccI2(input.angAccI2), angAccI3(input.angAccI3), angAccI4(input.angAccI4), angAccI5(input.angAccI5),
	angAccI6(input.angAccI6), angAccI7(input.angAccI7), angAccI8(input.angAccI8), angAccI9(input.angAccI9), angAccI10(input.angAccI10),
	tlo(input.tlo), thi(input.thi), expL_stored(input.expL_stored), expH_stored(input.expH_stored), expSin_stored(input.expSin_stored),
	expCos_stored(input.expCos_stored), intExpL_stored(input.intExpL_stored), intExpH_stored(input.intExpH_stored),
	intExpSin_stored(input.intExpSin_stored), intExpCos_stored(input.intExpCos_stored), _useTimeAcceptance(input._useTimeAcceptance),
	_numericIntegralForce(input._numericIntegralForce), _numericIntegralTimeOnly(input._numericIntegralTimeOnly), 
	_useCosAndSin(input._useCosAndSin), allowNegativeAsSq(input.allowNegativeAsSq), timeAcc(input.timeAcc),
	normalisationCacheValid(input.normalisationCacheValid)
{
	timeAcc = new SlicedAcceptance( *(input.timeAcc) );
}

//........................................................
//Destructor
Bs2Jpsifzero_SignalAlt_BaseClass_dev::~Bs2Jpsifzero_SignalAlt_BaseClass_dev() {}


//--------------------------------------------------------------------------
// Time primitives. These now interface to an external helper library

//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2Jpsifzero_SignalAlt_BaseClass_dev::preCalculateTimeFactors( ) const
{
	expL_stored = Mathematics::Exp( t, gamma_l(), resolution ) ;
	expH_stored = Mathematics::Exp( t, gamma_h(), resolution ) ;
	expSin_stored = Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
	expCos_stored = Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
	return ;
}


//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2Jpsifzero_SignalAlt_BaseClass_dev::preCalculateTimeIntegrals( ) const
{
	intExpL_stored = Mathematics::ExpInt( tlo, thi, gamma_l(), resolution )  ;
	intExpH_stored = Mathematics::ExpInt( tlo, thi, gamma_h(), resolution )  ;
	intExpSin_stored = Mathematics::ExpSinInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
	intExpCos_stored = Mathematics::ExpCosInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
	return ;
}


//-------------------------------------------------------------------------
// Differential crtoss sections and normalisations
//

//...................................
// Main Diff cross section

double Bs2Jpsifzero_SignalAlt_BaseClass_dev::diffXsec(  )  const
{   
	preCalculateTimeFactors() ;
	
	double xsec = 
	
	//0.5 * A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) +
	//0.5 * AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) +
	0.5 * timeFactorATAT(  )  ;
	

	//If this is the very first time here, and it is before a normalisation call has been made, then the time acceptance 
	// object cannot exist.  Only option is to ignore it until it is created - but this should happen within one iteration
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue(t);
				
	if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2Jpsifzero_SignalAlt_BaseClass_dev_v1::diffXsec( ) : return value < 0 = ", xsec ) ;
	
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
// Integral over angles only for a fixed time.

double Bs2Jpsifzero_SignalAlt_BaseClass_dev::diffXsecTimeOnly(  ) const
{          
	preCalculateTimeFactors() ;
	
	double xsec = 
	
	//0.5 * A0()*A0() * timeFactorA0A0(  ) * angAccI1 +
	//0.5 * AP()*AP() * timeFactorAPAP(  ) * angAccI2 +
	0.5 *  timeFactorATAT(  )  ;
	
	
	
	//If this is the very first time here, and it is before a normalisation call has been made, then the time acceptance 
	// object cannot exist.  Only option is to ignore it until it is created - but this should happen within one iteration
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue(t);
	
	if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2Jpsifzero_SignalAlt_BaseClass_dev_v1::diffXsecTimeOnly( ) : return value < 0 = ", xsec ) ;
	
	return xsec ;
}




//...................................
// Integral over all variables: t + angles

double Bs2Jpsifzero_SignalAlt_BaseClass_dev::diffXsecNorm1(  ) const
{ 
	preCalculateTimeIntegrals() ;
	
	double norm =
	
	//0.5 * A0()*A0() * timeFactorA0A0Int(  ) * angAccI1   +  
	//0.5 * AP()*AP() * timeFactorAPAPInt(  ) * angAccI2   +  
	0.5  * timeFactorATATInt(  ) ; 
	
	
	if( DEBUGFLAG && (norm < 0) ) this->DebugPrintNorm( " Bs2Jpsifzero_SignalAlt_BaseClass_dev_v1::diffXsecNorm1( ) : return value < 0 = ", norm ) ;
	 
	return norm ;
}



//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,

double Bs2Jpsifzero_SignalAlt_BaseClass_dev::diffXsecCompositeNorm1(  )  
{   
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;
	double returnValue = 0;
	
	if( useTimeAcceptance() ) {		 
		//This loops over each time slice, does the normalisation between the limits, and accumulates
		for( unsigned int islice = 0; islice < timeAcc->numberOfSlices(); ++islice )
		{
			tlo = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow() ;
			thi = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh() ;
			if( thi > tlo ) returnValue+= this->diffXsecNorm1(  ) * timeAcc->getSlice(islice)->height() ;
		}
	}	
	else {
		returnValue = this->diffXsecNorm1() ;
	}

	tlo = tlo_boundary;
	thi = thi_boundary ;
	return returnValue ;
}


//...............................................................................
// Debug printout
void Bs2Jpsifzero_SignalAlt_BaseClass_dev::DebugPrintXsec( string message, double value )  const
{   
    cout << "*************DEBUG OUTPUT FROM Bs2Jpsifzero_SignalAlt_BaseClass_dev::DebugPrintXsec ***************************" << endl ;
	cout << message << value << endl <<endl ;
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
	


void Bs2Jpsifzero_SignalAlt_BaseClass_dev::DebugPrintNorm( string message, double value )  const
{   
    cout << "*************DEBUG OUTPUT FROM Bs2Jpsifzero_SignalAlt_BaseClass_dev::DebugPrintNorm ***************************" << endl ;
	cout << message << value << endl <<endl ;

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
	
