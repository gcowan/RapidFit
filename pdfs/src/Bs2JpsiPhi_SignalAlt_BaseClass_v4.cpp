// $Id: Bs2JpsiPhi_SignalAlt_BaseClass_v4.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_BaseClass_v4 Bs2JpsiPhi_SignalAlt_BaseClass_v4.cpp
 *
 *  Base class for Bs2JpsiPhi_SignalAlt..... PDFs
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-12
 */

//	ROOT Headers
#include "TMath.h"
#include "RooMath.h"
//	RapidFit Headers
#include "Bs2JpsiPhi_SignalAlt_BaseClass_v4.h"
#include "SlicedAcceptance.h"
//	System Headers
#include <iostream>
#include <math.h>

//......................................
//Constructor(s)

//.....................................
// New Constructor which takes configuration object
Bs2JpsiPhi_SignalAlt_BaseClass_v4::Bs2JpsiPhi_SignalAlt_BaseClass_v4(PDFConfigurator configurator ) : 
// Physics parameters
	  gammaName     		( configurator.getName("gamma") )
	, deltaGammaName		( configurator.getName("deltaGamma") )
	, deltaMName			( configurator.getName("deltaM") )
	, Phi_sName				( configurator.getName("Phi_s") )
	, Azero_sqName			( configurator.getName("Azero_sq") )
	, Apara_sqName			( configurator.getName("Apara_sq") )
	, Aperp_sqName			( configurator.getName("Aperp_sq") )
	, As_sqName				( configurator.getName("As_sq") )
	, delta_zeroName		( configurator.getName("delta_zero") )
	, delta_paraName		( configurator.getName("delta_para") )
	, delta_perpName		( configurator.getName("delta_perp") )
	, delta_sName			( configurator.getName("delta_s") )
	// PELC NEW additions for v2
	, cosphisName			( configurator.getName("cosphis") )
	, sinphisName			( configurator.getName("sinphis") )
	// Detector parameters
	, mistagName			( configurator.getName("mistag") )
	, mistagP1Name			( configurator.getName("mistagP1") )
	, mistagP0Name			( configurator.getName("mistagP0") )
	, mistagSetPointName	( configurator.getName("mistagSetPoint") )
	// Detector parameters
	, resScaleName			( configurator.getName("timeResolutionScale") )
	, res1Name				( configurator.getName("timeResolution1") )
	, res2Name				( configurator.getName("timeResolution2") )
	, res3Name				( configurator.getName("timeResolution3") )
	, res2FractionName		( configurator.getName("timeResolution2Fraction") )
	, res3FractionName		( configurator.getName("timeResolution3Fraction") )
	, timeOffsetName		( configurator.getName("timeOffset") )
	// Angular acceptance factors
	, angAccI1Name			( configurator.getName("angAccI1") )
	, angAccI2Name			( configurator.getName("angAccI2") )
	, angAccI3Name			( configurator.getName("angAccI3") )
	, angAccI4Name			( configurator.getName("angAccI4") )
	, angAccI5Name			( configurator.getName("angAccI5") )
	, angAccI6Name			( configurator.getName("angAccI6") )
	, angAccI7Name			( configurator.getName("angAccI7") )
	, angAccI8Name			( configurator.getName("angAccI8") )
	, angAccI9Name			( configurator.getName("angAccI9") )
	, angAccI10Name			( configurator.getName("angAccI10") )
	// Observables
	, timeName				( configurator.getName("time") )
	, cosThetaName			( configurator.getName("cosTheta") )
	, phiName				( configurator.getName("phi") )
	, cosPsiName			( configurator.getName("cosPsi") )
	, tagName				( configurator.getName("tag") )
	//X, timeAcceptanceCategoryName	( configurator.getName("timeAcceptanceCategory") )
	, timeConstraintName	( configurator.getName("time") )
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
	intExpL_stored(), intExpH_stored(), intExpSin_stored(), intExpCos_stored(), timeAcc(NULL), normalisationCacheValid(false)
{

	//...........................................
	// Configure to use time acceptance machinery 
	_useTimeAcceptance = configurator.isTrue( "UseTimeAcceptance" ) ;
	
	if( useTimeAcceptance() ) {
		if( configurator.hasConfigurationValue( "TimeAcceptanceType", "Upper" ) ) {
			timeAcc = new SlicedAcceptance( 0., 14.0, 0.0157 ) ;
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4:: Constructing timeAcc: Upper time acceptance beta=0.0157 [0 < t < 14] " << endl ;
		}
		else if( configurator.hasConfigurationValue( "TimeAcceptanceType", "Lower2010" ) ) {
			timeAcc = new SlicedAcceptance( "Lower2010" ) ;
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4:: Constructing timeAcc: Lower time acceptance 2010  [0 < t < 14] " << endl ;
		}
		else if( configurator.getConfigurationValue( "TimeAcceptanceFile" ) != "" ) {
			timeAcc = new SlicedAcceptance( "File" , configurator.getConfigurationValue( "TimeAcceptanceFile" ) ) ;
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4:: Constructing timeAcc: using file: " << configurator.getConfigurationValue( "TimeAcceptanceFile" ) << endl ;
		}
	}
	else {
			timeAcc = new SlicedAcceptance( 0., 14. ) ;
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl ;
	}
	
	
	//Make empty Cache for the time integrals. This has to be done now after the SlicedAcceptance is created
	vector<double> empty ;
	for( int islice = 0; islice < timeAcc->numberOfSlices(); ++islice ) empty.push_back(0.0) ;
	for( int ires=0; ires < 4 ; ++ires ) {			
		storeExpL.push_back( empty ) ;
		storeExpH.push_back( empty ) ;
		storeExpSin.push_back( empty ) ;
		storeExpCos.push_back( empty ) ;
	}
	
	
	
	
	//...........................................
	// Configure numerical integration options 
	_numericIntegralForce    = configurator.isTrue( "NumericIntegralForce") ;
	_numericIntegralTimeOnly = configurator.isTrue( "NumericIntegralTimeOnly" ) ;
	
	//...........................................
	// Configure fit parameter options 
	_useCosAndSin = configurator.isTrue( "UseCosAndSin" ) ;
	allowNegativeAsSq = configurator.isTrue( "AllowNegativeAsSq" ) ;
	

	//PELC  - debug to plot the distribution of PDF values for each event 
	//histOfPdfValues = new TH1D( "HistOfPdfValue" ,  "HistOfPdfValue" , 110, -0.00001, 0.00001 ) ;
	//c0  = new TCanvas;
	//histCounter = 0;
	//~PELC
}


//........................................................
//Destructor
Bs2JpsiPhi_SignalAlt_BaseClass_v4::~Bs2JpsiPhi_SignalAlt_BaseClass_v4() {}


//........................................................
//Set the physics parameters into member variables

bool Bs2JpsiPhi_SignalAlt_BaseClass_v4::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;  //This is only used for the untagged events.
	timeIntegralCacheValid = false;  
	
	bool result = allParameters.SetPhysicsParameters(NewParameterSet);
	
	// Physics parameters. 
	_gamma  = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	dgam      = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	
	Azero_sq = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
	if( (Azero_sq < 0.) || (Azero_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MO_v4::SetPhysicsParameters: Azero_sq <0 or >1 but left as is" <<  endl ;	}	
	Aperp_sq = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
	if( (Aperp_sq < 0.) || (Aperp_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MO_v4::SetPhysicsParameters: Aperp_sq <0 or >1 but left as is" <<  endl ;	}	
	As_sq = allParameters.GetPhysicsParameter( As_sqName )->GetValue();
	
	if( allowNegativeAsSq ) {
		if( (As_sq > 1.) ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_dev::SetPhysicsParameters: As_sq >1 but left as is" <<  endl ;	}
		Apara_sq = (1. - Azero_sq - Aperp_sq ) ;
	}
	else {
		if( (As_sq < 0.) || (As_sq > 1.) ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_dev::SetPhysicsParameters: As_sq <0 or >1 but left as is" <<  endl ;	}	
		Apara_sq = (1. - Azero_sq - Aperp_sq  - As_sq) ;
	}
	
	if( Apara_sq < 0. ) {
		cout << "Warning in Bs2JpsiPhi_SignalAlt_MO_v4::SetPhysicsParameters: derived parameter Apara_sq <0  and so set to zero" <<  endl ;
		Apara_sq = 0. ;
	}	
	
	delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	delta_s	   = allParameters.GetPhysicsParameter( delta_sName )->GetValue();
	delta1 = delta_perp -  delta_para ;
	delta2 = delta_perp -  delta_zero ;
	
	delta_ms		= allParameters.GetPhysicsParameter( deltaMName )->GetValue();	
	
	if(_useCosAndSin){
		_cosphis = allParameters.GetPhysicsParameter( cosphisName )->GetValue();
		_sinphis = allParameters.GetPhysicsParameter( sinphisName )->GetValue();
	}
	else{
		phi_s     = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
		_cosphis = cos(phi_s) ;
		_sinphis = sin(phi_s) ;
	}
	
	// Mistag parameters
	_mistagP1		= allParameters.GetPhysicsParameter( mistagP1Name )->GetValue();
	_mistagP0		= allParameters.GetPhysicsParameter( mistagP0Name )->GetValue();
	_mistagSetPoint = allParameters.GetPhysicsParameter( mistagSetPointName )->GetValue();
	
	// Detector parameters
	resolutionScale		= allParameters.GetPhysicsParameter( resScaleName )->GetValue();
	resolution1         = allParameters.GetPhysicsParameter( res1Name )->GetValue();
	resolution2         = allParameters.GetPhysicsParameter( res2Name )->GetValue();
	resolution3         = allParameters.GetPhysicsParameter( res3Name )->GetValue();
	resolution2Fraction = allParameters.GetPhysicsParameter( res2FractionName )->GetValue();
	resolution3Fraction = allParameters.GetPhysicsParameter( res3FractionName )->GetValue();
	timeOffset          = allParameters.GetPhysicsParameter( timeOffsetName )->GetValue();
	
	// Angular acceptance factors
	angAccI1 = allParameters.GetPhysicsParameter( angAccI1Name )->GetValue();
	angAccI2 = allParameters.GetPhysicsParameter( angAccI2Name )->GetValue();
	angAccI3 = allParameters.GetPhysicsParameter( angAccI3Name )->GetValue();
	angAccI4 = allParameters.GetPhysicsParameter( angAccI4Name )->GetValue();
	angAccI5 = allParameters.GetPhysicsParameter( angAccI5Name )->GetValue();
	angAccI6 = allParameters.GetPhysicsParameter( angAccI6Name )->GetValue();
	angAccI7 = allParameters.GetPhysicsParameter( angAccI7Name )->GetValue();
	angAccI8 = allParameters.GetPhysicsParameter( angAccI8Name )->GetValue();
	angAccI9 = allParameters.GetPhysicsParameter( angAccI9Name )->GetValue();
	angAccI10 = allParameters.GetPhysicsParameter( angAccI10Name )->GetValue();
	
	return result;
}






//--------------------------------------------------------------------------
// Time primitives. These now interface to an external helper library

//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_SignalAlt_BaseClass_v4::preCalculateTimeFactors( ) const
{
	expL_stored = Mathematics::Exp( t, gamma_l(), resolution ) ;
	expH_stored = Mathematics::Exp( t, gamma_h(), resolution ) ;
	expSin_stored = Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
	expCos_stored = Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
	return ;
}


//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_SignalAlt_BaseClass_v4::preCalculateTimeIntegrals( ) const
{
	intExpL_stored = Mathematics::ExpInt( tlo, thi, gamma_l(), resolution )  ;
	intExpH_stored = Mathematics::ExpInt( tlo, thi, gamma_h(), resolution )  ;
	intExpSin_stored = Mathematics::ExpSinInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
	intExpCos_stored = Mathematics::ExpCosInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
	return ;
}


//-------------------------------------------------------------------------
// Differential cross sections and normalisations
//

//...................................
// Main Diff cross section

double Bs2JpsiPhi_SignalAlt_BaseClass_v4::diffXsec(  )  const
{   
	preCalculateTimeFactors() ;
	
	double xsec = 
	
	/*
	A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) +
	AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) +
	AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) +
	
	AP()*AT() * timeFactorImAPAT(  ) * angleFactorImAPAT( ) +
	A0()*AP() * timeFactorReA0AP(  ) * angleFactorReA0AP( ) +
	A0()*AT() * timeFactorImA0AT(  ) * angleFactorImA0AT( ) +
	
	AS()*AS() * timeFactorASAS(  ) * angleFactorASAS( ) +
	
	AS()*AP() * timeFactorReASAP(  ) * angleFactorReASAP( ) +
	AS()*AT() * timeFactorImASAT(  ) * angleFactorImASAT( ) +
	AS()*A0() * timeFactorReASA0(  ) * angleFactorReASA0( ) ;	
	*/
	
	CachedA1 * timeFactorA0A0(  )  +
	CachedA2 * timeFactorAPAP(  )  +
	CachedA3 * timeFactorATAT(  )  +
	
	CachedA4 * timeFactorImAPAT(  )  +
	CachedA5 * timeFactorReA0AP(  )  +
	CachedA6 * timeFactorImA0AT(  )  +
	
	CachedA7 * timeFactorASAS(  )  +
	
	CachedA8 * timeFactorReASAP(  )  +
	CachedA9 * timeFactorImASAT(  )  +
	CachedA10 * timeFactorReASA0(  )  ;	
	
	
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue(t);
				
	if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2JpsiPhi_SignalAlt_BaseClass_v4_v1::diffXsec( ) : return value < 0 = ", xsec ) ;
	
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

double Bs2JpsiPhi_SignalAlt_BaseClass_v4::diffXsecTimeOnly(  ) const
{          
	preCalculateTimeFactors() ;
	
	double xsec = 
	
	A0()*A0() * timeFactorA0A0(  ) * angAccI1 +
	AP()*AP() * timeFactorAPAP(  ) * angAccI2 +
	AT()*AT() * timeFactorATAT(  ) * angAccI3 +
	
	AP()*AT() * timeFactorImAPAT(  ) * angAccI4 +
	A0()*AP() * timeFactorReA0AP(  ) * angAccI5 +
	A0()*AT() * timeFactorImA0AT(  ) * angAccI6 +
	
	AS()*AS() * timeFactorASAS(  ) * angAccI7 +
	
	AS()*AP() * timeFactorReASAP(  ) * angAccI8 +
	AS()*AT() * timeFactorImASAT(  ) * angAccI9 +
	AS()*A0() * timeFactorReASA0(  ) * angAccI10 ;
	
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue(t);
	
	if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2JpsiPhi_SignalAlt_BaseClass_v4_v1::diffXsecTimeOnly( ) : return value < 0 = ", xsec ) ;
	
	return xsec ;
}




//...................................
// Integral over all variables: t + angles

double Bs2JpsiPhi_SignalAlt_BaseClass_v4::diffXsecNorm1(  ) const
{ 
	//preCalculateTimeIntegrals() ;  Replaced by new Caching mechanism  
	
	double norm =
	
	A0()*A0() * timeFactorA0A0Int(  ) * angAccI1   +  
	AP()*AP() * timeFactorAPAPInt(  ) * angAccI2   +  
	AT()*AT() * timeFactorATATInt(  ) * angAccI3   +  
	
	AP()*AT() * timeFactorImAPATInt(  ) * angAccI4 +  
	A0()*AP() * timeFactorReA0APInt(  ) * angAccI5 +  
	A0()*AT() * timeFactorImA0ATInt(  ) * angAccI6 +  
	
	AS()*AS() * timeFactorASASInt(  ) * angAccI7   +  
	
	AS()*AP() * timeFactorReASAPInt(  ) * angAccI8 +  
	AS()*AT() * timeFactorImASATInt(  ) * angAccI9 +  
	AS()*A0() * timeFactorReASA0Int(  ) * angAccI10 ;  
	
	if( DEBUGFLAG && (norm < 0) ) this->DebugPrintNorm( " Bs2JpsiPhi_SignalAlt_BaseClass_v4_v1::diffXsecNorm1( ) : return value < 0 = ", norm ) ;
	 
	return norm ;
}



//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,

double Bs2JpsiPhi_SignalAlt_BaseClass_v4::diffXsecCompositeNorm1( int resolutionIndex )  
{   
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;
	double returnValue = 0;
	
	
	if( true /*useTimeAcceptance()*/ ) {		    // Set to true because seleting false makes a single slice for 0 --> 14. 
		//This loops over each time slice, does the normalisation between the limits, and accumulates
		for( int islice = 0; islice < timeAcc->numberOfSlices(); ++islice )
		{
			//Set the time integrals
			this->deCacheTimeIntegrals( resolutionIndex, islice ) ;
			
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



//.......................................................
// New speed up method to Cache time integrals
void Bs2JpsiPhi_SignalAlt_BaseClass_v4::CacheAmplitudesAndAngles() {

	CachedA1 = A0()*A0() * angleFactorA0A0( ) ;
	CachedA2 = AP()*AP() * angleFactorAPAP( ) ;
	CachedA3 = AT()*AT() * angleFactorATAT( ) ;
	
	CachedA4 = AP()*AT() * angleFactorImAPAT( ) ;
	CachedA5 = A0()*AP() * angleFactorReA0AP( ) ;
	CachedA6 = A0()*AT() * angleFactorImA0AT( ) ;
	
	CachedA7 = AS()*AS() * angleFactorASAS( ) ;
	
	CachedA8 = AS()*AP() * angleFactorReASAP( ) ;
	CachedA9 = AS()*AT() * angleFactorImASAT( ) ;
	CachedA10= AS()*A0() * angleFactorReASA0( ) ;	
	
}

//.......................................................
// New speed up method to Cache time integrals
void Bs2JpsiPhi_SignalAlt_BaseClass_v4::CacheTimeIntegrals() {
	
	// This need to know  (and be modified for)
	//  --> Number of resolutions
	//  --> Number fo slices
	//  --->  tlo and thi for each slice
	
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;
	
	for( int ires=0; ires < 4 ; ++ires ) {
		
		if( ires==0 ) resolution = 0.0 ;
		if( ires==1 ) resolution = resolution1 * resolutionScale ;
		if( ires==2 ) resolution = resolution2 * resolutionScale ;
		if( ires==3 ) resolution = resolution3 * resolutionScale ;
		
		for( int islice = 0; islice < timeAcc->numberOfSlices(); ++islice ) {
			tlo = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow() ;
			thi = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh() ;
			if( thi > tlo ) {
				this->preCalculateTimeIntegrals() ;
				//cout << " >>>>> caching time integrals / " << intExpL_stored << "  /  "<< intExpH_stored << "  /  "<< intExpSin_stored << "  /  "<< intExpCos_stored << "  /  " << endl ;
				
				storeExpL[ires][islice] = intExpL_stored ;
				storeExpH[ires][islice] = intExpH_stored ;
				storeExpSin[ires][islice] = intExpSin_stored ;
				storeExpCos[ires][islice] = intExpCos_stored ;
			}
			else {
				storeExpL[ires][islice] = 0 ;
				storeExpH[ires][islice] = 0 ;
				storeExpSin[ires][islice] = 0 ;
				storeExpCos[ires][islice] = 0 ;
			}
			
		}
	}
	
	tlo = tlo_boundary;
	thi = thi_boundary ;
	
	
}		


//.......................................................
// New speed up method to Cache time integrals
void Bs2JpsiPhi_SignalAlt_BaseClass_v4::deCacheTimeIntegrals( int ires, int islice ) {
	
	//Time integrals are stored under 
	// ires =  resolution integral
	// islice = acceptance slice
	
	intExpL_stored   = storeExpL[ires][islice]  ;
	intExpH_stored   = storeExpH[ires][islice]  ;
	intExpSin_stored = storeExpSin[ires][islice]  ;
	intExpCos_stored = storeExpCos[ires][islice]  ;
	
	//cout << " <<<<< de-caching time integrals / " << intExpL_stored << "  /  "<< intExpH_stored << "  /  "<< intExpSin_stored << "  /  "<< intExpCos_stored << "  /  " << endl ;
	
}		



//...............................................................................
// Debug printout
void Bs2JpsiPhi_SignalAlt_BaseClass_v4::DebugPrintXsec( string message, double value )  const
{   
    cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_SignalAlt_BaseClass_v4::DebugPrintXsec ***************************" << endl ;
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
	A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) +
	AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) +
	AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) +		
	AP()*AT() * timeFactorImAPAT(  ) * angleFactorImAPAT( ) +
	A0()*AP() * timeFactorReA0AP(  ) * angleFactorReA0AP( ) +
	A0()*AT() * timeFactorImA0AT(  ) * angleFactorImA0AT( ) ;
	
	double SwaveAdditions =
	AS()*AS() * timeFactorASAS(  ) * angleFactorASAS( ) +
	AS()*AP() * timeFactorReASAP(  ) * angleFactorReASAP( ) +
	AS()*AT() * timeFactorImASAT(  ) * angleFactorImASAT( ) +
	AS()*A0() * timeFactorReASA0(  ) * angleFactorReASA0( ) ;
	
	cout << "   Pwave Only : " << PwaveTot << endl ;
	cout << "   Swave add : " <<  SwaveAdditions << endl ;
	
}
	


void Bs2JpsiPhi_SignalAlt_BaseClass_v4::DebugPrintNorm( string message, double value )  const
{   
    cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_SignalAlt_BaseClass_v4::DebugPrintNorm ***************************" << endl ;
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





	
