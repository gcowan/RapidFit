// $Id: Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4 Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4.cpp
 *
 *  Base class for Bs2JpsiPhi_SignalAlt..... PDFs
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-12
 */

//	ROOT Headers
#include "TMath.h"
//	RapidFit Headers
#include "Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4.h"
#include "SlicedAcceptance.h"
//	System Headers
#include <iostream>
#include <math.h>

//......................................
//Constructor(s)

//.....................................
// New Constructor which takes configuration object
Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4( PDFConfigurator* configurator ) : 
// Physics parameters
	  gammaName     		( configurator->getName("gamma") )
	, deltaGammaName		( configurator->getName("deltaGamma") )
	, deltaMName			( configurator->getName("deltaM") )
	, Phi_sName				( configurator->getName("Phi_s") )
	, Aeven_sqName			( configurator->getName("Aeven_sq") )
	, Aodd_sqName			( configurator->getName("Aodd_sq") )
	, As_sqName				( configurator->getName("As_sq") )
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
	, angAccEvenName		( configurator->getName("angAccEven") )
	, angAccOddName			( configurator->getName("angAccOdd") )
	// Observables
	, timeName				( configurator->getName("time") )
	, cosThetaName			( configurator->getName("cosTheta") )
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
	,t(), ctheta_tr(), tag(), /*timeAcceptanceCategory(),*/ _gamma(), dgam(), Aeven_sq(), Aodd_sq(), As_sq(), 
	delta_ms(), phi_s(), _cosphis(), _sinphis(), _mistag(), _mistagP1(), _mistagP0(), _mistagSetPoint(),
	resolutionScale(), resolution1(), resolution2(), resolution3(), resolution2Fraction(), resolution3Fraction(), timeOffset(), 
	angAccEven(), angAccOdd(),  
	tlo(), thi(), expL_stored(), expH_stored(), expSin_stored(), expCos_stored(),
	intExpL_stored(), intExpH_stored(), intExpSin_stored(), intExpCos_stored(), timeAcc(NULL), normalisationCacheValid(false)
{

	//...........................................
	// Configure to use time acceptance machinery 
	_useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;
	
	if( useTimeAcceptance() ) {
		if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) ) {
			timeAcc = new SlicedAcceptance( 0., 14.0, 0.0157 ) ;
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4:: Constructing timeAcc: Upper time acceptance beta=0.0157 [0 < t < 14] " << endl ;
		}
		else if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Lower2010" ) ) {
			timeAcc = new SlicedAcceptance( "Lower2010" ) ;
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4:: Constructing timeAcc: Lower time acceptance 2010  [0 < t < 14] " << endl ;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" ) {
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ) ) ;
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl ;
		}
	}
	else {
			timeAcc = new SlicedAcceptance( 0., 14. ) ;
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl ;
	}
	
	
	//Make empty Cache for the time integrals. This has to be done now after the SlicedAcceptance is created
	vector<double> empty ;
	for( unsigned int islice = 0; islice < timeAcc->numberOfSlices(); ++islice ) empty.push_back(0.0) ;
	for( unsigned int ires=0; ires < 4 ; ++ires ) {			
		storeExpL.push_back( empty ) ;
		storeExpH.push_back( empty ) ;
		storeExpSin.push_back( empty ) ;
		storeExpCos.push_back( empty ) ;
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


//........................................................
//Destructor
Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::~Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4() {}


//........................................................
//Set the physics parameters into member variables

bool Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;  //This is only used for the untagged events.
	timeIntegralCacheValid = false;  
	
	bool result = allParameters.SetPhysicsParameters(NewParameterSet);
	
	// Physics parameters. 
	_gamma  = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	dgam      = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	
	//Aeven_sq = allParameters.GetPhysicsParameter( Aeven_sqName )->GetValue();
	//if( (Aeven_sq < 0.) || (Aeven_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MO_1angle_v4::SetPhysicsParameters: Aeven_sq <0 or >1 but left as is" <<  endl ;	}	
	Aodd_sq = allParameters.GetPhysicsParameter( Aodd_sqName )->GetValue();
	if( (Aodd_sq < 0.) || (Aodd_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MO_1angle_v4::SetPhysicsParameters: Aodd_sq <0 or >1 but left as is" <<  endl ;	}	
	As_sq = allParameters.GetPhysicsParameter( As_sqName )->GetValue();	
	if( allowNegativeAsSq ) {
		if( (As_sq > 1.) ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MO_1angle_v4::SetPhysicsParameters: As_sq >1 but left as is" <<  endl ;	}
	}
	else {
		if( (As_sq < 0.) || (As_sq > 1.) ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MO_1angle_v4::SetPhysicsParameters: As_sq <0 or >1 but left as is" <<  endl ;	}	
	}
	Aeven_sq = 1.0 - Aodd_sq - As_sq ;
	
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
	angAccEven = allParameters.GetPhysicsParameter( angAccEvenName )->GetValue();
	angAccOdd = allParameters.GetPhysicsParameter( angAccOddName )->GetValue();
	
	return result;
}






//--------------------------------------------------------------------------
// Time primitives. These now interface to an external helper library

//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::preCalculateTimeFactors( ) const
{
	expL_stored = Mathematics::Exp( t, gamma_l(), resolution ) ;
	expH_stored = Mathematics::Exp( t, gamma_h(), resolution ) ;
	expSin_stored = Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
	expCos_stored = Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
	return ;
}


//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::preCalculateTimeIntegrals( ) const
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

double Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::diffXsec(  )  const
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
	
	CachedAEven * timeFactorEven(  )  +
	CachedAOdd  * timeFactorOdd(  )  +
	CachedAs    * timeFactorOdd(  )  ;
	
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue(t);
				
	if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::diffXsec( ) : return value < 0 = ", xsec ) ;
	
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

double Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::diffXsecTimeOnly(  ) const
{          
	preCalculateTimeFactors() ;
	
	double xsec = 
	
	Aeven_sq * timeFactorEven(  ) * angAccEven +
	Aodd_sq  * timeFactorOdd(  ) * angAccOdd +
	As_sq	 * timeFactorOdd(  ) * angAccEven ;    // yes this really is an even angle term
		
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue(t);
	
	if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4_v4::diffXsecTimeOnly( ) : return value < 0 = ", xsec ) ;
	
	return xsec ;
}




//...................................
// Integral over all variables: t + angles

double Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::diffXsecNorm1(  ) const
{ 
	//preCalculateTimeIntegrals() ;  Replaced by new Caching mechanism  
	
	double norm =
	
	Aeven_sq * timeFactorEvenInt(  ) * angAccEven   +  
	Aodd_sq  * timeFactorOddInt(  ) * angAccOdd   +  
	As_sq	* timeFactorOddInt(  ) * angAccEven   ; 
		
	if( DEBUGFLAG && (norm < 0) ) this->DebugPrintNorm( " Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4_v1::diffXsecNorm1( ) : return value < 0 = ", norm ) ;
	 
	return norm ;
}



//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,

double Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::diffXsecCompositeNorm1( unsigned int resolutionIndex )  
{   
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;
	double returnValue = 0;
	
	
	if( true /*useTimeAcceptance()*/ ) {		    // Set to true because seleting false makes a single slice for 0 --> 14. 
		//This loops over each time slice, does the normalisation between the limits, and accumulates
		for( unsigned int islice = 0; islice < timeAcc->numberOfSlices(); ++islice )
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
void Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::CacheAmplitudesAndAngles() {

	CachedAEven = Aeven_sq * angleFactorEven( ) ;
	CachedAOdd  = Aodd_sq  * angleFactorOdd( ) ;
	CachedAs    = As_sq    * angleFactorEven( ) ;
	
}

//.......................................................
// New speed up method to Cache time integrals
void Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::CacheTimeIntegrals() {
	
	// This need to know  (and be modified for)
	//  --> Number of resolutions
	//  --> Number fo slices
	//  --->  tlo and thi for each slice
	
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;
	
	for( unsigned int ires=0; ires < 4 ; ++ires ) {
		
		if( ires==0 ) resolution = 0.0 ;
		if( ires==1 ) resolution = resolution1 * resolutionScale ;
		if( ires==2 ) resolution = resolution2 * resolutionScale ;
		if( ires==3 ) resolution = resolution3 * resolutionScale ;
		
		for( unsigned int islice = 0; islice < timeAcc->numberOfSlices(); ++islice ) {
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
void Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::deCacheTimeIntegrals( unsigned int ires, unsigned int islice ) {
	
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
void Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::DebugPrintXsec( string message, double value )  const
{   
    cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::DebugPrintXsec ***************************" << endl ;
	cout << message << value << endl <<endl ;
/*
	cout << "   A0()*A0() term: " <<  A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) << endl ;
	cout << "   AP()*AP() term: " <<AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) << endl ;
	cout << "   AT()*AT() term: " <<AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) << endl << endl ;
	
	
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
*/	
}
	


void Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::DebugPrintNorm( string message, double value )  const
{   
    cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::DebugPrintNorm ***************************" << endl ;
	cout << message << value << endl <<endl ;
/*
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
 */
}





	
