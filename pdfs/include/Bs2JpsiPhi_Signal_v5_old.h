// $Id: Bs2JpsiPhi_Signal_v5_old.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2JpsiPhi_Signal_v5_old 
 *
 *  Bs2JpsiPhi_SignalAlt series with mistag as observable
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2011-11-05
 */

#ifndef Bs2JpsiPhi_Signal_v5_old_H
#define Bs2JpsiPhi_Signal_v5_old_H

#include <exception>

#include "BasePDF.h"
#include "PDFConfigurator.h"
#include "SlicedAcceptance.h"
#include "AngularAcceptance.h"
#include "Mathematics.h"
#include "RooComplex.h"
#include <iostream>
#include <cstdlib>
#include <float.h>
#include <vector>

#include "TH1.h"
#include "TCanvas.h"

#define DOUBLE_TOLERANCE 1E-6
#define DEBUGFLAG true


//============================================================================================
class Bs2JpsiPhi_Signal_v5_old : public BasePDF 
{
public:
	Bs2JpsiPhi_Signal_v5_old( PDFConfigurator* ); 
	//Bs2JpsiPhi_Signal_v5_old( const Bs2JpsiPhi_Signal_v5& );
	~Bs2JpsiPhi_Signal_v5_old();

	//Mandatory RapidFit Methods
	virtual double EvaluateForNumericIntegral(DataPoint*) ;
	virtual double Evaluate(DataPoint*);
	virtual double EvaluateTimeOnly(DataPoint*) ;
	virtual bool SetPhysicsParameters(ParameterSet*);
	virtual vector<string> GetDoNotIntegrateList();

	vector<string> PDFComponents();

	double EvaluateComponent( DataPoint* input, ComponentRef* );
	
private:
	Bs2JpsiPhi_Signal_v5_old& operator=( const Bs2JpsiPhi_Signal_v5_old& );
	void MakePrototypes();
	double normalisationCacheUntagged ;
	void prepareCDS() ;
	
protected:
	//Calculate the PDF normalisation
	virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	int componentIndex;	
	
	//PELC For debugging purposes
	//TH1D * histOfPdfValues ;
	//TCanvas * c0 ; 
	//mutable int histCounter ;
	//~PELC

	// Parameters
	ObservableRef gammaName;		// gamma
	ObservableRef deltaGammaName;	// delta gamma
	ObservableRef deltaMName;		// delta mass
	ObservableRef Azero_sqName;	// amplitude
	ObservableRef Apara_sqName;	// amplitude
	ObservableRef Aperp_sqName;	// amplitude
	ObservableRef As_sqName;		// amplitude
	ObservableRef CspName;
	ObservableRef delta_zeroName;	// strong phase, set to 0
	ObservableRef delta_paraName;	// strong phase
	ObservableRef delta_perpName;	// strong phase
	ObservableRef delta_sName;		// strong phase for S-wave
	ObservableRef cosdparName;		//PELC-COSDPAR Special for fitting cosdpar separately
	
	
	ObservableRef Phi_sName;		// what we want to measure!
	ObservableRef cosphisName;		// fitting cosphis and sinphis independently
	ObservableRef sinphisName;		// fitting cosphis and sinphis independently	
	ObservableRef lambdaName;		// magnitude of lambda	
	
	ObservableRef mistagName;		// mistag fraction  - may be used as observable also
	ObservableRef mistagP1Name;		// mistag calib
	ObservableRef mistagP0Name;		// mistag calib
	ObservableRef mistagSetPointName;// mistag calib
	ObservableRef mistagDeltaP1Name;		// mistag calib
	ObservableRef mistagDeltaP0Name;		// mistag calib
	ObservableRef mistagDeltaSetPointName;// mistag calib
	
	ObservableRef eventResolutionName;			// Scale to multiply all Gaussians with 
	ObservableRef resScaleName;			// Scale to multiply all Gaussians with 
	ObservableRef res1Name;				// time resolution narrow
	ObservableRef res2Name;				// time resolution wide
	ObservableRef res3Name;				// time resolution tail
	ObservableRef res2FractionName;		// fraction of wide
	ObservableRef res3FractionName;		// fraction of tail
	ObservableRef timeOffsetName;		// time offset
	
	ObservableRef angAccI1Name ;  
	ObservableRef angAccI2Name ;
	ObservableRef angAccI3Name ;
	ObservableRef angAccI4Name ;
	ObservableRef angAccI5Name ;
	ObservableRef angAccI6Name ;
	ObservableRef angAccI7Name ;
	ObservableRef angAccI8Name ;
	ObservableRef angAccI9Name ;
	ObservableRef angAccI10Name ;
	
	// Observables 
	ObservableRef timeName;		// proper time
	ObservableRef tagName;			// B tag
	ObservableRef cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
	ObservableRef cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
	ObservableRef phiName;			// azimuthal angle of the mu+ in Jpsi frame
	ObservableRef cthetakName;	
	ObservableRef cthetalName;		
	ObservableRef phihName;			
	
	// Measured Event Observables
	double t ;
	int tag ;
	// Transversity angles
	double ctheta_tr ;
	double phi_tr ;
	double ctheta_1 ;
	// Helicity angles
	double ctheta_k;
	double ctheta_l;
	double phi_h;
	bool _useHelicityBasis ;
	
	
	// Physics Fit Parameters 
	double _gamma ;
	double dgam ;

	double Aperp_sq ;
	double Apara_sq ;
	double Azero_sq ;
	double As_sq ;
	double CachedA1 ;
	double CachedA2 ;
	double CachedA3 ;	
	double CachedA4 ;
	double CachedA5 ;
	double CachedA6 ;
	double CachedA7 ;
	double CachedA8 ;
	double CachedA9 ;
	double CachedA10 ;
	void CacheAmplitudesAndAngles() ;
	
	double delta_para ;
	double delta_perp ;
	double delta_zero ;
	double delta_s ;
	double delta1 ;
	double delta2 ;
	double cosdpar ; //PELC-COSDPAR Special for fitting cosdpar separately
	
	double delta_ms ;
	double phi_s ;
	double _cosphis ;
	double _sinphis ;
	double lambda ;
	double _CC, _DD, _SS;

	double _mistag ;
	double _mistagP1 ;
	double _mistagP0 ;
	double _mistagSetPoint ;
	double _mistagDeltaP1 ;
	double _mistagDeltaP0 ;
	double _mistagDeltaSetPoint ;
	;
	double resolution ;
	double eventResolution ;
	double resolutionScale ;
	double resolution1 ;
	double resolution2 ;
	double resolution3 ;
	double resolution2Fraction ;
	double resolution3Fraction ;
	double timeOffset ;
	bool _useEventResolution ;
	inline bool useEventResolution() const {return _useEventResolution ; }
	inline bool useTimeAcceptance() const { return _useTimeAcceptance ; }		
	

	double angAccI1 ;
	double angAccI2 ;
	double angAccI3 ;
	double angAccI4 ;
	double angAccI5 ;
	double angAccI6 ;
	double angAccI7 ;
	double angAccI8 ;
	double angAccI9 ;
	double angAccI10 ;
	AngularAcceptance * angAcc ;
	bool _angAccIgnoreNumerator ;
	
	// Othere things calculated later on the fly
	double tlo, thi ;
	
	// stored time primitives
	mutable double expL_stored ;
	mutable double expH_stored ;
	mutable double expSin_stored ;
	mutable double expCos_stored ;
	mutable double intExpL_stored ;
	mutable double intExpH_stored ;
	mutable double intExpSin_stored ;
	mutable double intExpCos_stored ;
	void preCalculateTimeFactors() const ;
	void preCalculateTimeIntegrals() const ;
		
	bool timeIntegralCacheValid ;
	vector< vector<double> > storeExpL;
	vector< vector<double> > storeExpH;
	vector< vector<double> > storeExpSin;
	vector< vector<double> > storeExpCos;
	void CacheTimeIntegrals() ;
	void deCacheTimeIntegrals( unsigned int ires, unsigned int islice ) ;
	
	//Time acceptance 
	SlicedAcceptance * timeAcc ;
	
	//Configurationparameters
	bool _useTimeAcceptance ;	
	bool _numericIntegralForce ;
	bool _numericIntegralTimeOnly ;
	bool _useCosAndSin ;
	bool _useCosDpar ;
	bool _usePunziSigmat ;
	bool _usePunziMistag ;
	bool allowNegativeAsSq ;
	
	//....................................
	//Internal helper functions
	
	inline double AT() const { 
		if( Aperp_sq <= 0. ) return 0. ;
		else return sqrt(Aperp_sq) ; 
	}
	inline double AP() const { 
		if( Apara_sq <= 0. ) return 0. ;
		else return sqrt(Apara_sq) ; 
	}
	inline double A0() const { 
		if( Azero_sq <= 0. ) return 0. ;
		else return sqrt(Azero_sq) ; 
	}
	inline double AS() const { 
		if( As_sq <= 0. ) return 0. ;
		else return sqrt(As_sq) ; 
	}
	double Csp;
	inline double ASint() const {
		if( As_sq <= 0. ) return 0. ;
		else return sqrt(As_sq)*Csp;
	}
	
	// Functions required for transversity*******************
	inline double ctrsq() const { return (ctheta_tr*ctheta_tr) ; }
	inline double strsq() const { return (1.0 - ctrsq()) ; }
	inline double theta_tr() const { return acos(ctheta_tr) ; }	
	inline double ctr() const { return ctheta_tr ; }
	inline double str() const { return sin(theta_tr()) ; }
	inline double s2tr() const { return sin(2.0*theta_tr()) ; }
	
	inline double ct1sq() const { return (ctheta_1*ctheta_1) ; }
	inline double st1sq() const { return (1.0 - ct1sq()) ; }
	inline double theta_1() const { return acos(ctheta_1) ; }	
	inline double ct1() const { return ctheta_1 ; }
	inline double st1() const { return sin(theta_1()) ; }
	inline double s2t1() const { return sin(2.0*theta_1()) ; }
	
	inline double cph() const {  return cos(phi_tr) ; }
	inline double sph() const {  return sin(phi_tr) ; }
	inline double cphsq() const { return (cos(phi_tr)*cos(phi_tr)) ; }
	inline double sphsq() const { return (sin(phi_tr)*sin(phi_tr)) ; }
	inline double s2ph() const { return sin(2.0*phi_tr) ; }	
	
	// Functions required for helicity*******************
	inline double cHtksq() const { return (ctheta_k*ctheta_k) ; }
	inline double sHtksq() const { return (1.0 - cHtksq() ) ; }
	inline double cHtk()	const { return ctheta_k ; }
	inline double sHtk()	const { return sin(acos(ctheta_k)) ; }
	inline double c2Htk() const { return (cos(2.*acos(ctheta_k))) ; }
	inline double s2Htk() const { return (sin(2.*acos(ctheta_k))) ; }

	inline double cHtlsq() const { return (ctheta_l*ctheta_l) ; }
	inline double sHtlsq() const { return (1.0 - cHtlsq()) ; }
	inline double c2Htl() const { return (cos(2.*acos(ctheta_l))) ; }
	inline double s2Htl() const { return (sin(2.*acos(ctheta_l))) ; }

	inline double cHphi() const { return (cos(phi_h)); }
	inline double sHphi() const { return (sin(phi_h)); }
	inline double c2Hphi() const { return (cos(2.*phi_h)); }
	inline double s2Hphi() const { return (sin(2.*phi_h)); }
	
	//....................
	// Safety for gammas
	inline double gamma_l() const { 
		const double gl = gamma() + ( dgam *0.5 ) ;
		if( gl < 0. ) {
			cerr << " In Bs2JpsiPhi_SignalAlt_BaseClass_v4 : gamma_l() < 0 so setting it to 0.0000001 " << endl ;
			return 0.0000001 ;
		}
		else
			return gl ; 
	}
	
	inline double gamma_h() const { 
		const double gh = gamma() - ( dgam *0.5 ) ;
		if( gh < 0. ) {
			cerr << " In Bs2JpsiPhi_SignalAlt_BaseClass_v4 : gamma_h() < 0 so setting it to 0.0000001 " << endl ;
			return 0.0000001 ;
		}
		else
			return gh ;   
	}
	
	inline double gamma() const { return _gamma ; }
	
	//...............
	//tagging
	
	
	inline double q() const { return tag ;}
	
	inline double mistag() const { 
		double returnValue = -1000.;
		
		if( fabs((q()-0.0)) < DOUBLE_TOLERANCE ) {
			returnValue = 0.5 ;
		}			
		else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
			//Normal case
			returnValue =  _mistagP0 + _mistagP1*(_mistag - _mistagSetPoint ) ;
			if( returnValue < 0 )  returnValue = 0 ;
			if( returnValue > 0.5) returnValue = 0.5 ; 
		}			
		else if( _mistag < 0.0 ) {
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistag() : _mistag < 0 so set to 0 " << endl ;
			returnValue = 0 ;
		}
		else if( _mistag > 0.5 ) {
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistag() : _mistag > 0.5 so set to 0.5 "  << endl ;
			returnValue = 0.5 ;
		}
		else {
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistag() : WARNING ******If you got here you dont know what you are doing  "  << endl ;
			exit(1);
		}
		return returnValue ;			
	}

	inline double mistagB() const { 
		double returnValue = -1000.;
		
		if( fabs(q()) < 0.5 ) {
			returnValue = 0.5 ;
		}			
		else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
			//Normal case
			returnValue =  _mistagP0+(_mistagDeltaP0/2.0) + (_mistagP1+(_mistagDeltaP1/2.0))*(_mistag - (_mistagSetPoint+(_mistagDeltaSetPoint/2.0)) ) ;
			//if( true ) returnValue =  _mistagP0 + (_mistagP1)*(_mistag - (_mistagSetPoint) ) ;  // to mock up independent P1/P0 for each tag
			if( returnValue < 0 )  returnValue = 0 ;
			if( returnValue > 0.5) returnValue = 0.5 ; 
		}			
		else if( _mistag < 0.0 ) {
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistagB() : _mistag < 0 so deltaMistag set to 0 also " << endl ;
			returnValue = 0 ;
		}
		else if( _mistag > 0.5 ) {
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistagB() : _mistag > 0.5 so so deltaMistag set to 0.5 also "  << endl ;
			returnValue = 0.5 ;
		}
		else {
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistagB() : WARNING ******If you got here you dont know what you are doing  "  << endl ;
			exit(1);
		}
		return returnValue ;			
	}
	
	inline double mistagBbar() const { 
		double returnValue = -1000.;
		
		if( fabs(q()) < 0.5 ) {
			returnValue = 0.5 ;
		}			
		else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
			//Normal case
			returnValue =  _mistagP0-(_mistagDeltaP0/2.0) + (_mistagP1-(_mistagDeltaP1/2.0))*(_mistag - (_mistagSetPoint-(_mistagDeltaSetPoint/2.0)) ) ;
			//if( true ) returnValue =   _mistagDeltaP0 + (_mistagDeltaP1)*(_mistag - (_mistagDeltaSetPoint) ) ;// to mock up independent P1/P0 for each tag
			if( returnValue < 0 )  returnValue = 0 ;
			if( returnValue > 0.5) returnValue = 0.5 ; 
		}			
		else if( _mistag < 0.0 ) {
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistagBbar() : _mistag < 0 so deltaMistag set to 0 also " << endl ;
			returnValue = 0 ;
		}
		else if( _mistag > 0.5 ) {
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistagBbar() : _mistag > 0.5 so so deltaMistag set to 0.5 also "  << endl ;
			returnValue = 0.5 ;
		}
		else {
			cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistagBbar() : WARNING ******If you got here you dont know what you are doing  "  << endl ;
			exit(1);
		}
		return returnValue ;			
	}
	
	inline double D1() const {  return 1.0 - q()*(mistagB()-mistagBbar()) ; }  
	inline double D2() const {  return q()*( 1.0 - mistagB() -mistagBbar() ) ; }
	//inline double D1() const {  return 1.0 ; }  
	//inline double D2() const {  return  ( q()*( 1.0 - mistagB() -mistagBbar() ) ) / ( 1.0 - q()*(mistagB()-mistagBbar())  )  ; }
	
	
	//.....................
	// C, D, S
	inline double cosphis() const { return _DD ; } //  _cosphis ; }
	inline double sinphis() const { return _SS ; } //  _sinphis ; }
	inline double CC() const { return _CC ; } //  _sinphis ; }
	
		
	//......................................................
	// Time primitives
	
	inline double expL() const { return expL_stored ; }
	inline double intExpL( ) const { return intExpL_stored ; }
	
	inline double expH() const { return expH_stored ; }
	inline double intExpH( ) const { return intExpH_stored ; }
	
	inline double expSin() const  { return expSin_stored ; }
	inline double intExpSin( ) const { return intExpSin_stored ;  }
	
	inline double expCos() const { return expCos_stored ; }
	inline double intExpCos( ) const { return intExpCos_stored ; }
	
	
	
	//---------------------------------------------------------
	//............. Differential cross sections and normalisations
	double diffXsec(  )  const ;   	
	double diffXsecTimeOnly(  ) const ;
	double diffXsecNorm1(  ) const ;
	double diffXsecCompositeNorm1( int resolutionIndex )  ;
	
	bool normalisationCacheValid ;
	double normalisationCacheValue[3] ;
	//double normalisationCacheValueRes2[3] ;
	
	void DebugPrint( string , double ) const ;
	void DebugPrintXsec( string , double ) const ;
	void DebugPrintNorm( string , double ) const ;
	
	
	//------------------------------------------------------------------------------
	// These are the time factors and their analytic integrals for the one angle PDF
	
	//..................................
	inline double timeFactorEven(  )  const
	{
		//if( t < 0.0 ) return 0.0 ;
		const double result = 
		D1() * (
			   ( 1.0 + cosphis() ) * expL( ) 
			 + ( 1.0 - cosphis() ) * expH( ) 
		) +
		D2() * ( 
				( 2.0 * sinphis() ) * expSin( )  
			  + ( 2.0 * CC()      ) * expCos( ) 
		);
		
		//DEBUG
		if( DEBUGFLAG && (result < 0) ) {
			cout << " Bs2JpsiPhi_SignalAlt_BaseClass_v4::timeFactorEven() : result < 0 " << endl ;
			cout << " ->term1 " << ( 1.0 + cosphis() ) * expL( ) << endl ;
			cout << " ->term2 " << ( 1.0 - cosphis() ) * expH( ) << endl ;
			cout << " ->term3 " << q() * ( 2.0 * sinphis()   ) * expSin( ) * (1.0 - 2.0*mistag()) << endl ;
			cout << "   -->sin(phis) "  << sinphis() << endl ;
			cout << "   -->expSin    "  << expSin() << endl ;
			cout << "   -->tagFrac   "  << mistag() << endl ;
			cout << "   -->delta_ms  "  << delta_ms << endl ;
		}
		return result ;
	}
	
	inline double timeFactorEvenInt(  )  const
	{
		return
		D1() * (
				( 1.0 + cosphis() )  * intExpL()     
			  + ( 1.0 - cosphis() )  * intExpH()
		) +
		D2() * (
				 ( 2.0 * sinphis() ) * intExpSin( ) 
			  +  ( 2.0 * CC()      ) * intExpCos( ) 
		) ;
	}
	
	
	//..................................
	inline double timeFactorOdd(  )   const
	{
		//if( t < 0.0 ) return 0.0 ;
		return
		D1() * (
				( 1.0 - cosphis() ) * expL( ) 
			  + ( 1.0 + cosphis() ) * expH( ) 
		) +
		D2() * (
				-  ( 2.0 * sinphis() ) * expSin( ) 
				+  ( 2.0 * CC()      ) * expCos( )  
		) ;
	}
	
	inline double timeFactorOddInt(  )  const
	{
		return
		D1() * (
				( 1.0 - cosphis() ) * intExpL()
			  + ( 1.0 + cosphis() ) * intExpH() 
		) +
		D2() * (
				-  ( 2.0 * sinphis() ) * intExpSin( ) 
				+  ( 2.0 * CC()      ) * intExpCos( ) 
		) ;
	}
	
	
	//----------------------------------------------------------
	// These are the time factors and their analytic integrals for the three angle PDF
	
	//...........................
	inline double timeFactorA0A0( )    const { return timeFactorEven( ) ; }     
	inline double timeFactorA0A0Int( ) const { return timeFactorEvenInt( ) ; }
	
	//...........................
	inline double timeFactorAPAP( )    const { return timeFactorEven( ) ; }
	inline double timeFactorAPAPInt( ) const { return timeFactorEvenInt( ) ; }
	
	//...........................
	inline double timeFactorATAT( )    const { return timeFactorOdd( ) ; }
	inline double timeFactorATATInt( ) const { return timeFactorOddInt( ) ; }
	
	//...........................
	inline double timeFactorImAPAT( ) const
	{
		return
		D1() * (
				  ( expL( ) - expH( ) ) * cos(delta1) * sinphis()  
				+ ( expL( ) + expH( ) ) * sin(delta1) * CC()  
		) +
		D2() * (
				 2.0  * ( sin(delta1)*expCos( ) - cos(delta1)*cosphis()*expSin( ) ) 
		) ;
	}
	
	inline double timeFactorImAPATInt( ) const
	{
		
		return
		D1() * (
				  ( intExpL() - intExpH() ) * cos(delta1) * sinphis() 	
				+ ( intExpL() + intExpH() ) * sin(delta1) * CC() 
		) +
		D2() * (
				 2.0  * ( sin(delta1)*intExpCos() - cos(delta1)*cosphis()*intExpSin() )
		) ;
	}
	
	
	//...........................
	inline double timeFactorReA0AP( )  const
	{
		if( _useCosDpar ) return cosdpar * this->timeFactorEven(  ) ;//PELC-COSDPAR Special for fitting cosdpar separately
		else return cos(delta2-delta1) * this->timeFactorEven(  ) ;
	}
	
	inline double timeFactorReA0APInt( ) const
	{
		if( _useCosDpar ) return cosdpar * this->timeFactorEvenInt( ) ;//PELC-COSDPAR Special for fitting cosdpar separately
		else return cos(delta2-delta1) * this->timeFactorEvenInt( ) ;
	}
	
	
	//...........................
	inline double timeFactorImA0AT(  ) const
	{
		return 
		D1() * (
				  ( expL( ) - expH( ) ) * cos(delta2) * sinphis() 
				+ ( expL( ) + expH( ) ) * sin(delta2) * CC() 
		) +
		D2() * (
				 2.0  * ( sin(delta2)*expCos( ) - cos(delta2)*cosphis()*expSin( ) ) 
		) ;
	}
	
	inline double timeFactorImA0ATInt( ) const
	{
		
		return
		D1() * (
				  ( intExpL() - intExpH()  ) * cos(delta2) * sinphis() 
				+ ( intExpL() + intExpH()  ) * sin(delta2) * CC() 
		) +
		D2() * (
				 2.0  * ( sin(delta2)*intExpCos() - cos(delta2)*cosphis()*intExpSin()  ) 
		) ;
	}
	
	//.... S wave additions.......
	
	//...........................
	inline double timeFactorASAS( )    const { return timeFactorOdd( ) ; }
	inline double timeFactorASASInt( ) const { return timeFactorOddInt( ) ; }
	
	
	//...........................
	inline double timeFactorReASAP( ) const
	{
		double delta = delta_para - delta_s ;
		return
		D1() * (
				  ( expL( ) - expH( ) ) * sin(delta) * sinphis()  
				+ ( expL( ) + expH( ) ) * cos(delta) * CC()  
		) +
		D2() * (
				 2.0  * ( cos(delta)*expCos( ) - sin(delta)*cosphis()*expSin( ) ) 
		) ;
	}
	
	inline double timeFactorReASAPInt( ) const
	{
		
		double delta = delta_para - delta_s ;
		
		return
		D1() * (	
				  ( intExpL() - intExpH() ) * sin(delta) * sinphis() 	    
				+ ( intExpL() + intExpH() ) * cos(delta) * CC() 
		) +
		D2() * (
				 2.0  * ( cos(delta)*intExpCos() - sin(delta)*cosphis()*intExpSin() ) 
		) ;
	}
	
	
	//...........................
	inline double timeFactorImASAT( )  const
	{
		return sin(delta_perp-delta_s) * this->timeFactorOdd(  ) ;
	}
	
	inline double timeFactorImASATInt( ) const
	{
		return sin(delta_perp-delta_s) * this->timeFactorOddInt( ) ;
	}
	
	
	//...........................
	inline double timeFactorReASA0( ) const
	{			
		double delta = delta_zero - delta_s ;
		return
		D1() * (
				  ( expL( ) - expH( ) ) * sin(delta) * sinphis()  
				+ ( expL( ) + expH( ) ) * cos(delta) * CC()  
		) +
		D2() * (
				 2.0  * ( cos(delta)*expCos( ) - sin(delta)*cosphis()*expSin( ) ) 
		) ;
	}
	
	inline double timeFactorReASA0Int( ) const
	{
		
		double delta = delta_zero - delta_s ;
		
		return
		D1() * (
				  ( intExpL() - intExpH() ) * sin(delta) * sinphis() 	    
				+ ( intExpL() + intExpH() ) * cos(delta) * CC() 
		) +
		D2() * (
				 2.0  * ( cos(delta)*intExpCos() - sin(delta)*cosphis()*intExpSin() ) 
		) ;
	}
	
	
	//------------------------------------------------------
	// Angle factors for three angle PDFs  - with transversity/helicity switch
	
	inline double angleFactorEven(  )   const { return _useHelicityBasis ?	HangleFactorEven(  )	: TangleFactorEven(  ) ; }
	inline double angleFactorOdd(  )	const { return _useHelicityBasis ?	HangleFactorOdd(  )		: TangleFactorOdd(  ) ; }

	
	inline double angleFactorA0A0( )	const { return _useHelicityBasis ?  HangleFactorA0A0(  )	: TangleFactorA0A0(  ) ;	}
	inline double angleFactorAPAP( )    const { return _useHelicityBasis ?  HangleFactorAPAP(  )	: TangleFactorAPAP(  ) ;	}
	inline double angleFactorATAT( )	const { return _useHelicityBasis ?  HangleFactorATAT(  )	: TangleFactorATAT(  ) ;	}
	inline double angleFactorImAPAT( )	const { return _useHelicityBasis ?  HangleFactorImAPAT(  )	: TangleFactorImAPAT(  ) ;	}
	inline double angleFactorReA0AP( )  const { return _useHelicityBasis ?  HangleFactorReA0AP(  )	: TangleFactorReA0AP(  ) ;	}
	inline double angleFactorImA0AT( )	const { return _useHelicityBasis ?  HangleFactorImA0AT(  )	: TangleFactorImA0AT(  ) ;	}
	inline double angleFactorASAS( )	const { return _useHelicityBasis ?  HangleFactorASAS(  )	: TangleFactorASAS(  ) ;	}
	inline double angleFactorReASAP( )	const {	return _useHelicityBasis ?  HangleFactorReASAP(  )	: TangleFactorReASAP(  ) ;	}
	inline double angleFactorImASAT( )	const { return _useHelicityBasis ?  HangleFactorImASAT(  )	: TangleFactorImASAT(  ) ;	}
	inline double angleFactorReASA0( )	const { return _useHelicityBasis ?  HangleFactorReASA0(  )	: TangleFactorReASA0(  ) ;	}

/*	This is for comparing transversity and helicity when needed.  To do this you have to also make sure both sets are declared and read in in the pdf - so its not a trivial change.
	inline double angleFactorA0A0( )	const { { cout << " AF A0A0 T   = "<<TangleFactorA0A0(  )<<  "  /  H = "<<HangleFactorA0A0(  ) << endl; } return _useHelicityBasis ?  HangleFactorA0A0(  )	: TangleFactorA0A0(  ) ;	}
	inline double angleFactorAPAP( )    const { { cout << "    APAP T   = "<<TangleFactorAPAP(  )<<  "  /  H = "<<HangleFactorAPAP(  ) << endl; }return _useHelicityBasis ?  HangleFactorAPAP(  )	: TangleFactorAPAP(  ) ;	}
	inline double angleFactorATAT( )	const { { cout << "    ATAT T   = "<<TangleFactorATAT(  )<<  "  /  H = "<<HangleFactorATAT(  ) << endl; }return _useHelicityBasis ?  HangleFactorATAT(  )	: TangleFactorATAT(  ) ;	}
	inline double angleFactorImAPAT( )	const { { cout << "    ImAPAT T = "<<TangleFactorImAPAT(  )<<  "  /  H = "<<HangleFactorImAPAT(  ) << endl; }return _useHelicityBasis ?  HangleFactorImAPAT(  )	: TangleFactorImAPAT(  ) ;	}
	inline double angleFactorReA0AP( )  const { { cout << "    ReA0AP T = "<<TangleFactorReA0AP(  )<<  "  /  H = "<<HangleFactorReA0AP(  ) << endl; }return _useHelicityBasis ?  HangleFactorReA0AP(  )	: TangleFactorReA0AP(  ) ;	}
	inline double angleFactorImA0AT( )	const { { cout << "    ImA0AT T = "<<TangleFactorImA0AT(  )<<  "  /  H = "<<HangleFactorImA0AT(  ) << endl; }return _useHelicityBasis ?  HangleFactorImA0AT(  )	: TangleFactorImA0AT(  ) ;	}
	inline double angleFactorASAS( )	const { { cout << "    ASAS T   = "<<TangleFactorASAS(  )<<  "  /  H = "<<HangleFactorASAS(  ) << endl; }return _useHelicityBasis ?  HangleFactorASAS(  )	: TangleFactorASAS(  ) ;	}
	inline double angleFactorReASAP( )	const {	{ cout << "    ReASAP T = "<<TangleFactorReASAP(  )<<  "  /  H = "<<HangleFactorReASAP(  ) << endl; }return _useHelicityBasis ?  HangleFactorReASAP(  )	: TangleFactorReASAP(  ) ;	}
	inline double angleFactorImASAT( )	const { { cout << "    ImASAT T = "<<TangleFactorImASAT(  )<<  "  /  H = "<<HangleFactorImASAT(  ) << endl; }return _useHelicityBasis ?  HangleFactorImASAT(  )	: TangleFactorImASAT(  ) ;	}
	inline double angleFactorReASA0( )	const { { cout << "    ReASA0 T = "<<TangleFactorReASA0(  )<<  "  /  H = "<<HangleFactorReASA0(  ) << endl; }return _useHelicityBasis ?  HangleFactorReASA0(  )	: TangleFactorReASA0(  ) ;	}
*/	
	
	//------------------------------------------------------
	// Angle factors for three angle PDFs  in transversity basis
	
	//........ P Wave ..........
	
	//........ for one anfgle tests ...................
	inline double TangleFactorEven(  )   const { return (1.0 + ctrsq()) * Mathematics::Global_Frac()  * 4.0/3.0*TMath::Pi() ; }
	inline double TangleFactorOdd(  )   const { return  strsq() * Mathematics::Global_Frac() * 8.0/3.0*TMath::Pi(); }
	
	
	//...........................
	inline double TangleFactorA0A0(  )   const { return 2.0 * ct1sq() * (1.0 - strsq()*cphsq() ) * Mathematics::Global_Frac(); }
	
	//...........................
	inline double TangleFactorAPAP( )    const { return  st1sq() * (1.0 - strsq()*sphsq() ) * Mathematics::Global_Frac(); }
	
	//...........................
	inline double TangleFactorATAT(  )   const { return st1sq() * strsq() * Mathematics::Global_Frac(); }
	
	//...........................
	inline double TangleFactorImAPAT(  ) const { return  -1. * st1sq() * s2tr() * sph() * Mathematics::Global_Frac(); }
	
	//...........................
	inline double TangleFactorReA0AP( )  const { return   Mathematics::_Over_SQRT_2() * s2t1() * strsq() * s2ph() * Mathematics::Global_Frac(); }
	
	//...........................
	inline double TangleFactorImA0AT(  ) const { return   Mathematics::_Over_SQRT_2() * s2t1() * s2tr() * cph() * Mathematics::Global_Frac(); }
	
	//......  S wave  ....
	
	//.............................
	inline double TangleFactorASAS(  ) const   { return  2.0*Mathematics::Third() * (1.0 - strsq()*cphsq() ) * Mathematics::Global_Frac(); }
	
	//...........................
	inline double TangleFactorReASAP(  ) const {return Mathematics::Root_6()*Mathematics::Third() * st1() * strsq() * s2ph() *  Mathematics::Global_Frac(); }
	
	//...........................
	// This appreas to be ifferent to the LHCB note, but on inspection it is not. It is the difference in sign of ImASAT <=> ImATAS
	inline double TangleFactorImASAT(  ) const { return  Mathematics::Root_6()*Mathematics::Third() *  st1() * s2tr() *  cph() *  Mathematics::Global_Frac(); }
	
	
	//...........................
	inline double TangleFactorReASA0(  ) const
	{
		//There was a -1.0 here
		//Then I looked in the LHCb note being drafted and found this disagreed with it
		//For now i have made it +1.0 to agree with draft LHCb note.
		// Since then ive proved that this set of signs is consistent by my "pdfvalue < 0 test"
		return  4.0*Mathematics::Root_3()*Mathematics::Third() * ct1() *  ( 1.0 - strsq()* cphsq() ) * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi()) ;
	}
	
	
	
	//------------------------------------------------------
	// Angle factors for three angle PDFs  in helicity basis
	
	//........ for one anfgle tests ...................
	inline double HangleFactorEven(  )  const { cout<<"No Helicity One Angle Formula Yet"<<endl; exit(1); return 0; }
	inline double HangleFactorOdd(  )   const { cout<<"No Helicity One Angle Formula Yet"<<endl; exit(1); return 0; }
	
	
	//........ P Wave ..........
	
	// .................
	inline double HangleFactorA0A0( ) const { return 4. * cHtksq() * sHtlsq() * Mathematics::Global_Frac() * 0.5; }
	
	//..................
	inline double HangleFactorAPAP( ) const { return ( sHtksq()*(1.0+cHtlsq()) - sHtlsq()*sHtksq()*c2Hphi() ) * Mathematics::Global_Frac() * 0.5; }
	
	//..................
	inline double HangleFactorATAT( ) const { return ( sHtksq()*(1.0+cHtlsq()) + sHtlsq()*sHtksq()*c2Hphi() ) * Mathematics::Global_Frac() * 0.5; }
	
	//..................
	inline double HangleFactorImAPAT( ) const { return 2. * sHtksq() * sHtlsq() * s2Hphi() * Mathematics::Global_Frac()* 0.5 ; }
	
	//..................
	inline double HangleFactorReA0AP( )  const { return -sqrt(2.) * s2Htk() * s2Htl() * cHphi() * Mathematics::Global_Frac()* 0.5 ; }  
	
	//..................
	inline double HangleFactorImA0AT( )  const { return  sqrt(2.) * s2Htk() * s2Htl() * sHphi() * Mathematics::Global_Frac()* 0.5 ; }  
	
	//......  S wave  ....
	
	//.............................
	inline double HangleFactorASAS(  ) const   { return  4.0/3.0 * sHtlsq() * Mathematics::Global_Frac() * 0.5 ; }
	
	//...........................
	inline double HangleFactorReASAP(  ) const {return -2.0/3.0 * sqrt(6.0) * s2Htl() * sHtk() * cHphi() * Mathematics::Global_Frac() * 0.5 ; }
	
	//...........................
	inline double HangleFactorImASAT(  ) const { return 2.0/3.0 * sqrt(6.0) * s2Htl() * sHtk() * sHphi() * Mathematics::Global_Frac()* 0.5 ; }
	
	//...........................
	inline double HangleFactorReASA0(  ) const { return  8.0/3.0*sqrt(3.0) * sHtlsq() * cHtk() * Mathematics::Global_Frac() * 0.5 ;	}
	
	

};

#endif
