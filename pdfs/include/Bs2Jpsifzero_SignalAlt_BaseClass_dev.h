// $Id: Bs2Jpsifzero_SignalAlt_BaseClass_dev.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2Jpsifzero_SignalAlt_BaseClass_dev.h
 *
 *  Base Class for Bs2Jpsifzero_SignalAlt....  PDFs
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-12
 *
 */

#ifndef Bs2Jpsifzero_SignalAlt_BaseClass_dev_H
#define Bs2Jpsifzero_SignalAlt_BaseClass_dev_H

#include "BasePDF.h"
#include "PDFConfigurator.h"
#include "SlicedAcceptance.h"

#include "Mathematics.h"

#include <iostream>
#include <cstdlib>
#include <float.h>

//PELC
#include "TH1.h"
#include "TCanvas.h"
//~PELC

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6
#define DEBUGFLAG true

//========================================================================================================================



class Bs2Jpsifzero_SignalAlt_BaseClass_dev
{
	public:
		Bs2Jpsifzero_SignalAlt_BaseClass_dev ( PDFConfigurator* );
		Bs2Jpsifzero_SignalAlt_BaseClass_dev ( const Bs2Jpsifzero_SignalAlt_BaseClass_dev& );
		virtual ~Bs2Jpsifzero_SignalAlt_BaseClass_dev();

	private:
		Bs2Jpsifzero_SignalAlt_BaseClass_dev& operator=( const Bs2Jpsifzero_SignalAlt_BaseClass_dev& );

	protected:

		//PELC For debugging purposes
		//TH1D * histOfPdfValues ;
		//TCanvas * c0 ; 
		//mutable int histCounter ;
		//~PELC
	
		// These contain the ObservableRefs that correspond to the physics parameter names and references
		ObservableRef gammaName;		// gamma
		ObservableRef deltaGammaName;	// delta gamma
		ObservableRef deltaMName;		// delta mass
		ObservableRef Phi_sName;		// what we want to measure!
		ObservableRef Azero_sqName;	// amplitude
		ObservableRef Apara_sqName;	// amplitude
		ObservableRef Aperp_sqName;	// amplitude
		ObservableRef As_sqName;		// amplitude
		ObservableRef delta_zeroName;	// strong phase, set to 0
		ObservableRef delta_paraName;	// strong phase
		ObservableRef delta_perpName;	// strong phase
		ObservableRef delta_sName;		// strong phase for S-wave

		//PELC NEW : These are new physics parameters which might be used instead of some of those above later
		ObservableRef cosphisName;		// fitting cosphis and sinphis independently
		ObservableRef sinphisName;		// fitting cosphis and sinphis independently	

		// Mistag parameters  
		ObservableRef mistagName;		// mistag fraction  - may be used as observable also
		ObservableRef mistagP1Name;		// mistag calib
		ObservableRef mistagP0Name;		// mistag calib
		ObservableRef mistagSetPointName;// mistag calib

		// Time resolution
		ObservableRef resScaleName;			// Scale to multiply all Gaussians with 
		ObservableRef res1Name;				// time resolution narrow
		ObservableRef res2Name;				// time resolution wide
		ObservableRef res3Name;				// time resolution tail
		ObservableRef res2FractionName;		// fraction of wide
		ObservableRef res3FractionName;		// fraction of tail
		ObservableRef timeOffsetName;		// time offset

		// These are the angular accceptance factors. The first 6 are P-wave, the second 4 are S-wave
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
		ObservableRef cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		ObservableRef phiName;			// azimuthal angle of the mu+ in Jpsi frame
		ObservableRef cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
		ObservableRef tagName;			// B tag

		//	Constraints
		ObservableRef timeConstraintName;  // what is this ????????????

		// Measured Event Observables
		double t ;
		double ctheta_tr ;
		double phi_tr ;
		double ctheta_1 ;
		int tag ;
		//X int timeAcceptanceCategory ;
	
		// Physics Fit Parameters 
		double _gamma ;
		double dgam ;

		double Aperp_sq ;
		double Apara_sq ;
		double Azero_sq ;
		double As_sq ;

		double delta_para ;
		double delta_perp ;
		double delta_zero ;
		double delta_s ;
		double delta1 ;
		double delta2 ;
	
		double delta_ms ;
		double phi_s ;
		double _cosphis ;
		double _sinphis ;
	
		// Mistag parameters
		double _mistag ;
		double _mistagP1 ;
		double _mistagP0 ;
		double _mistagSetPoint ;
	
		// Time resolution
		double resolution ;
		double resolutionScale ;
		double resolution1 ;
		double resolution2 ;
		double resolution3 ;
		double resolution2Fraction ;
		double resolution3Fraction ;
		double timeOffset ;

                // Angular acceptance factors
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
	
		//Time acceptance 
		SlicedAcceptance * timeAcc ;
	
		//Configurationparameters
		bool _useTimeAcceptance ;
	
		bool _numericIntegralForce ;
		bool _numericIntegralTimeOnly ;
	
		bool _useCosAndSin ;
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

		
		inline double gamma_l() const { 
			const double gl = gamma() + ( dgam *0.5 ) ;
			if( gl < 0. ) {
				cerr << " In Bs2Jpsifzero_SignalAlt_BaseClass_dev : gamma_l() < 0 so setting it to 0.0000001 " << endl ;
				return 0.0000001 ;
			}
			else
				return gl ; 
		}

		inline double gamma_h() const { 
			const double gh = gamma() - ( dgam *0.5 ) ;
			if( gh < 0. ) {
				cerr << " In Bs2Jpsifzero_SignalAlt_BaseClass_dev : gamma_h() < 0 so setting it to 0.0000001 " << endl ;
				return 0.0000001 ;
			}
			else
				return gh ;   
		}

		inline double gamma() const { return _gamma ; }

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
				cout << "Bs2Jpsifzero_SignalAlt_BaseClass_dev::mistag() : _mistag < 0 so set to 0 " << endl ;
				returnValue = 0 ;
			}
			else if( _mistag > 0.5 ) {
				cout << "Bs2Jpsifzero_SignalAlt_BaseClass_dev::mistag() : _mistag > 0.5 so set to 0.5 "  << endl ;
				returnValue = 0.5 ;
			}
			else {
				cout << "Bs2Jpsifzero_SignalAlt_BaseClass_dev::mistag() : WARNING ******If you got here you dont know what you are doing  "  << endl ;
				exit(1);
			}
			return returnValue ;			
		}
	
		inline double cosphis() const { return _cosphis ; }
		inline double sinphis() const { return _sinphis ; }
	
		inline bool useTimeAcceptance() const { return _useTimeAcceptance ; }		
	
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
		double diffXsecCompositeNorm1(  )  ;

		bool normalisationCacheValid ;
		double normalisationCacheValue[3] ;
		//double normalisationCacheValueRes2[3] ;
	
		void DebugPrintXsec( string , double ) const ;
		void DebugPrintNorm( string , double ) const ;
	
		
		//------------------------------------------------------------------------------
		// These are the time factors and their analytic integrals for the one angle PDF

		//..................................
		inline double timeFactorEven(  )  const
		{
			//if( t < 0.0 ) return 0.0 ;
			const double result = 
			( 1.0 + cosphis() ) * expL( ) 
			+ ( 1.0 - cosphis() ) * expH( ) 
			+ q() * ( 2.0 * sinphis()   ) * expSin( ) * (1.0 - 2.0*mistag()) ;
		  
			//DEBUG
			if( DEBUGFLAG && (result < 0) ) {
				cout << " Bs2Jpsifzero_SignalAlt_BaseClass_dev::timeFactorEven() : result < 0 " << endl ;
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
			( 1.0 + cosphis() )  * intExpL()     
			+ ( 1.0 - cosphis() )  * intExpH()          
			+ q() * ( 2.0 * sinphis()   ) * intExpSin( ) * (1.0 - 2.0*mistag()) ;
		}


		//..................................
		inline double timeFactorOdd(  )   const
		{
			//if( t < 0.0 ) return 0.0 ;
			return
			( 1.0 - cosphis() ) * expL( ) 
			+ ( 1.0 + cosphis() ) * expH( ) 
			- q() * ( 2.0 * sinphis()   ) * expSin( ) * (1.0 - 2.0*mistag()) ;
		}

		inline double timeFactorOddInt(  )  const
		{
			return
			( 1.0 - cosphis() ) * intExpL()
			+ ( 1.0 + cosphis() ) * intExpH() 
			- q() * ( 2.0 * sinphis()   ) * intExpSin( ) * (1.0 - 2.0*mistag()) ;
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
			q() * 2.0  * ( sin(delta1)*expCos( ) - cos(delta1)*cosphis()*expSin( ) ) * (1.0 - 2.0*mistag())
			- ( expH( ) - expL( ) ) * cos(delta1) * sinphis()  ;
		}
		
		inline double timeFactorImAPATInt( ) const
		{
			//double _tlo = tlo ;
			//if(_tlo < 0.) _tlo = 0. ;
			
			return
			q() * 2.0  * ( sin(delta1)*intExpCos() - cos(delta1)*cosphis()*intExpSin() ) * (1.0 - 2.0*mistag())
			- ( intExpH() - intExpL() ) * cos(delta1) * sinphis() ;	
		}


		//...........................
		inline double timeFactorReA0AP( )  const
		{
			return cos(delta2-delta1) * this->timeFactorEven(  ) ;
		}

		inline double timeFactorReA0APInt( ) const
		{
			return cos(delta2-delta1) * this->timeFactorEvenInt( ) ;
		}


		//...........................
		inline double timeFactorImA0AT(  ) const
		{
			return 
			q() * 2.0  * ( sin(delta2)*expCos( ) - cos(delta2)*cosphis()*expSin( ) ) * (1.0 - 2.0*mistag())	
			- ( expH( ) - expL( ) ) * cos(delta2) * sinphis() ;
		}

		inline double timeFactorImA0ATInt( ) const
		{
			//double _tlo = tlo ;
			//if(_tlo < 0.) _tlo = 0. ;
	
			return 
			q() * 2.0  * ( sin(delta2)*intExpCos() - cos(delta2)*cosphis()*intExpSin()  ) * (1.0 - 2.0*mistag())
			- ( intExpH() - intExpL()  ) * cos(delta2) * sinphis() ;
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
			q() * 2.0  * ( cos(delta)*expCos( ) - sin(delta)*cosphis()*expSin( ) ) * (1.0 - 2.0*mistag())
			- ( expH( ) - expL( ) ) * sin(delta) * sinphis()  ;
		}

		inline double timeFactorReASAPInt( ) const
		{
			//double _tlo = tlo ;
			//if(_tlo < 0.) _tlo = 0. ;

			double delta = delta_para - delta_s ;

			return
			q() * 2.0  * ( cos(delta)*intExpCos() - sin(delta)*cosphis()*intExpSin() ) * (1.0 - 2.0*mistag())
			- ( intExpH() - intExpL() ) * sin(delta) * sinphis() ;	    
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
			q() * 2.0  * ( cos(delta)*expCos( ) - sin(delta)*cosphis()*expSin( ) ) * (1.0 - 2.0*mistag())
			- ( expH( ) - expL( ) ) * sin(delta) * sinphis()  ;
		}

		inline double timeFactorReASA0Int( ) const
		{
			//double _tlo = tlo ;
			//if(_tlo < 0.) _tlo = 0. ;
			
			double delta = delta_zero - delta_s ;
			
			return
			q() * 2.0  * ( cos(delta)*intExpCos() - sin(delta)*cosphis()*intExpSin() ) * (1.0 - 2.0*mistag())
			- ( intExpH() - intExpL() ) * sin(delta) * sinphis() ;	    
		}


		

		//------------------------------------------------------
		// Angle factors for three angle PDFs

		//........ P Wave ..........

                //...........................
		inline double angleFactorA0A0(  )   const { return 2.0 * ct1sq() * (1.0 - strsq()*cphsq() ) * Mathematics::Global_Frac(); }

		//...........................
		inline double angleFactorAPAP( )    const { return  st1sq() * (1.0 - strsq()*sphsq() ) * Mathematics::Global_Frac(); }

		//...........................
		inline double angleFactorATAT(  )   const { return st1sq() * strsq() * Mathematics::Global_Frac(); }

		//...........................
		inline double angleFactorImAPAT(  ) const { return  -1. * st1sq() * s2tr() * sph() * Mathematics::Global_Frac(); }

		//...........................
		inline double angleFactorReA0AP( )  const { return   Mathematics::_Over_SQRT_2() * s2t1() * strsq() * s2ph() * Mathematics::Global_Frac(); }

		//...........................
		inline double angleFactorImA0AT(  ) const { return   Mathematics::_Over_SQRT_2() * s2t1() * s2tr() * cph() * Mathematics::Global_Frac(); }

		//......  S wave additions ....

		//.............................
		inline double angleFactorASAS(  ) const   { return  2.0*Mathematics::Third() * (1.0 - strsq()*cphsq() ) * Mathematics::Global_Frac(); }

		//...........................
		inline double angleFactorReASAP(  ) const {return Mathematics::Root_6()*Mathematics::Third() * st1() * strsq() * s2ph() *  Mathematics::Global_Frac(); }

		//...........................
		// This appreas to be ifferent to the LHCB note, but on inspection it is not. It is the difference in sign of ImASAT <=> ImATAS
		inline double angleFactorImASAT(  ) const { return  Mathematics::Root_6()*Mathematics::Third() *  st1() * s2tr() *  cph() *  Mathematics::Global_Frac(); }


		//...........................
		inline double angleFactorReASA0(  ) const
		{
			//There was a -1.0 here
			//Then I looked in the LHCb note being drafted and found this disagreed with it
			//For now i have made it +1.0 to agree with draft LHCb note.
			// Since then ive proved that this set of signs is consistent by my "pdfvalue < 0 test"
			return  4.0*Mathematics::Root_3()*Mathematics::Third() * ct1() *  ( 1.0 - strsq()* cphsq() ) * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi()) ;
		}

};

#endif
