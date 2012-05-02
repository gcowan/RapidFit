// $Id: Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4.h
 *
 *  Base Class for Bs2JpsiPhi_SignalAlt....  PDFs
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-12
 *
 */

#ifndef Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4_H
#define Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4_H

#include "BasePDF.h"
#include "PDFConfigurator.h"
#include "SlicedAcceptance.h"

#include "Mathematics.h"

#include <iostream>
#include <cstdlib>
#include <float.h>
#include <vector>

//PELC
#include "TH1.h"
#include "TCanvas.h"
//~PELC

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6
#define DEBUGFLAG true

//========================================================================================================================



class Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4  :  public BasePDF 
{
	public:
		Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4(PDFConfigurator* );
		virtual ~Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4();
		virtual bool SetPhysicsParameters(ParameterSet*);

	protected:

	        //      Uncopyable!
		//Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4 ( const Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4& );
		//Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4& operator = ( const Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4& );

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
		ObservableRef Aeven_sqName;	// amplitude
		ObservableRef Aodd_sqName;	// amplitude
		ObservableRef As_sqName;		// amplitude

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

		// These are the angular accceptance factors. 
		ObservableRef angAccEvenName ;  
		ObservableRef angAccOddName ;

		// Observables 
		ObservableRef timeName;		// proper time
		ObservableRef cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		ObservableRef tagName;			// B tag

		//	Constraints
		ObservableRef timeConstraintName;  // what is this ????????????

		// Measured Event Observables
		double t ;
		double ctheta_tr ;
		int tag ;
	
		// Physics Fit Parameters 
		double _gamma ;
		double dgam ;

		double Aeven_sq ;
		double Aodd_sq ;
		double As_sq ;
	
		double CachedAEven ;
		double CachedAOdd ;
		double CachedAs ;
		void CacheAmplitudesAndAngles() ;
	
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
		double angAccEven ;
		double angAccOdd ;
	
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
		bool allowNegativeAsSq ;
	
		//....................................
		//Internal helper functions
	
		inline double ctrsq() const { return (ctheta_tr*ctheta_tr) ; }
		inline double strsq() const { return (1.0 - ctrsq()) ; }

		
		inline double gamma_l() const { 
			const double gl = gamma() + ( dgam *0.5 ) ;
			if( gl < 0. ) {
				cerr << " In Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4 : gamma_l() < 0 so setting it to 0.0000001 " << endl ;
				return 0.0000001 ;
			}
			else
				return gl ; 
		}

		inline double gamma_h() const { 
			const double gh = gamma() - ( dgam *0.5 ) ;
			if( gh < 0. ) {
				cerr << " In Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4 : gamma_h() < 0 so setting it to 0.0000001 " << endl ;
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
				cout << "Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::mistag() : _mistag < 0 so set to 0 " << endl ;
				returnValue = 0 ;
			}
			else if( _mistag > 0.5 ) {
				cout << "Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::mistag() : _mistag > 0.5 so set to 0.5 "  << endl ;
				returnValue = 0.5 ;
			}
			else {
				cout << "Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::mistag() : WARNING ******If you got here you dont know what you are doing  "  << endl ;
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
		double diffXsecCompositeNorm1( unsigned int resolutionIndex )  ;

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
				cout << " Bs2JpsiPhi_SignalAlt_BaseClass_1angle_v4::timeFactorEven() : result < 0 " << endl ;
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




		//------------------------------------------------------
		// Angle factors

		//........ P Wave ..........

		//........ for one anfgle tests ...................
		inline double angleFactorEven(  )   const { return (1.0 + ctrsq()) * 3.0/8.0 /*(Mathematics::Global_Frac()  * 4.0/3.0*TMath::Pi()*/ ; }
	
		//......... for one angle tests ...................
		inline double angleFactorOdd(  )   const { return  strsq() * 3.0/4.0 /*Mathematics::Global_Frac() * 8.0/3.0*TMath::Pi()*/ ; }

};

#endif
