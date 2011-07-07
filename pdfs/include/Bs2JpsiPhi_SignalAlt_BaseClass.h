// $Id: Bs2JpsiPhi_SignalAlt_BaseClass.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_BaseClass.h
 *
 *  Base Class for Bs2JpsiPhi_SignalAlt....  PDFs
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-12
 *
 */

#ifndef Bs2JpsiPhi_SignalAlt_BaseClass_H
#define Bs2JpsiPhi_SignalAlt_BaseClass_H

#ifndef __CINT__
#include "BasePDF.h"
#include "PDFConfigurator.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#include "framework/include/PDFConfigurator.h"
#endif

#ifndef __CINT__
#include "Mathematics.h"
#endif
#ifdef __CINT__
#include "framework/include/Mathematics.h"
#endif

#include "RooComplex.h"

#include <iostream>
#include <cstdlib>

//PELC
#include <TH1D.h>
#include <TCanvas.h>
//~PELC

#include "float.h"

#define DOUBLE_TOLERANCE DBL_MIN
#define DEBUGFLAG true

// Some hard code switches - 
#define USE_LOWER_TIME_ACCEPTANCE true   // If true then lower time acceptance is used **IF* the event is ALSO marked as time biased 
#define USE_UPPER_TIME_ACCEPTANCE false   // If true then upper time acceptance is used  

#define UPPER_TIME_ACCEPTANCE_FACTOR 0.025


//================================================================================================
//--------------------------------------------
// This is included here to avoid checking it in as a separate file to SVN - as this may
// not be the final solution

class TimeAcceptanceFunction {
	
	private :
	
	vector<double> sliceTimeStart ;
	vector<double> sliceAcceptance ;
	vector<double> sliceFraction ;
	
	int nslices ;
	int startSlice;
	
public:
	
	//.................................
	// Constructor
	TimeAcceptanceFunction() :  sliceTimeStart(), sliceAcceptance(), sliceFraction(), nslices(), startSlice() {
		
		// The data ....
		
		//nslices = 10 ;
		//double ts[10] = {0.05,  0.15,  0.25,  0.40,   0.55,  0.70,   0.90,   1.25,  1.75,  8.00 };
		//double ac[10] = {0.001, 0.040, 0.176, 0.497, 0.732, 0.861, 0.933, 0.976, 0.993, 1.0 };
		
        nslices = 31;
		
		double ts[31] = {0, 0.05,  0.1, 0.15, 0.2, 0.25, 0.3, 0.34, 0.38, 0.42,
						0.46, 0.5, 0.54, 0.58, 0.62, 0.66, 0.7,  0.75,  0.8, 0.85 , 
						0.9, 0.95, 1.0, 1.05, 1.1,  1.2, 1.3, 1.4,1.6 , 1.8, 
						2};
		//  MC
//		double ac[31] = {0.000116891 , 0.00515391 , 0.0288208 , 0.0860796 , 0.181686 , 0.30586 , 0.425054 , 0.525121 , 
//	                             0.613731 , 0.687362 , 0.747977 , 0.79574 , 0.834891 , 0.865584 , 0.888149 , 0.90725 , 0.924371 , 
//	                             0.938193 , 0.95011 , 0.960408 , 0.966892 , 0.97265 , 0.975763 , 0.979905 , 0.98415 , 0.988728 , 
//	                             0.990672 , 0.994114 , 0.996717 , 0.999169 , 1 };

		//  DATA
		double ac[31] = {0.000122069 , 0.00436829 , 0.022288 , 0.0633692 , 0.132094 , 0.225401 , 0.321979 , 0.409769 , 
			0.494438 , 0.570185 , 0.637997 , 0.695043 , 0.744576 , 0.784799 , 0.816805 , 0.844606 , 0.869917 , 
			0.891761 , 0.910512 , 0.926418 , 0.937439 , 0.947105 , 0.953599 , 0.960726 , 0.968278 , 0.976149 , 
			0.980624 , 0.98685 , 0.991765 , 0.995575 , 1 };
		
		
		//double ac[10] = {1.0, 1.000, 1.00, 1.00, 1.00, 1.00, 1.000, 1.00, 1.00, 1.00 };
		
		//..............
		
		double previous = 0 ;
		
		for( int is=0; is < nslices ; ++is)
		{
			sliceTimeStart.push_back( ts[is] ) ;
			sliceAcceptance.push_back( ac[is] ) ;
			sliceFraction.push_back(ac[is]-previous) ;
			previous = ac[is] ;
		}
		
		startSlice = -1 ;    //Uninitialised
		
	}
	
	//.................................
	void configure ( double tstart ) 
	{
		if( tstart < sliceTimeStart[0] ) {cout << " InTimeAcceptance::configure : start time less than bins specified " << endl ; exit(1) ;}	
		
		// If in last time slice, or later, then slice = nslices-1
		if ( sliceTimeStart[nslices-1] <= tstart ) { startSlice = nslices-1 ; }
		
		// Find the slice which corresponds to the lower chosen limit of integration.
		else {
			for( int is = 0; is < nslices-1; ++is ) { 
				if( (sliceTimeStart[is] <= tstart ) && ( tstart < sliceTimeStart[is+1] )  )  
				{
					startSlice = is ;
					break ;
				}
			}
		}
		
		if( fabs(startSlice + 1. ) < DOUBLE_TOLERANCE )  {cout << " InTimeAcceptance::configure : times in slices are screwed up " << endl ; exit(1) ;}	
		
		return ;
	}
	
	//.................................
	// Return the acceptance corresponding to this slice
	double acceptance ( double time ) 
	{
	        
		// If less than first time slice then zero
		if( time < sliceTimeStart[0] ) { return 0 ; }
		
		// If in last time slice, or later, then slice = last slice = nslices-1
		else if ( sliceTimeStart[nslices-1] <= time ) { return sliceAcceptance[nslices-1] ; }
		
		// Find slice in range
		else {
			for( int is = 0; is < nslices-1; ++is ) { 
				if( (sliceTimeStart[is] <= time ) && ( time < sliceTimeStart[is+1] )  ) return sliceAcceptance[is];  
			}
		}
		
		cout << " InTimeAcceptance::acceptance : times in slices are screwed up - couldnt find slice - time requested was : " << time << endl ; 
		exit(1) ;	
		
		return 0 ;
	}
	
	
	//.....................................
	int numberOfSlices() { return nslices - startSlice ; }
	int firstSlice( ) { return startSlice ; }
	int lastSlice( ) { return nslices-1 ; }
	
	//.....................................
	// Return time of beginning of selected slice -  must be between startSlice and nslices-1
	double sliceStart( int islice ) {
		if( (islice < startSlice) || (islice >= nslices ) ) {cout << " InTimeAcceptance::sliceFraction : islice is crazy: " << islice << endl ; exit(1) ;} 
		return sliceTimeStart[islice] ;
	}
	
	
	//.....................................
	//The fraction to associate with this slice
	double fraction( int islice ) {
		
		if( (islice < startSlice) || (islice >= nslices ) ) {cout << " InTimeAcceptance::sliceFraction : islice is crazy: " << islice << endl ; exit(1) ;} 
		
		// If this is first slice, then you need to do something differently - yuo need cumulative fraction (i.e. acceptance)
		if( islice == startSlice ) {
			return sliceAcceptance[islice] ;
		}
		else{
			return sliceFraction[islice] ; 
		}
	}
	
	//.....................
	void print() {
		
		if( startSlice == -1 ) {
			cout << " InTimeAcceptance::print() - cant call this before configuring " << endl ;
			return ;
		}
		
		cout << endl ;
		cout  << "Time Acceptane Histogram " << endl ;
		for( int is=this->firstSlice() ; is <= this->lastSlice() ; ++is ) {
			cout << "    Time slice : "  << this->sliceStart(is) << "  Fraction : " << this->fraction(is) << endl ;
		}
		return ;
	}
	
};


//========================================================================================================================
//========================================================================================================================



class Bs2JpsiPhi_SignalAlt_BaseClass
{
	public:
		Bs2JpsiPhi_SignalAlt_BaseClass();
		Bs2JpsiPhi_SignalAlt_BaseClass(PDFConfigurator);
		virtual ~Bs2JpsiPhi_SignalAlt_BaseClass();

	protected:
	
		//PELC For debugging purposes
		//TH1D * histOfPdfValues ;
		//TCanvas * c0 ; 
		//mutable int histCounter ;
		//~PELC

		//Configurationparameters
		bool useCosAndSin ;
		bool allowNegativeAsSq ;
	
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
                ObservableRef mistagName;		// mistag fraction
		ObservableRef mistagP1Name;		// mistag calib
		ObservableRef mistagP0Name;		// mistag calib
		ObservableRef mistagSetPointName;// mistag calib

		// Time resolution
		ObservableRef res1Name;		  // time resolution core
		ObservableRef res2Name;		  // time resolution tail
		ObservableRef res1FractionName;  // fraction of core
		ObservableRef timeOffsetName;    // time offset

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

		ObservableRef timeName;		// proper time
		ObservableRef cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		ObservableRef phiName;			// azimuthal angle of the mu+ in Jpsi frame
		ObservableRef cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
		ObservableRef tagName;			// B tag
		ObservableRef timeAcceptanceCategoryName ; //Originally included for using unbiased and biased events.

		//	Constraints
		ObservableRef timeConstraintName;

		// Measured Event Observables
		double t ;
		double ctheta_tr ;
		double phi_tr ;
		double ctheta_1 ;
		int tag ;
		int timeAcceptanceCategory ;
	
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
		double resolution1 ;
		double resolution2 ;
		double resolution1Fraction ;
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
	
		// lower time acceptance object
		TimeAcceptanceFunction timeAcceptance ;

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
				cerr << " In Bs2JpsiPhi_SignalAlt_BaseClass : gamma_l() < 0 so setting it to 0.0000001 " << endl ;
				return 0.0000001 ;
			}
			else
				return gl ; 
		}

		inline double gamma_h() const { 
			const double gh = gamma() - ( dgam *0.5 ) ;
			if( gh < 0. ) {
				cerr << " In Bs2JpsiPhi_SignalAlt_BaseClass : gamma_h() < 0 so setting it to 0.0000001 " << endl ;
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
				cout << "Bs2JpsiPhi_SignalAlt_BaseClass::mistag() : _mistag < 0 so set to 0 " << endl ;
				returnValue = 0 ;
			}
			else if( _mistag > 0.5 ) {
				cout << "Bs2JpsiPhi_SignalAlt_BaseClass::mistag() : _mistag > 0.5 so set to 0.5 "  << endl ;
				returnValue = 0.5 ;
			}
			else {
				cout << "Bs2JpsiPhi_SignalAlt_BaseClass::mistag() : WARNING ******If you got here you dont know what you are doing  "  << endl ;
				exit(1);
			}
			return returnValue ;			
		}
	
		inline double cosphis() const { return _cosphis ; }
		inline double sinphis() const { return _sinphis ; }

		inline bool useLowerTimeAcceptance() const { return (USE_LOWER_TIME_ACCEPTANCE && (timeAcceptanceCategory > 0)) ; }

		inline bool useUpperTimeAcceptance() const { return (USE_UPPER_TIME_ACCEPTANCE && ( fabs(UPPER_TIME_ACCEPTANCE_FACTOR - 0.) > DOUBLE_TOLERANCE) ) ; }

		inline double upperTimeAcceptanceBeta() const { return UPPER_TIME_ACCEPTANCE_FACTOR ; }
		
	
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
		double diffXsecNorm1(  ) const ;
		double diffXsecCompositeNorm1(  )  ;
		double diffXsecNorm2(  ) const ;
	
	
		
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
				cout << " Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorEven() : result < 0 " << endl ;
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
