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
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
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

#define DOUBLE_TOLERANCE 1E-6
#define DEBUGFLAG true

// Some hard code switches - 
#define USE_LOWER_TIME_ACCEPTANCE true   // If true then lower time acceptance is used **IF* the event is ALSO marked as time biased 
#define USE_UPPER_TIME_ACCEPTANCE false   // If true then upper time acceptance is used  

#define UPPER_TIME_ACCEPTANCE_FACTOR 0.025


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
	TimeAcceptanceFunction() {
		
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
		
		if( (startSlice + -1. ) < DOUBLE_TOLERANCE )  {cout << " InTimeAcceptance::configure : times in slices are screwed up " << endl ; exit(1) ;}	
		
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


//---------------------------------------------



class Bs2JpsiPhi_SignalAlt_BaseClass
{
	public:
		Bs2JpsiPhi_SignalAlt_BaseClass();
		~Bs2JpsiPhi_SignalAlt_BaseClass();

	protected:
	
		// These contain the pair<string,int>s that correspond to the physics parameter names and references
		pair<string,int> gammaName;		// gamma
		pair<string,int> deltaGammaName;	// delta gamma
		pair<string,int> deltaMName;		// delta mass
		pair<string,int> Phi_sName;		// what we want to measure!
		pair<string,int> Azero_sqName;	// amplitude
		pair<string,int> Apara_sqName;	// amplitude
		pair<string,int> Aperp_sqName;	// amplitude
		pair<string,int> As_sqName;		// amplitude
		pair<string,int> delta_zeroName;	// strong phase, set to 0
		pair<string,int> delta_paraName;	// strong phase
		pair<string,int> delta_perpName;	// strong phase
		pair<string,int> delta_sName;		// strong phase for S-wave

	    // These are detector parameters
	    pair<string,int> mistagName;		// mistag fraction
		pair<string,int> res1Name;		  // time resolution core
		pair<string,int> res2Name;		  // time resolution tail
		pair<string,int> res1FractionName;  // fraction of core
		pair<string,int> timeOffsetName;    // time offset
	
	    // These are the angular accceptance factors. The first 6 are P-wave, the second 4 are S-wave
	    pair<string,int> angAccI1Name ;  
		pair<string,int> angAccI2Name ;
		pair<string,int> angAccI3Name ;
		pair<string,int> angAccI4Name ;
		pair<string,int> angAccI5Name ;
		pair<string,int> angAccI6Name ;
		pair<string,int> angAccI7Name ;
		pair<string,int> angAccI8Name ;
		pair<string,int> angAccI9Name ;
		pair<string,int> angAccI10Name ;
	
		// These contain the pair<string,int>s that correspond to the observable names
		pair<string,int> timeName;		// proper time
		pair<string,int> cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		pair<string,int> phiName;			// azimuthal angle of the mu+ in Jpsi frame
		pair<string,int> cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
		pair<string,int> tagName;			// B tag
		pair<string,int> timeAcceptanceCategoryName ;  //Originally included for using unbiased and biased events.

		// Measured Event Observables
		double t ;
		double ctheta_tr ;
		double phi_tr ;
		double ctheta_1 ;
		int tag ;
		int timeAcceptanceCategory ;
	
		// Physics Fit Parameters 
		double gamma_in ;
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
	
		// Other experimental parameters
		double tagFraction ;
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
	
		// lower and upper time acceptance things
		TimeAcceptanceFunction timeAcceptance ;
//		bool useLowerTimeAcceptance() const ;
//		bool useUpperTimeAcceptance() const ;
//		double upperTimeAcceptanceBeta() const ;

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

		
		//	It has been suggested that inlining as much of the code as possible will force the compiler
		//	to generate better code.
		//	This is going to be tested with a profile run.
		
		//Amplitudes Used in three angle PDF
//		inline double AT() const ;
//		inline double AP() const ;
//		inline double A0() const ;
//		inline double AS() const ;
	
//		inline double ctrsq() const ;
//		inline double strsq() const ;
//		inline double ct1sq() const ;
//		inline double st1sq() const ;
//		inline double cphsq() const ;
//		inline double sphsq() const ;
	
		// Widths
//		inline double gamma_l() const ;
//		inline double gamma_h() const ;
//		inline double gamma() const ;
		
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
		inline double ct1sq() const { return (ctheta_1*ctheta_1) ; }
		inline double st1sq() const { return (1.0 - ct1sq()) ; }
		inline double cphsq() const { return (cos(phi_tr)*cos(phi_tr)) ; }
		inline double sphsq() const { return (1.0 - cphsq()) ; }

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

		inline double gamma() const { return gamma_in ; }

		inline double q() const { return tag ;}

		inline bool useLowerTimeAcceptance() const { return (USE_LOWER_TIME_ACCEPTANCE && (timeAcceptanceCategory > 0)) ; }

		inline bool useUpperTimeAcceptance() const { return (USE_UPPER_TIME_ACCEPTANCE && ( fabs(UPPER_TIME_ACCEPTANCE_FACTOR - 0) < DOUBLE_TOLERANCE) ) ; }

		inline double upperTimeAcceptanceBeta() const { return UPPER_TIME_ACCEPTANCE_FACTOR ; }

		
	
		// Time primitives
//		inline double expL() const ;
//		inline double expH() const ;
//		inline double expCos() const ;
//		inline double expSin() const ;
//		inline double intExpL() const ;
//		inline double intExpH() const ;
//		inline double intExpSin() const ;
//		inline double intExpCos() const ;

		inline double expL() const 
		{
			return expL_stored ; //Mathematics::Exp( t, gamma_l(), resolution ) ;
		}

		inline double expH() const 
		{
			return expH_stored ; //Mathematics::Exp( t, gamma_h(), resolution ) ;
		}

		inline double intExpL( ) const {
			return intExpL_stored ;
			//if( useUpperTimeAcceptance() ) return Mathematics::ExpInt_betaAcceptance( tlo, thi, gamma_l(), resolution, upperTimeAcceptanceBeta() )  ;
			//else return Mathematics::ExpInt( tlo, thi, gamma_l(), resolution )  ;
		}

		inline double intExpH( ) const {
			return intExpH_stored ;
			//if( useUpperTimeAcceptance() )return Mathematics::ExpInt_betaAcceptance( tlo, thi, gamma_h(), resolution, upperTimeAcceptanceBeta() )  ;
			//else return Mathematics::ExpInt( tlo, thi, gamma_h(), resolution )  ;
		}

		//......................................................
		// Exponential x sine  and cosine

		inline double expSin() const  
		{
			return expSin_stored ; // Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
		}

		inline double expCos() const 
		{
			return expCos_stored ; // Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
		}

		inline double intExpSin( ) const 
		{
			return intExpSin_stored ;   // Mathematics::ExpSinInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
		}

		// Integral of exp( - G * t ) * cos( dm * t )  
		inline double intExpCos( ) const 
		{
			return intExpCos_stored ;   // Mathematics::ExpCosInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
		}



	
		//--------------------
		// Tag category, i.e B, Bbar or untagged.
//		double q() const ;
	
	
		//----------------------------------------------------------
		// These are the time factors and their analytic integrals 
	
//		//..................................
//		inline double timeFactorEven(  )  const ;
//		inline double timeFactorEvenInt(  )  const ;
	
		//..................................
//		inline double timeFactorOdd(  )   const ;
//		inline double timeFactorOddInt(  )  const ;
	
		//.......... P Wave .........
		//...........................
//		inline double timeFactorA0A0( ) const ;      
//		inline double timeFactorA0A0Int( ) const ;
	
		//...........................
//		inline double timeFactorAPAP( ) const ;
//		inline double timeFactorAPAPInt( ) const ;
	
		//...........................
//		inline double timeFactorATAT( ) const ;
//		inline double timeFactorATATInt( ) const ;

		//...........................
//		inline double timeFactorImAPAT( ) const ; 
//		inline double timeFactorImAPATInt( ) const ;
	
		//...........................	
//		inline double timeFactorReA0AP( )  const ;		
//		inline double timeFactorReA0APInt( ) const ;
	
		//...........................
//		inline double timeFactorImA0AT(  ) const ;
//		inline double timeFactorImA0ATInt( ) const ;
    
		//.....  S Wave..............
		//...........................
//		inline double timeFactorASAS( ) const ;
//		inline double timeFactorASASInt( ) const ;
	
		//............................
//		inline double timeFactorReASAP( ) const ;
//		inline double timeFactorReASAPInt( ) const ;

		//............................
//		inline double timeFactorImASAT( ) const ;
//		inline double timeFactorImASATInt( ) const ;

		//............................
//		inline double timeFactorReASA0( ) const ;
//		inline double timeFactorReASA0Int( ) const ;
		
		
		//------------------------------------------------------------------------------
// These are the time factors and their analytic integrals for the one angle PDF

		//..................................
		inline double timeFactorEven(  )  const
		{
			//if( t < 0.0 ) return 0.0 ;
			const double result = 
			( 1.0 + cos(phi_s) ) * expL( ) 
			+ ( 1.0 - cos(phi_s) ) * expH( ) 
			+ q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) ;
		  
			//DEBUG
			if( DEBUGFLAG && (result < 0) ) {
				cout << " Bs2JpsiPhi_SignalAlt_BaseClass::timeFactorEven() : result < 0 " << endl ;
				cout << " ->term1 " << ( 1.0 + cos(phi_s) ) * expL( ) << endl ;
				cout << " ->term2 " << ( 1.0 - cos(phi_s) ) * expH( ) << endl ;
				cout << " ->term3 " << q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) << endl ;
				cout << "   -->sin(phis) "  << sin(phi_s) << endl ;
				cout << "   -->expSin    "  << expSin() << endl ;
				cout << "   -->tagFrac   "  << tagFraction << endl ;
				cout << "   -->delta_ms  "  << delta_ms << endl ;
			}
			return result ;
		}

		inline double timeFactorEvenInt(  )  const
		{
			return
			( 1.0 + cos(phi_s) )  * intExpL()     
			+ ( 1.0 - cos(phi_s) )  * intExpH()          
			+ q() * ( 2.0 * sin(phi_s)   ) * intExpSin( ) * (1.0 - 2.0*tagFraction) ;
		}


		//..................................
		inline double timeFactorOdd(  )   const
		{
			//if( t < 0.0 ) return 0.0 ;
			return
			( 1.0 - cos(phi_s) ) * expL( ) 
			+ ( 1.0 + cos(phi_s) ) * expH( ) 
			- q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) ;
		}

		inline double timeFactorOddInt(  )  const
		{
			return
			( 1.0 - cos(phi_s) ) * intExpL()
			+ ( 1.0 + cos(phi_s) ) * intExpH() 
			- q() * ( 2.0 * sin(phi_s)   ) * intExpSin( ) * (1.0 - 2.0*tagFraction) ;
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
			//if( t < 0.0 ) return 0.0 ;
			return
			q() * 2.0  * ( sin(delta1)*expCos( ) - cos(delta1)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)
			- /*1.0 **/ ( expH( ) - expL( ) ) * cos(delta1) * sin(phi_s)  ;
		}
		
		inline double timeFactorImAPATInt( ) const
		{
			double _tlo = tlo ;
			if(_tlo < 0.) _tlo = 0. ;
			
			return
			q() * 2.0  * ( sin(delta1)*intExpCos() - cos(delta1)*cos(phi_s)*intExpSin() ) * (1.0 - 		2.0*tagFraction)
			- /*1.0 **/ ( intExpH() - intExpL() ) * cos(delta1) * sin(phi_s) ;	    
		}


		//...........................
		inline double timeFactorReA0AP( )  const
		{
			//if( t < 0.0 ) return 0.0 ;
			return cos(delta2-delta1) * this->timeFactorEven(  ) ;
		}

		inline double timeFactorReA0APInt( ) const
		{
			return cos(delta2-delta1) * this->timeFactorEvenInt( ) ;
		}


		//...........................
		inline double timeFactorImA0AT(  ) const
		{
			//if( t < 0.0 ) return 0.0 ;
			return 
			q() * 2.0  * ( sin(delta2)*expCos( ) - cos(delta2)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)	
			-/*1.0 **/ ( expH( ) - expL( ) ) * cos(delta2) * sin(phi_s) ;
		}

		inline double timeFactorImA0ATInt( ) const
		{
			double _tlo = tlo ;
			if(_tlo < 0.) _tlo = 0. ;
	
			return 
			q() * 2.0  * ( sin(delta2)*intExpCos() - cos(delta2)*cos(phi_s)*intExpSin()  ) * (1.0 - 2.0*tagFraction)
			-/*1.0 **/ ( intExpH() - intExpL()  ) * cos(delta2) * sin(phi_s) ;
		}

		//.... S wave additions.......

		//...........................
		inline double timeFactorASAS( )    const { return timeFactorOdd( ) ; }
		inline double timeFactorASASInt( ) const { return timeFactorOddInt( ) ; }


		//...........................
		inline double timeFactorReASAP( ) const
		{
			//if( t < 0.0 ) return 0.0 ;
			
			double delta = delta_para - delta_s ;
			return
			q() * 2.0  * ( cos(delta)*expCos( ) - sin(delta)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)
			- /*1.0 **/ ( expH( ) - expL( ) ) * sin(delta) * sin(phi_s)  ;
		}

		inline double timeFactorReASAPInt( ) const
		{
			double _tlo = tlo ;
			if(_tlo < 0.) _tlo = 0. ;

			double delta = delta_para - delta_s ;

			return
			q() * 2.0  * ( cos(delta)*intExpCos() - sin(delta)*cos(phi_s)*intExpSin() ) * (1.0 - 2.0*tagFraction)
			- /*1.0 **/ ( intExpH() - intExpL() ) * sin(delta) * sin(phi_s) ;	    
		}


		//...........................
		inline double timeFactorImASAT( )  const
		{
			//if( t < 0.0 ) return 0.0 ;
			return sin(delta_perp-delta_s) * this->timeFactorOdd(  ) ;
		}

		inline double timeFactorImASATInt( ) const
		{
			return sin(delta_perp-delta_s) * this->timeFactorOddInt( ) ;
		}


		//...........................
		inline double timeFactorReASA0( ) const
		{
			//if( t < 0.0 ) return 0.0 ;
			
			double delta = delta_zero - delta_s ;
			return
			q() * 2.0  * ( cos(delta)*expCos( ) - sin(delta)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)
			- /*1.0 **/ ( expH( ) - expL( ) ) * sin(delta) * sin(phi_s)  ;
		}

		inline double timeFactorReASA0Int( ) const
		{
			double _tlo = tlo ;
			if(_tlo < 0.) _tlo = 0. ;
			
			double delta = delta_zero - delta_s ;
			
			return
			q() * 2.0  * ( cos(delta)*intExpCos() - sin(delta)*cos(phi_s)*intExpSin() ) * (1.0 - 2.0*tagFraction)
			- /*1.0 **/ ( intExpH() - intExpL() ) * sin(delta) * sin(phi_s) ;	    
		}


		
		
		
		
		
	
		//-----------------------------------------------------------------
		// Angle factors 
	
//		double angleFactorA0A0(  ) const ;
//		double angleFactorAPAP(  ) const ;
//		double angleFactorATAT(  ) const ;
//		double angleFactorImAPAT(  ) const ;
//		double angleFactorReA0AP( ) const ;
//		double angleFactorImA0AT(  ) const ;

//		double angleFactorASAS(  ) const ;
//		double angleFactorReASAP( ) const ;
//		double angleFactorImASAT(  ) const ;
//		double angleFactorReASA0( ) const ;


		//------------------------------------------------------
		// Angle factors for three angle PDFs

		//........ P Wave ..........

		//...........................
		inline double angleFactorA0A0(  ) const
		{
			// Normalised to  1	
			return 2.0 * ct1sq() * (1.0 - strsq()*cphsq() ) * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi());
		}

		//...........................
		inline double angleFactorAPAP(  ) const
		{
			// Normalised to  1
			return  st1sq() * (1.0 - strsq()*sphsq() ) * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi());
		}

		//...........................
		inline double angleFactorATAT(  ) const
		{
			// Normalised to  1
			return st1sq() * strsq() * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi());
		}

		//...........................
		inline double angleFactorImAPAT(  ) const
		{
			// Normalised to  0
			double theta_tr = acos(ctheta_tr) ;		
			return   -/*1.0 **/  st1sq() * sin(2.0*theta_tr) * sin(phi_tr) * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi()) ;
		}

		//...........................
		inline double angleFactorReA0AP( ) const
		{
			// Normalised to  0
			double theta_1 = acos(ctheta_1) ;	
			return    sin(2.0*theta_1) * strsq() * sin(2.0*phi_tr) *Mathematics::_Over_SQRT_2()/*/ sqrt(2.0)*/ * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi());
		}

		//...........................
		inline double angleFactorImA0AT(  ) const
		{
			// Normalised to  0
			double theta_tr = acos(ctheta_tr) ;		
			double theta_1 = acos(ctheta_1) ;		
			return  +/*1.0**/   sin(2.0*theta_1) * sin(2.0*theta_tr) * cos(phi_tr) *Mathematics::_Over_SQRT_2()/*/ sqrt(2.0)*/ * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi());
		}

		//......  S wave additions ....

		//.............................
		inline double angleFactorASAS(  ) const
		{
			return  (1.0 - strsq()*cphsq() ) * 2*Mathematics::Third()/*(2./3.)*/ * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi());
		}

		//...........................
		inline double angleFactorReASAP(  ) const
		{
			double stheta_1 =  sqrt(st1sq());		
			return   strsq() * stheta_1 * sin(2.0*phi_tr) * Mathematics::Root_6()*Mathematics::Third()/*(sqrt(6.)/3.)*/ * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi()) ;
		}

		//...........................
		inline double angleFactorImASAT(  ) const
		{
			double theta_tr = acos(ctheta_tr) ;		
			double stheta_1 =  sqrt(st1sq());		
			return -/*1.0**/  sin(2.0*theta_tr) * stheta_1 * cos(phi_tr) * Mathematics::Root_6()*Mathematics::Third()/*(sqrt(6.)/3.)*/ * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi()) ;
		}


		//...........................
		inline double angleFactorReASA0(  ) const
		{
			return -/*1.0 **/  ( 1.0 -  strsq()* cphsq() ) * ctheta_1 *  4*Mathematics::Third()/*(4.0*sqrt(3.)/3.)*/ * Mathematics::Global_Frac();//(9.0/32.0/TMath::Pi()) ;
		}




};

#endif
