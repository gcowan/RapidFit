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

#include "BasePDF.h"
#include "RooComplex.h"

#include <iostream>
#include <cstdlib>

#define DOUBLE_TOLERANCE 1E-6

// Some hard code switches - 
#define USE_LOWER_TIME_ACCEPTANCE false   // If true then lower time acceptance is used **IF* the event is ALSO marked as time biased 
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
		//  Data
		double ac[31] = {0.000116891 , 0.00515391 , 0.0288208 , 0.0860796 , 0.181686 , 0.30586 , 0.425054 , 0.525121 , 
	                             0.613731 , 0.687362 , 0.747977 , 0.79574 , 0.834891 , 0.865584 , 0.888149 , 0.90725 , 0.924371 , 
	                             0.938193 , 0.95011 , 0.960408 , 0.966892 , 0.97265 , 0.975763 , 0.979905 , 0.98415 , 0.988728 , 
	                             0.990672 , 0.994114 , 0.996717 , 0.999169 , 1 };		

		//  MC
//		double ac[31] = {0.000122069 , 0.00436829 , 0.022288 , 0.0633692 , 0.132094 , 0.225401 , 0.321979 , 0.409769 , 
//			0.494438 , 0.570185 , 0.637997 , 0.695043 , 0.744576 , 0.784799 , 0.816805 , 0.844606 , 0.869917 , 
//			0.891761 , 0.910512 , 0.926418 , 0.937439 , 0.947105 , 0.953599 , 0.960726 , 0.968278 , 0.976149 , 
//			0.980624 , 0.98685 , 0.991765 , 0.995575 , 1 };
		
		
		//double ac[10] = {1.0, 1.000, 1.00, 1.00, 1.00, 1.00, 1.000, 1.00, 1.00, 1.00 };
		
		//..............
		
		double previous = 0 ;
		
		for( int is=0; is < nslices ; is++)
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
			for( int is = 0; is < nslices-1; is++ ) { 
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
			for( int is = 0; is < nslices-1; is++ ) { 
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
		for( int is=this->firstSlice() ; is <= this->lastSlice() ; is++ ) {
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
	
		// These contain the strings that correspond to the physics parameter names.
		string gammaName;		// gamma
		string deltaGammaName;	// delta gamma
		string deltaMName;		// delta mass
		string Phi_sName;		// what we want to measure!
		string Azero_sqName;	// amplitude
		string Apara_sqName;	// amplitude
		string Aperp_sqName;	// amplitude
		string As_sqName;		// amplitude
		string delta_zeroName;	// strong phase, set to 0
		string delta_paraName;	// strong phase
		string delta_perpName;	// strong phase
		string delta_sName;		// strong phase for S-wave

	    // These are detector parameters
	    string mistagName;		// mistag fraction
		string res1Name;		  // time resolution core
		string res2Name;		  // time resolution tail
		string res1FractionName;  // fraction of core
		string timeOffsetName;    // time offset
	
	    // These are the angular accceptance factors. The first 6 are P-wave, the second 4 are S-wave
	    string angAccI1Name ;  
		string angAccI2Name ;
		string angAccI3Name ;
		string angAccI4Name ;
		string angAccI5Name ;
		string angAccI6Name ;
		string angAccI7Name ;
		string angAccI8Name ;
		string angAccI9Name ;
		string angAccI10Name ;
	
		// These contain the strings that correspond to the observable names
		string timeName;		// proper time
		string cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		string phiName;			// azimuthal angle of the mu+ in Jpsi frame
		string cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
		string tagName;			// B tag
		string timeAcceptanceCategoryName ;  //Originally included for using unbiased and biased events.

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
		bool useLowerTimeAcceptance() const ;
		bool useUpperTimeAcceptance() const ;
		double upperTimeAcceptanceBeta() const ;

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
			    	
		//Amplitudes Used in three angle PDF
		double AT() const ;
		double AP() const ;
		double A0() const ;
		double AS() const ;
	
		double ctrsq() const ;
		double strsq() const ;
		double ct1sq() const ;
		double st1sq() const ;
		double cphsq() const ;
		double sphsq() const ;
	
		// Widths
		double gamma_l() const ;
		double gamma_h() const ;
		double gamma() const ;
	
		// Time primitives
		double expL() const ;
		double expH() const ;
		double expCos() const ;
		double expSin() const ;
		double intExpL() const ;
		double intExpH() const ;
		double intExpSin() const ;
		double intExpCos() const ;
	
	
		//--------------------
		// Tag category, i.e B, Bbar or untagged.
		double q() const ;
	
	
		//----------------------------------------------------------
		// These are the time factors and their analytic integrals 
	
		//..................................
		double timeFactorEven(  )  const ;
		double timeFactorEvenInt(  )  const ;
	
		//..................................
		double timeFactorOdd(  )   const ;
		double timeFactorOddInt(  )  const ;
	
		//.......... P Wave .........
		//...........................
		double timeFactorA0A0( ) const ;      
		double timeFactorA0A0Int( ) const ;
	
		//...........................
		double timeFactorAPAP( ) const ;
		double timeFactorAPAPInt( ) const ;
	
		//...........................
		double timeFactorATAT( ) const ;
		double timeFactorATATInt( ) const ;

		//...........................
		double timeFactorImAPAT( ) const ; 
		double timeFactorImAPATInt( ) const ;
	
		//...........................	
		double timeFactorReA0AP( )  const ;		
		double timeFactorReA0APInt( ) const ;
	
		//...........................
		double timeFactorImA0AT(  ) const ;
		double timeFactorImA0ATInt( ) const ;
    
		//.....  S Wave..............
		//...........................
		double timeFactorASAS( ) const ;
		double timeFactorASASInt( ) const ;
	
		//............................
		double timeFactorReASAP( ) const ;
		double timeFactorReASAPInt( ) const ;

		//............................
		double timeFactorImASAT( ) const ;
		double timeFactorImASATInt( ) const ;

		//............................
		double timeFactorReASA0( ) const ;
		double timeFactorReASA0Int( ) const ;
	
		//-----------------------------------------------------------------
		// Angle factors 
	
		double angleFactorA0A0(  ) const ;
		double angleFactorAPAP(  ) const ;
		double angleFactorATAT(  ) const ;
		double angleFactorImAPAT(  ) const ;
		double angleFactorReA0AP( ) const ;
		double angleFactorImA0AT(  ) const ;

		double angleFactorASAS(  ) const ;
		double angleFactorReASAP( ) const ;
		double angleFactorImASAT(  ) const ;
		double angleFactorReASA0( ) const ;
	
};

#endif
