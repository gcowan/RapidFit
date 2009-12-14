// $Id: Bs2JpsiPhiClassic_Alt.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2JpsiPhiClassic_Alt RaPDF_.h
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2009-12-06
 */

#ifndef Bs2JpsiPhiClassic_Alt_H
#define Bs2JpsiPhiClassic_Alt_H

#include "BasePDF.h"
#include "RooComplex.h"

class Bs2JpsiPhiClassic_Alt : public BasePDF
{
	public:
		Bs2JpsiPhiClassic_Alt();
		~Bs2JpsiPhiClassic_Alt();

		//Mandatory method to evaluate the PDF value:
		virtual double Evaluate(DataPoint*);

		//Other opeating methods
		virtual bool SetPhysicsParameters(ParameterSet*);
		virtual vector<string> GetDoNotIntegrateList();

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
	
		void MakePrototypes();
		void getPhysicsParameters( double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double & );

		// These contain the strings that correspond to the physics parameter names.
		string gammaName;		// gamma
		string deltaGammaName;	// delta gamma
		string deltaMName;		// delta mass
		string Phi_sName;		// what we want to measure!
		string Azero_sqName;	// amplitude
		string Apara_sqName;	// amplitude
		string Aperp_sqName;	// amplitude
		string delta_zeroName;	// strong phase, set to 0
		string delta_paraName;	// strong phase
		string delta_perpName;	// strong phase
		string mistagName;		// mistag fraction
		string timeresName;		// time resolution
   
		// These contain the strings that correspond to the observable names
		string timeName;		// proper time
		string cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		string phiName;			// azimuthal angle of the mu+ in Jpsi frame
		string cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
		string tagName;			// B tag

		// Measured Event Observables
		double t ;
		double ctheta_tr ;
		double phi_tr ;
		double ctheta_1 ;
		int tag ;
	
		// Physics Fit Parameters 
		double gamma_in ;
		double dgam ;
		double Rt ;
		double Rp ;
		double delta1 ;
		double delta2 ;
		double delta_ms ;
		double phi_s ;
	
		// Other experimental parameters
		double tagFraction ;
		double resolution ;
	
		// Othere things calculated later on the fly
		double tlo, thi ;

		// Caching 
		bool normalisationCacheValid ;
		double normalisationCacheValue[3] ;


		//------ This is all stuff from Petes J/PsiPhi PDF ---------------------------
		    
		//Amplitudes Used in one angle PDF
		double AoAo() const ;   
		double AeAe() const ;
	
		//Amplitudes Used in three angle PDF
		double AT() const ;
		double AP() const ;
		double A0() const ;
	
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
	
		// Functions to help convolve a single gaussian into time primitives
		// DIDNT APPEAR TO BE USED RooComplex evalCerf( double, double, double ) const ;
		RooComplex evalCerfApprox( double, double, double ) const ;
		double evalCerfRe( double, double, double ) const ;
		double evalCerfIm( double, double, double ) const ;
	
		//---------------------
		// Some time primitive integrals
	
		// Integral of exp( - G * t ) from t1 to t2
		double intExp( double G, double t1, double t2 ) const ;
	
		// Integral of exp( - G * t ) * cos( dm * t )  from t1 to t2
		double intExpCos( double G, double dm, double t1, double t2 ) const ;
	
		// Integral of exp( - G * t ) * sin( dm * t )  from t1 to t2
		double intExpSin( double G, double dm, double t1, double t2 ) const ;
	
		//--------------------
		// Tag category, i.e B, Bbar or untagged.
		double q() const ;
	
	
		//------------------------------------------------------------------------------
		// These are the time factors and their analytic integrals for the one angle PDF
	
		//..................................
		double timeFactorEven(  )  const ;
		double timeFactorEvenInt(  )  const ;
	
		//..................................
		double timeFactorOdd(  )   const ;
		double timeFactorOddInt(  )  const ;
	
		//----------------------------------------------------------
		// These are the time factors and their analytic integrals for the three angle PDF
	
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
		double timeFactorReA0AP( )  const ;  
		double timeFactorReA0APInt( ) const ;
	
		//...........................
		double timeFactorImAPAT( ) const ; 
		double timeFactorImAPATInt( ) const ;
	
		//...........................
		double timeFactorImA0AT(  ) const ;
		double timeFactorImA0ATInt( ) const ;
    
	
		//------------------------------------------------------
		// Angle factors for one angle distributions
	
		double angleFactorEven(  )  const ;
		double angleFactorOdd(  )   const ;
	
		//------------------------------------------------------
		// Angle factors for three angle distributions
	
		double angleFactorA0A0(  ) const ;
		double angleFactorAPAP(  ) const ;
		double angleFactorATAT(  ) const ;
		double angleFactorReA0AP( ) const ;
		double angleFactorImAPAT(  ) const ;
		double angleFactorImA0AT(  ) const ;
	
		//-------------------------------------------------------
		// Putting it together for the differential cross section and its integrals.
	
		//...................................
		double diffXsec(  )  const ;    // 3 angles
		double diffXsecOne(  ) const ; // 1 angle
		
		//...................................
		// Integral over all variables: t + angles
		double diffXsecNorm1(  ) const ;
		double diffXsecOneNorm1(  ) const ;
	
		//...................................
		// Integral over angles only
		double diffXsecNorm2(  ) const ;
		double diffXsecOneNorm2(  ) const ;
	
	

};

#endif
