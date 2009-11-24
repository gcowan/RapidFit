/**
        @class RooPdf_Bs2JPsiPhi

        A RooFit format PDF describing Bs decay to J/Psi Phi

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//-----------------------------------------------------------------------------
// Header file for class : RooPdf_Bs2JPsiPhi
//
// 2007-11-06 : Peter Clarke
//-----------------------------------------------------------------------------


/* Documentation

This is a PDF for use in fitting Bs -> Vector Vector Decays. 

It has been written for the Bs -> J/Psi Phi decay to analyse the full three decay angle
distribution in order to separate the different CP eigenstate components.

It also implements a mode for analysing just the single transversity angle (through a different constructor)

The PDF can handle integrals in the following ranges
- time:  arbitrary t1 to t2  in range [0,inf]
- all angles:  full range, no smaller limits allowed

The reference for the three angle distributions is:

Dighe, Dunietz, Fleisher /Eur. Phys. Journal C6 , 647-662 (1999)

The definition of phi_s is such that in the SM it is negative, i.e phi_s = - 2 Chi

The PDF is instantiation is:

  RooPdf_Bs2JPsiPhi signal (

      These are fixed parametes used only as construction 
       "name of pdf", 
       "description",
       btype                // Sets whether PDF is for B (=1),  Bbar (=1) or Untagged (=0)

      These are RooRealVar parameters 
       time,                // proper time for event decay
       ctheta_tr,           // transversity angle for event decay
       phi_tr,              // phi angle for event decay 
       ctheta_1,            // other angle for event decay 
       gamma,               // average Bs lifetime
       dgam,                // Bs lifetime difference 
       Rt,                  // Fraction of "transverse" or CP -1 final state
       Rp,                  // Fraction of "parallel" CP final state
       delta1,              // First strong phase  (see reference)
       delta2,              // Second strong phase
       delta_ms,            // Bs mass difference
       phi_s,               // CP violating phase
       tagFraction,         // Mistag fraction in range [0, 0.5]
       resolution,          // Detector resolution
    ) ;


The relation between the Rt and Rp input parameters and the amplitudes is:

 AT^2 = Rt
 AP^2 = Rp
 AO^2 = 1-Rt-Rp

with AT^2 + AP^2 + AO^2 = 0


*/

#ifndef ROOPDFBS2JPSIPHI
#define ROOPDFBS2JPSIPHI

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooMath.h"
#include "RooComplex.h"



class RooPdf_Bs2JPsiPhi : public RooAbsPdf 
{

protected:
 
  // Measured Event Attributes
  RooRealProxy t ;
  RooRealProxy ctheta_tr ;
  RooRealProxy phi_tr ;
  RooRealProxy ctheta_1 ;

  // Physics Fit Parameters 
  RooRealProxy gamma_in ;
  RooRealProxy dgam ;
  RooRealProxy Rt ;
  RooRealProxy Rp ;
  RooRealProxy delta1 ;
  RooRealProxy delta2 ;
  RooRealProxy delta_ms ;
  RooRealProxy phi_s ;

  // Other experimental parameters
  RooRealProxy tagFraction ;
  RooRealProxy resolution ;
  
  //.......................................
  // Evaluate mandatory method		     
  Double_t evaluate( )  const ;

 public:

  //..................................
  //Constructor for three angle PDF
  RooPdf_Bs2JPsiPhi(const char *name, const char *title, Int_t _btype, 
		    RooAbsReal& _t,
		    RooAbsReal& _ctheta_tr,
		    RooAbsReal& _phi_tr,
		    RooAbsReal& _ctheta_1,
		    RooAbsReal& _gamma,
		    RooAbsReal& _dgam,
		    RooAbsReal& _Rt,
		    RooAbsReal& _Rp,
		    RooAbsReal& _delta1,
		    RooAbsReal& _delta2,
		    RooAbsReal& _delta_ms,
		    RooAbsReal& _phi_s,
		    RooAbsReal& _tagFraction,
		    RooAbsReal& _resolution
		    ) ;

  //..................................
  //Constructor for one angle PDF
  RooPdf_Bs2JPsiPhi(const char *name, const char *title, Int_t _btype, 
		    RooAbsReal& _t,
		    RooAbsReal& _ctheta_tr,
		    RooAbsReal& _gamma,
		    RooAbsReal& _dgam,
		    RooAbsReal& _Rt,
		    RooAbsReal& _delta_ms,
		    RooAbsReal& _phi_s,
		    RooAbsReal& _tagFraction,
		    RooAbsReal& _resolution
		    ) ;

  //....................................
  //Copy constructor  
  RooPdf_Bs2JPsiPhi(const RooPdf_Bs2JPsiPhi& other, const char* name=0) ;

  //.....................................
  // Mandatory method    
  //virtual TObject* RooPdf_Bs2JPsiPhi::clone(const char* newname) const { 
  //  return new RooPdf_Bs2JPsiPhi(*this,newname); 
  //};
  //Stupid GCC upgrade wants the class specifier removed. Don't know why.
  virtual TObject* clone(const char* newname) const { 
    return new RooPdf_Bs2JPsiPhi(*this,newname); 
  };

  //......................................
  //Destructor
  inline virtual ~RooPdf_Bs2JPsiPhi() {};


private:

  // State variables
  const Int_t BTYPE ;  // Whether constructed for  B (+1), Bbar (-1) or  Untagged (0)
  const Int_t MODE  ;  // Whether constructed for one angle (-1) or for three angle (+1)
    
  //Amplitudes Used in one angle PDF
  Double_t AoAo() const ;   
  Double_t AeAe() const ;

  //Amplitudes Used in three angle PDF
  Double_t AT() const ;
  Double_t AP() const ;
  Double_t A0() const ;

  Double_t ctrsq() const ;
  Double_t strsq() const ;
  Double_t ct1sq() const ;
  Double_t st1sq() const ;
  Double_t cphsq() const ;
  Double_t sphsq() const ;

  // Widths
  Double_t gamma_l() const ;
  Double_t gamma_h() const ;
  Double_t gamma() const ;

  // Time primitives
  Double_t expL() const ;
  Double_t expH() const ;
  Double_t expCos() const ;
  Double_t expSin() const ;

  // Functions to help convolve a single gaussian into time primitives
  // DIDNT APPEAR TO BE USED RooComplex evalCerf( Double_t, Double_t, Double_t ) const ;
  RooComplex evalCerfApprox( Double_t, Double_t, Double_t ) const ;
  Double_t evalCerfRe( Double_t, Double_t, Double_t ) const ;
  Double_t evalCerfIm( Double_t, Double_t, Double_t ) const ;
  
  //---------------------
  // Some time primitive integrals

  // Integral of exp( - G * t ) from t1 to t2
  Double_t intExp( Double_t G, Double_t t1, Double_t t2 ) const ;

  // Integral of exp( - G * t ) * cos( dm * t )  from t1 to t2
  Double_t intExpCos( Double_t G, Double_t dm, Double_t t1, Double_t t2 ) const ;

  // Integral of exp( - G * t ) * sin( dm * t )  from t1 to t2
  Double_t intExpSin( Double_t G, Double_t dm, Double_t t1, Double_t t2 ) const ;

  //--------------------
  // Tag category, i.e B, Bbar or untagged.
  Double_t q() const ;


  //------------------------------------------------------------------------------
  // These are the time factors and their analytic integrals for the one angle PDF

  //..................................
  Double_t timeFactorEven(  )  const ;
  Double_t timeFactorEvenInt(  )  const ;

  //..................................
  Double_t timeFactorOdd(  )   const ;
  Double_t timeFactorOddInt(  )  const ;

  //----------------------------------------------------------
  // These are the time factors and their analytic integrals for the three angle PDF

  //...........................
  Double_t timeFactorA0A0( ) const ;      
  Double_t timeFactorA0A0Int( ) const ;
      
  //...........................
  Double_t timeFactorAPAP( ) const ;
  Double_t timeFactorAPAPInt( ) const ;
      
  //...........................
  Double_t timeFactorATAT( ) const ;
  Double_t timeFactorATATInt( ) const ;
     
 //...........................
  Double_t timeFactorReA0AP( )  const ;  
  Double_t timeFactorReA0APInt( ) const ;
     
  //...........................
  Double_t timeFactorImAPAT( ) const ; 
  Double_t timeFactorImAPATInt( ) const ;
 
  //...........................
  Double_t timeFactorImA0AT(  ) const ;
  Double_t timeFactorImA0ATInt( ) const ;
    

  //------------------------------------------------------
  // Angle factors for one angle distributions
	
  Double_t angleFactorEven(  )  const ;
  Double_t angleFactorOdd(  )   const ;
	
  //------------------------------------------------------
  // Angle factors for three angle distributions

  Double_t angleFactorA0A0(  ) const ;
  Double_t angleFactorAPAP(  ) const ;
  Double_t angleFactorATAT(  ) const ;
  Double_t angleFactorReA0AP( ) const ;
  Double_t angleFactorImAPAT(  ) const ;
  Double_t angleFactorImA0AT(  ) const ;

  //-------------------------------------------------------
  // Putting it together for the differential cross section and its integrals.

  //...................................
  Double_t diffXsec(  )  const ;    // 3 angles
  Double_t diffXsecOne(  ) const ; // 1 angle

  //...................................
  // Integral over all variables: t + angles
  Double_t diffXsecNorm1(  ) const ;
  Double_t diffXsecOneNorm1(  ) const ;
 
 //...................................
  // Integral over angles only
  Double_t diffXsecNorm2(  ) const ;
  Double_t diffXsecOneNorm2(  ) const ;

 //........................................
 // Optional method to advertise analytic integral
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* ) const ;

  //......................................
 // Optional method to calculate analytic integral
  Double_t analyticalIntegral(Int_t code, const char* ) const ;


  ClassDef(RooPdf_Bs2JPsiPhi, 0 ) // B2-> JPsi Phi three angle PDF


};


#endif

