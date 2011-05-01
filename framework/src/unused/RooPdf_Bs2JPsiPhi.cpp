// $Id: RooPdf_Bs2JPsiPhi.cpp,v 1.1 2009/11/10 10:35:48 gcowan Exp $
/**
        @class RooPdf_Bs2JPsiPhi

        A RooFit format PDF describing Bs decay to J/Psi Phi

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/
//-----------------------------------------------------------------------------
// Implementation file for class : RooPdf_Bs2JPsiPhi
//
// 2007-11-05 : Peter Clarke
//-----------------------------------------------------------------------------

#include "RooPdf_Bs2JPsiPhi.h"  

#include <TMath.h>
#include <cmath>
#include <iostream>

using namespace std;


//..................................
// Constructor for three angle
RooPdf_Bs2JPsiPhi::RooPdf_Bs2JPsiPhi(const char *name, const char *title, Int_t btype, 
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
    ) : 
    RooAbsPdf( name, title ), 
    t( "t", "time", this, _t ),
    ctheta_tr( "ctheta_tr", "ctheta_tr", this, _ctheta_tr ),
    phi_tr( "phi_tr", "phi_tr", this, _phi_tr ),
    ctheta_1( "ctheta_1", "ctheta_1", this, _ctheta_1 ),
    gamma_in("gamma", "gamma", this, _gamma ), 
    dgam("dgam", "Delta Gamma", this, _dgam ), 
    Rt("Rt", "Rt", this, _Rt ), 
    Rp("Rp", "Rp", this, _Rp ), 
    delta1( "delta1", "delta1", this, _delta1 ),
    delta2( "delta2", "delta2", this, _delta2 ),
    delta_ms( "delta_ms", "delta_ms", this, _delta_ms ),
    phi_s( "phi_s", "phi_s", this, _phi_s ),
    tagFraction( "tagFraction", "tagFraction", this, _tagFraction ),
    resolution( "resolution", "resolution", this, _resolution ),
    BTYPE(btype),
    MODE(1)
    {  }

  //..................................
  // Constructor for one angle
RooPdf_Bs2JPsiPhi::RooPdf_Bs2JPsiPhi(const char *name, const char *title, Int_t btype,
		    RooAbsReal& _t,
		    RooAbsReal& _ctheta_tr,
		    RooAbsReal& _gamma,
		    RooAbsReal& _dgam,
		    RooAbsReal& _Rt,
		    RooAbsReal& _delta_ms,
		    RooAbsReal& _phi_s,
		    RooAbsReal& _tagFraction,
		    RooAbsReal& _resolution
   ) : 
    RooAbsPdf( name, title ), 
    t( "t", "time", this, _t ),
    ctheta_tr( "ctheta_tr", "ctheta_tr", this, _ctheta_tr ),
    gamma_in("gamma", "gamma", this, _gamma ), 
    dgam("dgam", "Delta Gamma", this, _dgam ), 
    Rt("Rt", "Rt", this, _Rt ), 
    delta_ms( "delta_ms", "delta_ms", this, _delta_ms ),
    phi_s( "phi_s", "phi_s", this, _phi_s ),
    tagFraction( "tagFraction", "tagFraction", this, _tagFraction ),
    resolution( "resolution", "resolution", this, _resolution ),
    BTYPE(btype),
    MODE(-1)
    { }


  //....................................
  //Copy constructor  
  RooPdf_Bs2JPsiPhi::RooPdf_Bs2JPsiPhi(const RooPdf_Bs2JPsiPhi& other, const char* name) :
    RooAbsPdf( other, name ), 
    t( "t", this, other.t ),
    ctheta_tr( "ctheta_tr", this, other.ctheta_tr ),
    phi_tr( "phi_tr", this, other.phi_tr ),
    ctheta_1( "ctheta_1", this, other.ctheta_1 ),
    gamma_in( "gamma", this, other.gamma_in), 
    dgam( "dgam", this, other.dgam), 
    Rt( "Rt", this, other.Rt), 
    Rp( "Rp", this, other.Rp), 
    delta1( "delta1", this, other.delta1 ),
    delta2( "delta2", this, other.delta2 ),
    delta_ms( "delta_ms", this, other.delta_ms ),
    phi_s( "phi_s", this, other.phi_s),
    tagFraction( "tagFraction", this, other.tagFraction),
    resolution( "resolution", this, other.resolution),
    BTYPE(other.BTYPE),
    MODE(other.MODE)
    { } 

 
  //....................................
  //Internal helper functions
  
  //Amplitudes Used in one angle PDF
  Double_t RooPdf_Bs2JPsiPhi::AoAo() const  { return Rt; };   
  Double_t RooPdf_Bs2JPsiPhi::AeAe() const { return 1-Rt ; };

  //Amplitudes Used in three angle PDF
  Double_t RooPdf_Bs2JPsiPhi::AT() const { return sqrt(Rt) ; };
  Double_t RooPdf_Bs2JPsiPhi::AP() const { return sqrt(Rp) ; };
  Double_t RooPdf_Bs2JPsiPhi::A0() const { if( (1-Rt-Rp) < 0 ) return 0; else return sqrt(1-Rt-Rp) ; };

  Double_t RooPdf_Bs2JPsiPhi::ctrsq() const { return (ctheta_tr*ctheta_tr) ; }
  Double_t RooPdf_Bs2JPsiPhi::strsq() const { return (1.0 - ctheta_tr*ctheta_tr) ; }
  Double_t RooPdf_Bs2JPsiPhi::ct1sq() const { return (ctheta_1*ctheta_1) ; }
  Double_t RooPdf_Bs2JPsiPhi::st1sq() const { return (1.0 - ctheta_1*ctheta_1) ; }
  Double_t RooPdf_Bs2JPsiPhi::cphsq() const { return (cos(phi_tr)*cos(phi_tr)) ; }
  Double_t RooPdf_Bs2JPsiPhi::sphsq() const { return (sin(phi_tr)*sin(phi_tr)) ; }

  Double_t RooPdf_Bs2JPsiPhi::gamma_l() const { return gamma() + ( dgam / 2.0 ) ; }
  Double_t RooPdf_Bs2JPsiPhi::gamma_h() const { return gamma() - ( dgam / 2.0 ) ; }
  Double_t RooPdf_Bs2JPsiPhi::gamma() const   { return gamma_in ; }

  Double_t RooPdf_Bs2JPsiPhi::q() const { return BTYPE ;}


  //-----------------------------------
  // Time primitives including single gaussian resolution
  //
  // Much of the single gausian resolutino code is copied from Yue Hongs pdf
  
  Double_t RooPdf_Bs2JPsiPhi::expL() const 
  {
    if(resolution>0.) {     

      Double_t theExp = exp( -t*gamma_l() + resolution*resolution * gamma_l()*gamma_l() / 2. ) ;
      Double_t theErfc = RooMath::erfc(  -( t - resolution*resolution*gamma_l() ) /sqrt(2.)/resolution )  ;
      return theExp * theErfc  / 2.0 ;

      // Yue hongs code
      //Double_t c = gamma_l() * resolution /sqrt(2.); 
      //Double_t u = t / resolution / sqrt(2.);
      //return exp( c*c - gamma_l()*t ) * RooMath::erfc(c-u) / 2.;
    }	
    else 
      if( t < 0.0 ) return 0.0 ;
      return exp( -gamma_l() * t ) ;
  }

  Double_t RooPdf_Bs2JPsiPhi::expH() const 
  {
    if(resolution>0.) {     

      Double_t theExp = exp( -t*gamma_h() + resolution*resolution * gamma_h()*gamma_h() / 2. ) ;
      Double_t theErfc = RooMath::erfc(  -( t - resolution*resolution*gamma_h() ) /sqrt(2.)/resolution )  ;
      return theExp * theErfc  / 2.0 ;

      //Yue Hongs code
      //Double_t c = gamma_h() * resolution /sqrt(2.); 
      //Double_t u = t / resolution / sqrt(2.);
      //return exp( c*c - gamma_h()*t ) * RooMath::erfc(c-u) / 2.;
    }	
    else 
      if( t < 0.0 ) return 0.0 ;
      return exp( -gamma_h() * t ) ;
  }

  Double_t RooPdf_Bs2JPsiPhi::expSin() const  
  {
    if( resolution > 0. ) {

      //Double_t theExp = exp( -t*gamma() + resolution*resolution * ( gamma()*gamma() - delta_ms*delta_ms ) / 2. ) ;
      //Double_t theCos = cos( delta_ms * ( t - resolution*resolution*gamma() ) ) ;
      //Double_t theSin = sin( delta_ms * ( t - resolution*resolution*gamma() ) ) ;
      //RooComplex z( -( t - resolution*resolution*gamma() )/sqrt(2.)/resolution,  - resolution*delta_ms/sqrt(2.) ) ;
      //Double_t theReErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncRe(z) ) : ( 1.0 - RooMath::FastComplexErrFuncRe(z) );
      //Double_t theImErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncIm(z) ) : ( 1.0 - RooMath::FastComplexErrFuncIm(z) );
      //return theExp * ( theCos*theImErfc + theSin*theReErfc ) / 2.0 ;

      // Yue Hongs code
      Double_t c = gamma() * resolution/sqrt(2.);
      Double_t u = t / resolution / sqrt(2.) ;
      Double_t wt = delta_ms / gamma() ;
      return ( evalCerfIm(wt,-u,c) - evalCerfIm(-wt,-u,c) ) /4. ;
	}
	else
      if( t < 0.0 ) return 0.0 ;
      return exp( -gamma() *t ) * sin( delta_ms * t )  ;

  }

  Double_t RooPdf_Bs2JPsiPhi::expCos() const 
  {
    if( resolution > 0. ) {
 
      //Double_t theExp = exp( -t*gamma() + resolution*resolution * ( gamma()*gamma() - delta_ms*delta_ms ) / 2. ) ;
      //Double_t theCos = cos( delta_ms * ( t - resolution*resolution*gamma() ) ) ;
      //Double_t theSin = sin( delta_ms * ( t - resolution*resolution*gamma() ) ) ;
      //RooComplex z( -( t - resolution*resolution*gamma() )/sqrt(2.)/resolution,  - resolution*delta_ms/sqrt(2.) ) ;
      //Double_t theReErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncRe(z) ) : ( 1.0 - RooMath::FastComplexErrFuncRe(z) );
      //Double_t theImErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncIm(z) ) : ( 1.0 - RooMath::FastComplexErrFuncIm(z) );
      //return theExp * ( theCos*theReErfc - theSin*theImErfc ) / 2.0 ;

	  //Yue Hongs code
      Double_t c = gamma() * resolution/sqrt(2.);
      Double_t u = t / resolution / sqrt(2.) ;
      Double_t wt = delta_ms / gamma() ;
      return ( evalCerfRe(wt,-u,c) + evalCerfRe(-wt,-u,c) ) / 4. ;
	}
	else
      if( t < 0.0 ) return 0.0 ;
      return exp( -gamma() *t ) * cos( delta_ms * t )  ;

  }

  // All of these functions were taken from Yue Hongs code
  
  // Calculate exp(-u^2) cwerf(swt*c + i(u+c)), taking care of numerical instabilities
  //RooComplex RooPdf_Bs2JPsiPhi::evalCerf(Double_t swt, Double_t u, Double_t c) const {
  //  RooComplex z(swt*c,u+c);
  //  return (z.im()>-4.0) ? RooMath::FastComplexErrFunc(z)*exp(-u*u) : evalCerfApprox(swt,u,c) ;
  ///}  DIDNT APPEAR TO BE USED
    
  // Calculate Re(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
  Double_t RooPdf_Bs2JPsiPhi::evalCerfRe(Double_t swt, Double_t u, Double_t c) const {
    RooComplex z(swt*c,u+c);
    return (z.im()>-4.0) ? RooMath::FastComplexErrFuncRe(z)*exp(-u*u) : evalCerfApprox(swt,u,c).re() ;
  }
  
  // Calculate Im(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
  Double_t RooPdf_Bs2JPsiPhi::evalCerfIm(Double_t swt, Double_t u, Double_t c) const {
    RooComplex z(swt*c,u+c);
    return (z.im()>-4.0) ? RooMath::FastComplexErrFuncIm(z)*exp(-u*u) : evalCerfApprox(swt,u,c).im() ;
  }

  // use the approximation: erf(z) = exp(-z*z)/(sqrt(pi)*z) to explicitly cancel the divergent exp(y*y) behaviour of
  // CWERF for z = x + i y with large negative y
  RooComplex RooPdf_Bs2JPsiPhi::evalCerfApprox(Double_t swt, Double_t u, Double_t c) const
  {
     static Double_t rootpi= sqrt(atan2(0.,-1.));
     RooComplex z(swt*c,u+c);  
     RooComplex zc(u+c,-swt*c);
     RooComplex zsq= z*z;
     RooComplex v= -zsq - u*u;
     return v.exp()*(-zsq.exp()/(zc*rootpi) + 1)*2 ; //why shoule be a 2 here?
     //   return v.exp()*(-zsq.exp()/(zc*rootpi) + 1);
 }


  //---------------------
  // Some time integrals

  // Integral of exp( - G * t ) from t1 to t2
  Double_t RooPdf_Bs2JPsiPhi::intExp( Double_t G, Double_t t1, Double_t t2 ) const {
    return (1/G) * (exp(-G*t1) - exp(-G*t2) ) ;
  }

  // Integral of exp( - G * t ) * cos( dm * t )  from t1 to t2
  Double_t RooPdf_Bs2JPsiPhi::intExpCos( Double_t G, Double_t dm, Double_t t1, Double_t t2 ) const {
    return (1/(G*G + dm*dm)) * (
				( exp(-G*t1)* (G*cos(dm*t1) - dm*sin(dm*t1)))
			       -( exp(-G*t2)* (G*cos(dm*t2) - dm*sin(dm*t2)))
				);
  }

  // Integral of exp( - G * t ) * sin( dm * t )  from t1 to t2
  Double_t RooPdf_Bs2JPsiPhi::intExpSin( Double_t G, Double_t dm, Double_t t1, Double_t t2 ) const {
    return (1/(G*G + dm*dm)) * (
				( exp(-G*t1)* (G*sin(dm*t1) + dm*cos(dm*t1)))
			       -( exp(-G*t2)* (G*sin(dm*t2) + dm*cos(dm*t2)))
				);
  }

//------------------------------------------------------------------------------
// These are the time factors and their analytic integrals for the one angle PDF

  //..................................
  Double_t RooPdf_Bs2JPsiPhi::timeFactorEven(  )  const
    {
      //if( t < 0.0 ) return 0.0 ;
      Double_t result = 
	  ( 1.0 + cos(phi_s) ) * expL( ) 
        + ( 1.0 - cos(phi_s) ) * expH( ) 
        + q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) ;
      return result ;
    };

  Double_t RooPdf_Bs2JPsiPhi::timeFactorEvenInt(  )  const
    {
      Double_t tlo = t.min() ;
      Double_t thi = t.max() ;
      if(tlo < 0.) tlo = 0. ;

      Double_t result = 
	  ( 1.0 + cos(phi_s) )  * intExp( gamma_l(), tlo, thi )     
        + ( 1.0 - cos(phi_s) )  * intExp( gamma_h(), tlo, thi )          
	+ q() * ( 2.0 * sin(phi_s)   ) * intExpSin( gamma(), delta_ms, tlo, thi ) * (1.0 - 2.0*tagFraction) ;
      return result ;
    };


  //..................................
  Double_t RooPdf_Bs2JPsiPhi::timeFactorOdd(  )   const
    {
      //if( t < 0.0 ) return 0.0 ;
      Double_t result = 
          ( 1.0 - cos(phi_s) ) * expL( ) 
        + ( 1.0 + cos(phi_s) ) * expH( ) 
        - q() * ( 2.0 * sin(phi_s)   ) * expSin( ) * (1.0 - 2.0*tagFraction) ;
      return result ;
    };

  Double_t RooPdf_Bs2JPsiPhi::timeFactorOddInt(  )  const
    {
      Double_t tlo = t.min() ;
      Double_t thi = t.max() ;
      if(tlo < 0.) tlo = 0. ;

      Double_t result = 
	  ( 1.0 - cos(phi_s) ) * intExp( gamma_l(), tlo, thi )
        + ( 1.0 + cos(phi_s) ) * intExp( gamma_h(), tlo, thi ) 
        - q() * ( 2.0 * sin(phi_s)   ) * intExpSin( gamma(), delta_ms, tlo, thi ) * (1.0 - 2.0*tagFraction) ;
      return result ;
    };


//----------------------------------------------------------
// These are the time factors and their analytic integrals for the three angle PDF

  //...........................
  Double_t RooPdf_Bs2JPsiPhi::timeFactorA0A0( )    const { return timeFactorEven( ) ; } ;      
  Double_t RooPdf_Bs2JPsiPhi::timeFactorA0A0Int( ) const { return timeFactorEvenInt( ) ; } ;
      
  //...........................
  Double_t RooPdf_Bs2JPsiPhi::timeFactorAPAP( )    const { return timeFactorEven( ) ; } ;
  Double_t RooPdf_Bs2JPsiPhi::timeFactorAPAPInt( ) const { return timeFactorEvenInt( ) ; } ;
      
  //...........................
  Double_t RooPdf_Bs2JPsiPhi::timeFactorATAT( )    const { return timeFactorOdd( ) ; } ;
  Double_t RooPdf_Bs2JPsiPhi::timeFactorATATInt( ) const { return timeFactorOddInt( ) ; } ;
     
 //...........................
  Double_t RooPdf_Bs2JPsiPhi::timeFactorReA0AP( )  const
    {
      //if( t < 0.0 ) return 0.0 ;
      Double_t result = cos(delta2-delta1) * this->timeFactorEven(  ) ;
      return result ;
    } ;
     
  Double_t RooPdf_Bs2JPsiPhi::timeFactorReA0APInt( ) const
    {
      Double_t result = cos(delta2-delta1) * this->timeFactorEvenInt( ) ;
      return result ;
    } ;
     
//...........................
  Double_t RooPdf_Bs2JPsiPhi::timeFactorImAPAT( ) const
    {
      //if( t < 0.0 ) return 0.0 ;
      Double_t result = 
	    q() * 2.0  * ( sin(delta1)*expCos( ) - cos(delta1)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)
	  - 1.0 * ( expH( ) - expL( ) ) * cos(delta1) * sin(phi_s)  ;
	    
	  return result ;
    } ;
     
  Double_t RooPdf_Bs2JPsiPhi::timeFactorImAPATInt( ) const
    {
      Double_t tlo = t.min() ;
      Double_t thi = t.max() ;
      if(tlo < 0.) tlo = 0. ;

      Double_t result = 
	q() * 2.0  * ( sin(delta1)*intExpCos(gamma(),delta_ms,tlo,thi) - cos(delta1)*cos(phi_s)*intExpSin(gamma(),delta_ms,tlo,thi) ) * (1.0 - 2.0*tagFraction)
	- 1.0 * ( intExp(gamma_h(),tlo,thi) - intExp(gamma_l(),tlo,thi) ) * cos(delta1) * sin(phi_s) ;	    
         return result ;
    } ;
     

//...........................
  Double_t RooPdf_Bs2JPsiPhi::timeFactorImA0AT(  ) const
    {
      //if( t < 0.0 ) return 0.0 ;
      Double_t result =
	    q() * 2.0  * ( sin(delta2)*expCos( ) - cos(delta2)*cos(phi_s)*expSin( ) ) * (1.0 - 2.0*tagFraction)	
	   -1.0 * ( expH( ) - expL( ) ) * cos(delta2) * sin(phi_s) ;
       return result ;
    } ;

  Double_t RooPdf_Bs2JPsiPhi::timeFactorImA0ATInt( ) const
    {
      Double_t tlo = t.min() ;
      Double_t thi = t.max() ;
      if(tlo < 0.) tlo = 0. ;

      Double_t result = 
	q() * 2.0  * ( sin(delta2)*intExpCos(gamma(),delta_ms,tlo,thi) - cos(delta2)*cos(phi_s)*intExpSin(gamma(),delta_ms,tlo,thi)  ) * (1.0 - 2.0*tagFraction)
	-1.0 * ( intExp(gamma_h(),tlo,thi) - intExp(gamma_l(),tlo,thi)  ) * cos(delta2) * sin(phi_s) ;
       return result ;
    } ;
     

//------------------------------------------------------
// Angle factors for one angle PDF
	
  //.................................
  Double_t RooPdf_Bs2JPsiPhi::angleFactorEven(  )  const
    {
      // Note that this is normalised to 1
      Double_t result = 3.0/8.0 * (1.0 + ctrsq() ) ;
      return result ;
    };

  //.................................
  Double_t RooPdf_Bs2JPsiPhi::angleFactorOdd(  )   const
    {
      // Note that this is normalised to 1
      Double_t result = 3.0/4.0 * (1.0 - ctrsq() ) ;
      return result ;
    };

	
//------------------------------------------------------
// Angle factors for three angle PDFs
	
	
  //...........................
  Double_t RooPdf_Bs2JPsiPhi::angleFactorA0A0(  ) const
	{
          // Normalised to  1	
	  Double_t result = 2.0 * ct1sq() * (1.0 - strsq()*cphsq() ) * (9.0/32.0/TMath::Pi());
	  return result ;	
	};

  //...........................
  Double_t RooPdf_Bs2JPsiPhi::angleFactorAPAP(  ) const
	{
          // Normalised to  1
	  Double_t result =  st1sq() * (1.0 - strsq()*sphsq() ) * (9.0/32.0/TMath::Pi());
	  return result ;	
	};

  //...........................
  Double_t RooPdf_Bs2JPsiPhi::angleFactorATAT(  ) const
	{
          // Normalised to  1
	  Double_t result = st1sq() * strsq() * (9.0/32.0/TMath::Pi());
	  return result ;
	
	};

  //...........................
  Double_t RooPdf_Bs2JPsiPhi::angleFactorReA0AP( ) const
	{
          // Normalised to  0
          Double_t theta_1 = acos(ctheta_1) ;	
	  Double_t result =    sin(2.0*theta_1) * strsq() * sin(2.0*phi_tr) / sqrt(2.0) * (9.0/32.0/TMath::Pi());
	  return result ;	
	};

  //...........................
  Double_t RooPdf_Bs2JPsiPhi::angleFactorImAPAT(  ) const
	{
          // Normalised to  0
          Double_t theta_tr = acos(ctheta_tr) ;		
	  Double_t result =   -1.0 *  st1sq() * sin(2.0*theta_tr) * sin(phi_tr) * (9.0/32.0/TMath::Pi()) ;
	  return result ;	
	};

  //...........................
  Double_t RooPdf_Bs2JPsiPhi::angleFactorImA0AT(  ) const
	{
          // Normalised to  0
          Double_t theta_tr = acos(ctheta_tr) ;		
          Double_t theta_1 = acos(ctheta_1) ;		
	  Double_t result =  +1.0*   sin(2.0*theta_1) * sin(2.0*theta_tr) * cos(phi_tr) / sqrt(2.0) * (9.0/32.0/TMath::Pi());
	  return result ;	
	};


//-------------------------------------------------------------
// Putting it all together to make up the differential cross sections.

  //.......................................
  // Evaluate mandatory method		     
  Double_t RooPdf_Bs2JPsiPhi::evaluate( )  const
   {
     if( MODE == 1 ) 
       return this->diffXsec( );
      else 
        return this->diffXsecOne( );
   };

  //...................................
  // Diff cross sections

  Double_t RooPdf_Bs2JPsiPhi::diffXsec(  )  const
    {   
      Double_t xsec = 
        0.5 * A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) +
        0.5 * AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) +
        0.5 * AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) +
        0.5 * A0()*AP() * timeFactorReA0AP(  ) * angleFactorReA0AP( ) +
        0.5 * AP()*AT() * timeFactorImAPAT(  ) * angleFactorImAPAT( ) +
        0.5 * A0()*AT() * timeFactorImA0AT(  ) * angleFactorImA0AT( ) ;
      
      return xsec ;
    };

  Double_t  RooPdf_Bs2JPsiPhi::diffXsecOne(  ) const
    {
      Double_t result = 
	0.5 * AeAe() * timeFactorEven(  ) * angleFactorEven(  )  +
	0.5 * AoAo() * timeFactorOdd(  )  * angleFactorOdd(  ) ;
      return result ;
    };

  //...................................
  // Integral over all variables: t + angles

  Double_t RooPdf_Bs2JPsiPhi::diffXsecNorm1(  ) const
    {      
     Double_t norm = 
       0.5 * A0()*A0() * timeFactorA0A0Int(  ) +    // Angle factors normalised to 1
       0.5 * AP()*AP() * timeFactorAPAPInt(  ) +
       0.5 * AT()*AT() * timeFactorATATInt(  ) ;

      return norm ;
    };

  Double_t RooPdf_Bs2JPsiPhi::diffXsecOneNorm1(  ) const
    {      
     Double_t norm = 
	0.5 * AeAe() * timeFactorEvenInt(  )  +    // Angle factors normalised to 1
	0.5 * AoAo() * timeFactorOddInt(  )   ;
      return norm ;
    };


  //...................................
  // Integral over angles only 3 

  Double_t RooPdf_Bs2JPsiPhi::diffXsecNorm2(  ) const
    {          
     Double_t norm = 
       0.5 * A0()*A0() * timeFactorA0A0(  ) +    // Angle factors normalised to 1
       0.5 * AP()*AP() * timeFactorAPAP(  ) +
       0.5 * AT()*AT() * timeFactorATAT(  ) ;

      return norm ;
    };

  Double_t RooPdf_Bs2JPsiPhi::diffXsecOneNorm2(  ) const
    {          
     Double_t norm = 
	0.5 * AeAe() * timeFactorEven(  )  +     // Angle factors normalised to 1
	0.5 * AoAo() * timeFactorOdd(  )   ;
      return norm ;
    };

 //........................................
 // Optional method to advertise analytic integral
  Int_t RooPdf_Bs2JPsiPhi::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* ) const
 {
   if( MODE == 1 ) {
     if (matchArgs(allVars, analVars, t, ctheta_tr, phi_tr, ctheta_1)) return 1;
     if (matchArgs(allVars, analVars, ctheta_tr, phi_tr, ctheta_1)) return 2;
     return 0;
        }
   else {
     if (matchArgs(allVars, analVars, t, ctheta_tr)) return 1;
     if (matchArgs(allVars, analVars, ctheta_tr)) return 2;
     return 0;
     }
 }

  //......................................
 // Optional method to calculate analytic integral
  Double_t RooPdf_Bs2JPsiPhi::analyticalIntegral(Int_t code, const char* ) const
 {
   if( MODE == 1 ) {
     if (code==1)  return this->diffXsecNorm1( ) ;
     if (code==2)  return this->diffXsecNorm2( ) ;
     return 0 ;
     }
   else {
     if (code==1)  return this->diffXsecOneNorm1( ) ;
     if (code==2)  return this->diffXsecOneNorm2( ) ;
     return 0 ;
     }
 }


