/**
        @class RooBs2PhiPhiFullPdf

        A RooFit format PDF for Bs decay to Phi Phi

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOBS2PHIPHIFULLPDF
#define ROOBS2PHIPHIFULLPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooCategoryProxy.h"

#include "RooComplex.h"
#include "RooMath.h"
#include "RooRandom.h"//added include
#include "TRandom3.h"
#include "TRandom.h"

class RooBs2PhiPhiFullPdf : public RooAbsPdf {
public:
  RooBs2PhiPhiFullPdf(const char *name, const char *title,
                      RooAbsReal& _theta1,
                      RooAbsReal& _theta2,
                      RooAbsReal& _kai,
                      RooAbsReal& _time,
                      //RooRealVar& _tag,
                      RooAbsCategory& _tag,
                      //Int_t tag,
                      RooAbsReal& _Rt,
                      RooAbsReal& _Rp,
                      RooAbsReal& _sp1,
                      RooAbsReal& _sp2,
                      RooAbsReal& _Gamma_L,
                      RooAbsReal& _Gamma_H,
                      RooAbsReal& _dms,
                      RooAbsReal& _phi_s,
                      RooAbsReal& _efftag,
                      RooAbsReal& _defftag,
                      RooAbsReal& _wtag,
                      RooAbsReal& _dwtag,
                      RooAbsReal& _aprod,
                      RooAbsReal& _sigtime,
                      RooAbsReal& _sigtime2,
                      RooAbsReal& _sigtheta,
                      RooAbsReal& _sigphi,
                      RooAbsReal& _biastime,
                      RooAbsReal& _accpara1,
                      RooAbsReal& _accpara2);
  RooBs2PhiPhiFullPdf(const RooBs2PhiPhiFullPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooBs2PhiPhiFullPdf(*this,newname); }
  inline virtual ~RooBs2PhiPhiFullPdf() { }

protected:

  RooRealProxy theta1 ;
  RooRealProxy theta2 ;
  RooRealProxy kai ;
  RooRealProxy time ;
  RooCategoryProxy tag ;
  //RooRealProxy tag ;
  //Int_t tag;
  RooRealProxy Rt ;
  RooRealProxy Rp ;
  RooRealProxy sp1 ;
  RooRealProxy sp2 ;
  RooRealProxy Gamma_L ;
  RooRealProxy Gamma_H ;
  RooRealProxy dms ;
  RooRealProxy phi_s ;
  RooRealProxy efftag ;
  RooRealProxy defftag ;
  RooRealProxy wtag ;
  RooRealProxy dwtag ;
  RooRealProxy aprod ;
  RooRealProxy sigtime ;
  RooRealProxy sigtime2 ;
  RooRealProxy sigtheta;
  RooRealProxy sigphi;
  RooRealProxy biastime ;
  RooRealProxy accpara1 ;
  RooRealProxy accpara2 ;

  Double_t evaluate() const ;

private:

  Double_t A1() const;
  Double_t A2() const;
  Double_t A3() const;
  Double_t A4() const;
  Double_t A5() const;
  Double_t A6() const;

  Double_t ang1() const;
  Double_t ang2() const;
  Double_t ang3() const;
  Double_t ang4() const;
  Double_t ang5() const;
  Double_t ang6() const;

  Double_t ang1_Inte() const;
  Double_t ang2_Inte() const;
  Double_t ang3_Inte() const;
  Double_t ang4_Inte() const;
  Double_t ang5_Inte() const;
  Double_t ang6_Inte() const;


  Double_t t1() const;
  Double_t t2() const;
  Double_t t3() const;
  Double_t t4() const;
  Double_t t5() const;
  Double_t t6() const;
  Double_t eff_t() const;

  Double_t wtageff() const;
  Double_t tageff() const;

  Double_t frac_flav() const;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,const char*) const ;
  Double_t analyticalIntegral(Int_t code, const char*) const ;

  //copy from RooGaussModel.rdl

  RooComplex evalCerfApprox(Double_t swt, Double_t u, Double_t c) const ;

  // Calculate exp(-u^2) cwerf(swt*c + i(u+c)), taking care of numerical instabilities
  inline RooComplex evalCerf(Double_t swt, Double_t u, Double_t c) const {
    RooComplex z(swt*c,u+c);
    return (z.im()>-4.0) ? RooMath::FastComplexErrFunc(z)*exp(-u*u) : evalCerfApprox(swt,u,c) ;
  }
    
  // Calculate Re(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
  inline Double_t evalCerfRe(Double_t swt, Double_t u, Double_t c) const {
    RooComplex z(swt*c,u+c);
    return (z.im()>-4.0) ? RooMath::FastComplexErrFuncRe(z)*exp(-u*u) : evalCerfApprox(swt,u,c).re() ;
  }
  
  // Calculate Im(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
  inline Double_t evalCerfIm(Double_t swt, Double_t u, Double_t c) const {
    RooComplex z(swt*c,u+c);
    return (z.im()>-4.0) ? RooMath::FastComplexErrFuncIm(z)*exp(-u*u) : evalCerfApprox(swt,u,c).im() ;
  }

  //contribution from exp(-Gamma_H*t) 
  Double_t ExpL() const;
  //contribution from exp(-Gamma_H*t)
  Double_t ExpH() const;
  //contribution from  exp(-(Gamma_H+Gamma_L)/2*t)*sin(dms*t)
  Double_t ExpSin() const;
  //contribution from  exp(-(Gamma_H+Gamma_L)/2*t)*cos(dms*t)
  Double_t ExpCos() const;
  
  Double_t reso_t(Double_t dt) const;

  Double_t reso_theta(Double_t dt, Double_t res) const;

  Double_t reso_phi(Double_t dt, Double_t res) const;

  Double_t CosTheta1() const;
  
  Double_t SinTheta1() const;
  
  Double_t CosTheta2() const;
  
  Double_t SinTheta2() const;
  
  Double_t Sin2Theta1() const;
  
  Double_t Sin2Theta2() const;
  
  Double_t SinPhi() const;
  
  Double_t CosPhi() const;
  
  Double_t Sin2Phi() const;
  
  Double_t  Cos2Phi() const;

  ClassDef(RooBs2PhiPhiFullPdf,0) // Your description goes here...
};


#endif
