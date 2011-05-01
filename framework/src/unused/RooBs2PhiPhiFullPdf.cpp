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

 // -- CLASS DESCRIPTION [PDF] -- 
 // Your description goes here... 
//Analytic integrals removed...
#include <iostream> 

#include "RooBs2PhiPhiFullPdf.h" 
#include "RooAbsReal.h" 

#include "TF1.h"
#include "TRandom3.h"
#include "TRandom.h"

 ClassImp(RooBs2PhiPhiFullPdf) 

   //what to do about asymmetry in untagged events?


 RooBs2PhiPhiFullPdf::RooBs2PhiPhiFullPdf(const char *name, const char *title, 
                                          RooAbsReal& _theta1,
                                          RooAbsReal& _theta2,
                                          RooAbsReal& _kai,
                                          RooAbsReal& _time,
                                          //	RooRealVar& _tag,
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
                                          RooAbsReal& _defftag,//define defftag as [eff(Bbar) - eff(B)]/2
                                          RooAbsReal& _wtag,
                                          RooAbsReal& _dwtag,//define dwtag as [wtag(Bbar) - wtag(B)]/2
                                          RooAbsReal& _aprod,
                                          RooAbsReal& _sigtime,//add additional time resolutions 
                                          RooAbsReal& _sigtime2,
                                          RooAbsReal& _sigtheta,
                                          RooAbsReal& _sigphi,
                                          RooAbsReal& _biastime,//and generate a random number to choose
                                          RooAbsReal& _accpara1,//which gaussian is used - add for angles
                                          RooAbsReal& _accpara2) ://also
   RooAbsPdf(name,title), 
   theta1("theta1","theta1",this,_theta1),
   theta2("theta2","theta2",this,_theta2),
   kai("kai","kai",this,_kai),
   time("time","time",this,_time),
   tag("tag","tag",this,_tag),
   Rt("Rt","Rt",this,_Rt),
   Rp("Rp","Rp",this,_Rp),
   sp1("sp1","sp1",this,_sp1),
   sp2("sp2","sp2",this,_sp2),
   Gamma_L("Gamma_L","Gamma_L",this,_Gamma_L),
   Gamma_H("Gamma_H","Gamma_H",this,_Gamma_H),
   dms("dms","dms",this,_dms),
   phi_s("phi_s","phi_s",this,_phi_s),
   efftag("efftag","efftag",this,_efftag),
   defftag("defftag","defftag",this,_defftag),
   wtag("wtag","wtag",this,_wtag),
   dwtag("dwtag","dwtag",this,_dwtag),
   aprod("aprod","aprod",this,_aprod),
   sigtime("sigtime","sigtime",this,_sigtime),
   sigtime2("sigtime2","sigtime2",this,_sigtime2),
   sigtheta("sigtheta","sigtheta",this,_sigtheta),
   sigphi("sigphi", "sigphi",this, _sigphi),
   biastime("biastime","biastime",this,_biastime),
   accpara1("accpara1","accpara1",this,_accpara1),
   accpara2("accpara2","accpara2",this,_accpara2)
 { 
 } 


 RooBs2PhiPhiFullPdf::RooBs2PhiPhiFullPdf(const RooBs2PhiPhiFullPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   theta1("theta1",this,other.theta1),
   theta2("theta2",this,other.theta2),
   kai("kai",this,other.kai),
   time("time",this,other.time),
   tag("tag",this,other.tag),
   Rt("Rt",this,other.Rt),
   Rp("Rp",this,other.Rp),
   sp1("sp1",this,other.sp1),
   sp2("sp2",this,other.sp2),
   Gamma_L("Gamma_L",this,other.Gamma_L),
   Gamma_H("Gamma_H",this,other.Gamma_H),
   dms("dms",this,other.dms),
   phi_s("phi_s",this,other.phi_s),
   efftag("efftag",this,other.efftag),
   defftag("defftag",this,other.defftag),
   wtag("wtag",this,other.wtag),
   dwtag("dwtag",this,other.dwtag),
   aprod("aprod",this,other.aprod),
   sigtime("sigtime",this,other.sigtime),
    sigtime2("sigtime2",this,other.sigtime2),
   sigtheta("sigtheta",this,other.sigtheta),
   sigphi("sigphi",this,other.sigphi),
   biastime("biastime",this,other.biastime),
   accpara1("accpara1",this,other.accpara1),
   accpara2("accpara2",this,other.accpara2)
 { 
 } 


Double_t RooBs2PhiPhiFullPdf::wtageff() const
{
  Double_t v = tag;
  if(tag==0) v = 0.0;
  else v = wtag - (dwtag*tag);  
  return v;
}

Double_t RooBs2PhiPhiFullPdf::tageff() const
{
  Double_t v = efftag;
  if(tag==0) v = aprod - 2*defftag;
  else v = tag; 
  return v;
}


Double_t RooBs2PhiPhiFullPdf::evaluate() const 
{
  

   Double_t v = A1()*ang1()*t1()*eff_t()+
                A2()*ang2()*t2()*eff_t()+
                A3()*ang3()*t3()*eff_t()+
                A4()*ang4()*t4()*eff_t()+
                A5()*ang5()*t5()*eff_t()+
                A6()*ang6()*t6()*eff_t();
   v *= frac_flav();
   return v;
 } 

 Double_t RooBs2PhiPhiFullPdf::eff_t() const
 {
   if (time<0) return 0;
   else if ((accpara1>1)&&(accpara2>1)) 
     return  accpara1*time*time*time/(1 + accpara2*time*time*time);
   else return accpara1*time*time*time/(accpara2+time*time*time);

   //return 1.;
   //return (accpara1*time*time)/(1 + (accpara2*time)*(accpara2*time));
 }

 // |A0*A0|
 Double_t RooBs2PhiPhiFullPdf::A1() const
 {
//   return  sqrt(1-Rt-Rp);  //bug
   return  1-Rt-Rp;
 }

 // |ApAp|
 Double_t RooBs2PhiPhiFullPdf::A2() const
 {
   return  Rp;
 }

 // |AtAt|
 Double_t RooBs2PhiPhiFullPdf::A3() const
 {
   return Rt;
 }

 // |ApAt|
 Double_t RooBs2PhiPhiFullPdf::A4() const
 {
   return sqrt(Rp*Rt);
 }

 // |A0Ap|
 Double_t RooBs2PhiPhiFullPdf::A5() const
 {
   return sqrt(1-Rt-Rp)*sqrt(Rp);
 }

 // |A0At|
 Double_t RooBs2PhiPhiFullPdf::A6() const
 {
  return sqrt(1-Rt-Rp)*sqrt(Rt);
 }

 Double_t RooBs2PhiPhiFullPdf::ang1() const
 {
   Double_t v=4.*cos(theta1)*cos(theta1)*cos(theta2)*cos(theta2);
   return v*sin(theta1)*sin(theta2);
 }

 Double_t RooBs2PhiPhiFullPdf::ang2() const
 {
   Double_t v=sin(theta1)*sin(theta1)*sin(theta2)*sin(theta2)*
              (1+cos(2*kai));
   return v*sin(theta1)*sin(theta2);
 }

 Double_t RooBs2PhiPhiFullPdf::ang3() const
 {
   Double_t v=sin(theta1)*sin(theta1)*sin(theta2)*sin(theta2)*
              (1-cos(2*kai));
   return v*sin(theta1)*sin(theta2);
 }

 Double_t RooBs2PhiPhiFullPdf::ang4() const
 {
   Double_t v=-2*sin(theta1)*sin(theta1)*sin(theta2)*sin(theta2)*
              sin(2*kai);
   return v*sin(theta1)*sin(theta2);
 }

 Double_t RooBs2PhiPhiFullPdf::ang5() const
 {
   Double_t v=sqrt(2.)*sin(2*theta1)*sin(2*theta2)*cos(kai);
   return v*sin(theta1)*sin(theta2);
 }

 Double_t RooBs2PhiPhiFullPdf::ang6() const
 {
   Double_t v=-sqrt(2.)*sin(2*theta1)*sin(2*theta2)*sin(kai);
   return v*sin(theta1)*sin(theta2);
 }

 Double_t RooBs2PhiPhiFullPdf::ang1_Inte() const
 {
      return 3.14159*32./9.;
   //return 2*3.14159*64./9.;
}

 Double_t RooBs2PhiPhiFullPdf::ang2_Inte() const
 {
       return 3.14159*32./9.; 
   //return 2*3.14159*16./9.;
 }

 Double_t RooBs2PhiPhiFullPdf::ang3_Inte() const
 {
        return 3.14159*32./9.;
   //return 2*3.14159*16./9.;
 }

 Double_t RooBs2PhiPhiFullPdf::ang4_Inte() const
 {
   return 0.;
 }

 Double_t RooBs2PhiPhiFullPdf::ang5_Inte() const
 {
   return 0.;
 }

 Double_t RooBs2PhiPhiFullPdf::ang6_Inte() const
 {
   return 0.;
 }

 Int_t RooBs2PhiPhiFullPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,const char* ) const
 {
      if (matchArgs(allVars, analVars, theta1,theta2,kai)){
   return 1;
   }
      /*    
     if (matchArgs(allVars, analVars, theta2,kai)){
      cout<<"Using Analytic integral 2"<<endl;
      return 2;
      }
      if (matchArgs(allVars, analVars, theta1,kai)){
      cout<<"Using Analytic Integral 3"<<endl;
      return 3;
      }
      if (matchArgs(allVars, analVars, theta1,theta2)){
      cout<<"Using Analytic Integral 4"<<endl;
      return 4; 
      }
      */
      else return 0;
 }

 Double_t RooBs2PhiPhiFullPdf::analyticalIntegral(Int_t code,const char* ) const
 {
   //assert(code==1);
   Double_t y=0.;
   if(code==1) {
     cout<<"Using Analytic Integral 1"<<endl;
     y=A1()*ang1_Inte()*t1()*eff_t()+
       A2()*ang2_Inte()*t2()*eff_t()+
       A3()*ang3_Inte()*t3()*eff_t()+
       A4()*ang4_Inte()*t4()*eff_t()+
       A5()*ang5_Inte()*t5()*eff_t()+
       A6()*ang6_Inte()*t6()*eff_t();
     y *= frac_flav();
   }
   
   if(code==2){
     //another integral with anaytic parts for theta2, kai and time
     cout<<"Using Analytic Integral 2"<<endl;
      y=A1()*cos(theta1)*cos(theta1)*sin(theta1)*(8./3)*(6.28318)*t1()*eff_t()+
	A2()*sin(theta1)*sin(theta1)*sin(theta1)*(4./3)*(6.28318)*t2()*eff_t()+
	A3()*sin(theta1)*sin(theta1)*sin(theta1)*(4./3)*(6.28318)*t3()*eff_t()+
	A4()*(0)*t4()*eff_t()+
	A5()*(0)*t5()*eff_t()+
	A6()*(0)*t6()*eff_t();
      y *= frac_flav();
   }
   if(code==3){
     cout<<"Using Analytic Integral 3"<<endl;
     y=A1()*cos(theta2)*cos(theta2)*sin(theta2)*(8./3)*(6.28318)*t1()*eff_t()+
	A2()*sin(theta2)*sin(theta2)*sin(theta2)*(4./3)*(6.28318)*t2()*eff_t()+
	A3()*sin(theta2)*sin(theta2)*sin(theta2)*(4./3)*(6.28318)*t3()*eff_t()+
	A4()*(0)*t4()*eff_t()+
	A5()*(0)*t5()*eff_t()+
	A6()*(0)*t6()*eff_t();
      y *= frac_flav();
}
   if(code==4){
	   cout << "Using Analytic Integral 4" << endl;
     y=A1()*(16./9)*t1()*eff_t()+
       A2()*(16./9)*(1+cos(2*kai))*t2()*eff_t()+
       A3()*(16./9)*(1-cos(2*kai))*t3()*eff_t()+
       A4()*(-32./9)*(sin(2*kai))*t4()*eff_t()+
       A5()*(sqrt(2))*(0)*cos(kai)*t5()*eff_t()+
       A6()*(-sqrt(2))*(0)*sin(kai)*t6()*eff_t();
     y *= frac_flav();
   }

return y;
 }

//need a way to implement different misatg rates for B and Bbar.IF statement needed somewhere.

 Double_t RooBs2PhiPhiFullPdf::t1() const
 {

   return (1+cos(phi_s))*ExpL()+
          (1-cos(phi_s))*ExpH()+
          2*ExpSin()*sin(phi_s)*(1-wtageff()*2)*tageff();
 }


Double_t RooBs2PhiPhiFullPdf::t2() const
 {
   return (1+cos(phi_s))*ExpL()+
          (1-cos(phi_s))*ExpH()+
          2*ExpSin()*sin(phi_s)*(1-wtageff()*2)*tageff();
 }


Double_t RooBs2PhiPhiFullPdf::t3() const
 {
   return (1-cos(phi_s))*ExpL()+
          (1+cos(phi_s))*ExpH()-
          2*ExpSin()*sin(phi_s)*(1-wtageff()*2)*tageff();
 }


Double_t RooBs2PhiPhiFullPdf::t4() const
 {
   return 2* (sin(sp1)*ExpCos()-cos(sp1)*ExpSin()*cos(phi_s)) *(1-wtageff()*2)*tageff()-
          (ExpH()-ExpL()) *cos(sp1)*sin(phi_s);
 }


 Double_t RooBs2PhiPhiFullPdf::t5() const
 {
   return ( (1+cos(phi_s))*ExpL()+
            (1-cos(phi_s))*ExpH()+
            2*ExpSin()*sin(phi_s)*(1-wtageff()*2)*tageff()
          )*cos(sp1-sp2);
 }

 Double_t RooBs2PhiPhiFullPdf::t6() const
 {
   return 2* ( sin(sp2)*ExpCos()-cos(sp2)*ExpSin()*cos(phi_s)) *(1-wtageff()*2)*tageff() -
          (ExpH()-ExpL())*cos(sp2)*sin(phi_s);

 }

 RooComplex RooBs2PhiPhiFullPdf::evalCerfApprox(Double_t swt, Double_t u, Double_t c) const
 {
   // use the approximation: erf(z) = exp(-z*z)/(sqrt(pi)*z)
   // to explicitly cancel the divergent exp(y*y) behaviour of
   // CWERF for z = x + i y with large negative y

   static Double_t rootpi= sqrt(atan2(0.,-1.));
   RooComplex z(swt*c,u+c);  
   RooComplex zc(u+c,-swt*c);
   RooComplex zsq= z*z;
   RooComplex v= -zsq - u*u;

   return v.exp()*(-zsq.exp()/(zc*rootpi) + 1)*2 ; //why shoule be a 2 here?
   //   return v.exp()*(-zsq.exp()/(zc*rootpi) + 1);
 }

 Double_t RooBs2PhiPhiFullPdf::ExpH() const 
 {
   Double_t v=0;
   if(sigtime>0) {     
     Double_t x = Gamma_H*(time-biastime) ;
     Double_t c = Gamma_H*sigtime/sqrt(2.); 
     Double_t u = (time-biastime)/(sqrt(2.)*sigtime);
     v= exp(c*c-x) * RooMath::erfc(c-u);
     v /=2 ;
   } else if(sigtime<0) {
     double tl=time-biastime-5*sigtime;
     double th=time-biastime+5*sigtime;
     if(tl<0) tl=0;
     if(th<tl+5*sigtime) th=tl+5*sigtime;
     int n=100;
     double dx=(th-tl)/(n-1);
     for (int i=0; i<n;i++) {
       double xi=tl+i*dx;
       v+= exp(-Gamma_H*xi)*reso_t(time-xi);
     }
     v*=dx;
   } else {
     v = exp(-Gamma_H*(time-biastime));
   }

  Double_t v1=0;
   if(sigtime>0) {     
     Double_t x = Gamma_H*(time-biastime) ;
     Double_t c = Gamma_H*sigtime2/sqrt(2.); 
     Double_t u = (time-biastime)/(sqrt(2.)*sigtime2);
     v1= exp(c*c-x) * RooMath::erfc(c-u);
     v1 /=2 ;
   }

   //return v;
   return (0.833*v + 0.167*v1);
 }


 Double_t RooBs2PhiPhiFullPdf::ExpL() const
 {
   Double_t v=0;
   if(sigtime>0) {
     Double_t x = Gamma_L*(time-biastime) ;
     Double_t c = Gamma_L*sigtime/sqrt(2.);
     Double_t u = (time-biastime)/(sqrt(2.)*sigtime);
     v= exp(c*c-x) * RooMath::erfc(c-u);
     v /= 2;
   } else if(sigtime<0) {
     double tl=(time-biastime)-5*sigtime;
     double th=(time-biastime)+5*sigtime;
     if(tl<0) tl=0;
     if(th<tl+5*sigtime) th=tl+5*sigtime;
     int n=100;
     double dx=(th-tl)/(n-1);
     for (int i=0; i<n;i++) {
       double xi=tl+i*dx;
       v+= exp(-Gamma_L*xi)*reso_t(time-xi);
     }
     v/=n;
   } else {
     v = exp(-Gamma_L*(time-biastime));
   }

  Double_t v1=0;
   if(sigtime>0) {
     Double_t x = Gamma_L*(time-biastime) ;
     Double_t c = Gamma_L*sigtime2/sqrt(2.);
     Double_t u = (time-biastime)/(sqrt(2.)*sigtime2);
     v1= exp(c*c-x) * RooMath::erfc(c-u);
     v1 /= 2;
   }

   //   return v;
   return (0.833*v + 0.167*v1); 
 }


 Double_t RooBs2PhiPhiFullPdf::ExpSin() const
 {
   Double_t v=0;
   if(sigtime>0) {
     double g = (Gamma_L+Gamma_H)/2;
     Double_t x = g*(time-biastime) ;
     Double_t c = g*sigtime/sqrt(2.);
     Double_t u = (time-biastime)/(sqrt(2.)*sigtime);
     Double_t wt = dms/g ;
     v += +1*evalCerfIm(wt,-u,c) ;
     v += -1*evalCerfIm(-wt,-u,c) ;
     v /= 4;
   } else if(sigtime<0) {
     double tl=time-biastime-5*sigtime;
     double th=time-biastime+5*sigtime;
     if(tl<0) tl=0;
     if(th<tl+5*sigtime) th=tl+5*sigtime;
     int n=100;
     double dx=(th-tl)/(n-1);
     for (int i=0; i<n;i++) {
       double xi=tl+i*dx;
       v+= exp(-(Gamma_H+Gamma_L)/2*xi)*sin(dms*xi)*reso_t(time-xi);
     }
     v/=n;
   } else {
     v = exp(-(Gamma_H+Gamma_L)/2*(time-biastime))*sin(dms*(time-biastime));
   }

   Double_t v1=0;
   if(sigtime>0) {
     double g = (Gamma_L+Gamma_H)/2;
     Double_t x = g*(time-biastime) ;
     Double_t c = g*sigtime2/sqrt(2.);
     Double_t u = (time-biastime)/(sqrt(2.)*sigtime2);
     Double_t wt = dms/g ;
     v1 += +1*evalCerfIm(wt,-u,c) ;
     v1 += -1*evalCerfIm(-wt,-u,c) ;
     v1 /= 4;
   }

   return (0.833*v + 0.167*v1); 
   //   return v;
 }


Double_t RooBs2PhiPhiFullPdf::ExpCos() const
 {
   Double_t v=0;
   if(sigtime>0) {
     double g = (Gamma_L+Gamma_H)/2;
     Double_t x = g*(time-biastime) ;
     Double_t c = g*sigtime/sqrt(2.);
     Double_t u = (time-biastime)/(sqrt(2.)*sigtime);
     Double_t wt = dms/g ;
     v += evalCerfRe(wt,-u,c) ;
     v += evalCerfRe(-wt,-u,c) ;
     v /= 4;
   } else if(sigtime<0) {
     double tl=time-biastime-5*sigtime;
     double th=time-biastime+5*sigtime;
     if(tl<0) tl=0;
     if(th<tl+5*sigtime) th=tl+5*sigtime;
     int n=100;
     double dx=(th-tl)/(n-1);
     for (int i=0; i<n;i++) {
       double xi=tl+i*dx;
       v+= exp(-(Gamma_H+Gamma_L)/2*xi)*cos(dms*xi)*reso_t(time-xi);
     }
     v/=n;
   } else {
     v = exp(-(Gamma_H+Gamma_L)/2*(time-biastime))*cos(dms*(time-biastime));
   }

   Double_t v1=0;
   if(sigtime>0) {
     double g = (Gamma_L+Gamma_H)/2;
     Double_t x = g*(time-biastime) ;
     Double_t c = g*sigtime2/sqrt(2.);
     Double_t u = (time-biastime)/(sqrt(2.)*sigtime2);
     Double_t wt = dms/g ;
     v1 += evalCerfRe(wt,-u,c) ;
     v1 += evalCerfRe(-wt,-u,c) ;
     v1 /= 4;
   }
    return (0.833*v + 0.167*v1); 
   //  return v;
 }


Double_t RooBs2PhiPhiFullPdf::reso_t(Double_t dt) const
 {
   return exp(-(dt-biastime)*(dt-biastime)/2/sigtime/sigtime)/sigtime/sqrt(2*3.1415927);
 }

 Double_t RooBs2PhiPhiFullPdf::reso_theta(Double_t dt, Double_t res) const
 {
   return exp(-(dt)*(dt)/2/res/res)/res/sqrt(2*3.1415927);
 }

 Double_t RooBs2PhiPhiFullPdf::reso_phi(Double_t dt, Double_t res) const
 {
   return exp(-(dt)*(dt)/2/res/res)/res/sqrt(2*3.1415927);
 }

Double_t RooBs2PhiPhiFullPdf::frac_flav() const
 {
   //   Double_t v = 1;
    Double_t v = 1-efftag;
   if(tag==1)  v=(efftag - defftag)/2;
   if(tag==-1) v=(efftag + defftag)/2;
   return v;
 }

Double_t RooBs2PhiPhiFullPdf::CosTheta1() const
{
  Double_t v=0;
  if(sigtheta==0) v=cos(theta1);
  else{
    double tl=theta1-5*sigtheta;
    double th=theta1+5*sigtheta;
    if(tl<0) tl=0;
    if(th<tl+5*sigtheta) th=tl+5*sigtheta;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v+=cos(xi)*reso_theta(theta1-xi, sigtheta);
    }
    v/=n;
  }

  Double_t v1=0;
  if(sigtheta==0) v=cos(theta1);
  else{
    double tl=theta1-5*0.032;
    double th=theta1+5*0.032;
    if(tl<0) tl=0;
    if(th<tl+5*0.032) th=tl+5*0.032;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v1+=cos(xi)*reso_theta(theta1-xi,0.032);
    }
    v1/=n;
  }

  return (0.237*v1 + 0.763*v);
  //  return v;
}

Double_t RooBs2PhiPhiFullPdf::SinTheta1() const
{
  Double_t v=0;
  if(sigtheta==0) v=sin(theta1);
  else{
    double tl=theta1-5*sigtheta;
    double th=theta1+5*sigtheta;
    if(tl<0) tl=0;
    if(th<tl+5*sigtheta) th=tl+5*sigtheta;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v+=sin(xi)*reso_theta(theta1-xi,sigtheta);
    }
    v/=n;
  }

 Double_t v1=0;
  if(sigtheta==0) v=sin(theta1);
  else{
    double tl=theta1-5*0.032;
    double th=theta1+5*0.032;
    if(tl<0) tl=0;
    if(th<tl+5*0.032) th=tl+5*0.032;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v1+=sin(xi)*reso_theta(theta1-xi,0.032);
    }
    v1/=n;
  }
  
  return (0.237*v1 + 0.763*v);
  //  return v;
}

Double_t RooBs2PhiPhiFullPdf::CosTheta2() const
{
  Double_t v=0;
  if(sigtheta==0) v=cos(theta2);
  else{
  double tl=theta2-5*sigtheta;
  double th=theta2+5*sigtheta;
  if(tl<0) tl=0;
  if(th<tl+5*sigtheta) th=tl+5*sigtheta;
  int n=100;
  double dx=(th-tl)/(n-1);
  for (int i=0; i<n;i++) {
    double xi=tl+i*dx;
    v+=cos(xi)*reso_theta(theta2-xi,sigtheta);
  }
  v/=n;
  }

 Double_t v1=0;
  if(sigtheta==0) v=cos(theta2);
  else{
  double tl=theta2-5*0.032;
  double th=theta2+5*0.032;
  if(tl<0) tl=0;
  if(th<tl+5*0.032) th=tl+5*0.032;
  int n=100;
  double dx=(th-tl)/(n-1);
  for (int i=0; i<n;i++) {
    double xi=tl+i*dx;
    v1+=cos(xi)*reso_theta(theta2-xi,0.032);
  }
  v1/=n;
  }

 return (0.237*v1 + 0.763*v);
  //  return v;
}

Double_t RooBs2PhiPhiFullPdf::SinTheta2() const
{
  Double_t v=0;
  if(sigtheta==0) v=sin(theta2);
  else{
  double tl=theta2-5*sigtheta;
  double th=theta2+5*sigtheta;
  if(tl<0) tl=0;
  if(th<tl+5*sigtheta) th=tl+5*sigtheta;
  int n=100;
  double dx=(th-tl)/(n-1);
  for (int i=0; i<n;i++) {
    double xi=tl+i*dx;
    v+=sin(xi)*reso_theta(theta2-xi,sigtheta);
  }
  v/=n;
  }

  Double_t v1=0;
  if(sigtheta==0) v=sin(theta2);
  else{
    double tl=theta2-5*0.032;
    double th=theta2+5*0.032;
    if(tl<0) tl=0;
    if(th<tl+5*0.032) th=tl+5*0.032;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v1+=sin(xi)*reso_theta(theta2-xi,0.032);
    }
    v1/=n;
  }
 
  return (0.237*v1 + 0.763*v); 
  //return v;
}

Double_t RooBs2PhiPhiFullPdf::Sin2Theta1() const
{
  Double_t v=0;
  if(sigtheta==0) v=sin(2*theta1);
  else{
    double tl=theta1-5*sigtheta;
    double th=theta1+5*sigtheta;
    if(tl<0) tl=0;
    if(th<tl+5*sigtheta) th=tl+5*sigtheta;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v+=sin(2*xi)*reso_theta(theta1-xi,sigtheta);
    }
    v/=n;
  }

  Double_t v1=0;
  if(sigtheta==0) v=sin(2*theta1);
  else{
    double tl=theta1-5*0.032;
    double th=theta1+5*0.032;
    if(tl<0) tl=0;
    if(th<tl+5*0.032) th=tl+5*0.032;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v1+=sin(2*xi)*reso_theta(theta1-xi,sigtheta);
    }
    v1/=n;
  }
  
  return (0.237*v1 + 0.763*v); 
  //  return v;
}


Double_t RooBs2PhiPhiFullPdf::Sin2Theta2() const
{
  Double_t v=0;
  if(sigtheta==0) v=sin(2*theta2);
  else{
    double tl=theta2-5*sigtheta;
    double th=theta2+5*sigtheta;
    if(tl<0) tl=0;
    if(th<tl+5*sigtheta) th=tl+5*sigtheta;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v+=sin(2*xi)*reso_theta(theta2-xi,sigtheta);
    }
    v/=n;
  }

  Double_t v1=0;
  if(sigtheta==0) v=sin(2*theta2);
  else{
    double tl=theta2-5*0.032;
    double th=theta2+5*0.032;
    if(tl<0) tl=0;
    if(th<tl+5*0.032) th=tl+5*0.032;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v1+=sin(2*xi)*reso_theta(theta2-xi,0.032);
    }
    v1/=n;
  }

  return (0.237*v1 + 0.763*v); 
  //  return v;
}

Double_t RooBs2PhiPhiFullPdf::CosPhi() const
{
  Double_t v=0;
  if(sigphi==0) v=cos(kai);
  else{
  double tl=kai-5*sigphi;
  double th=kai+5*sigphi;
  if(tl<0) tl=0;
  if(th<tl+5*sigphi) th=tl+5*sigphi;
  int n=100;
  double dx=(th-tl)/(n-1);
  for (int i=0; i<n;i++) {
    double xi=tl+i*dx;
    v+=cos(xi)*reso_phi(kai-xi,sigphi);
  }
  v/=n;
  }

  Double_t v1=0;
  if(sigphi==0) v=cos(kai);
  else{
  double tl=kai-5*0.0625;
  double th=kai+5*0.0625;
  if(tl<0) tl=0;
  if(th<tl+5*0.0625) th=tl+5*0.0625;
  int n=100;
  double dx=(th-tl)/(n-1);
  for (int i=0; i<n;i++) {
    double xi=tl+i*dx;
    v1+=cos(xi)*reso_phi(kai-xi,0.0625);
  }
  v1/=n;
  }

  return (0.439*v1 + 0.561*v); 
  //  return v;
}

Double_t RooBs2PhiPhiFullPdf::SinPhi() const
{
  Double_t v=0;
   if(sigphi==0) v=sin(kai);
  else{
  double tl=kai-5*sigphi;
  double th=kai+5*sigphi;
  if(tl<0) tl=0;
  if(th<tl+5*sigphi) th=tl+5*sigphi;
  int n=100;
  double dx=(th-tl)/(n-1);
  for (int i=0; i<n;i++) {
    double xi=tl+i*dx;
    v+=sin(xi)*reso_phi(kai-xi,sigphi);
  }
  v/=n;
  }  

 Double_t v1=0;
   if(sigphi==0) v=sin(kai);
   else{
  double tl=kai-5*0.0625;
  double th=kai+5*0.0625;
  if(tl<0) tl=0;
  if(th<tl+5*0.0625) th=tl+5*0.0625;
  int n=100;
  double dx=(th-tl)/(n-1);
  for (int i=0; i<n;i++) {
    double xi=tl+i*dx;
    v1+=sin(xi)*reso_phi(kai-xi,0.0625);
  }
  v1/=n;
  }  
  
   return (0.439*v1 + 0.561*v); 
   //return v;
}

Double_t RooBs2PhiPhiFullPdf::Cos2Phi() const
{
  Double_t v=0;
  if(sigphi==0) v=cos(2*kai);
  else{
    double tl=kai-5*sigphi;
    double th=kai+5*sigphi;
    if(tl<0) tl=0;
    if(th<tl+5*sigphi) th=tl+5*sigphi;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v+=cos(2*xi)*reso_phi(kai-xi,sigphi);
    }
    v/=n;
  } 
  
 Double_t v1=0;
  if(sigphi==0) v=cos(2*kai);
  else{
    double tl=kai-5*0.0625;
    double th=kai+5*0.0625;
    if(tl<0) tl=0;
    if(th<tl+5*0.0625) th=tl+5*0.0625;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v1+=cos(2*xi)*reso_phi(kai-xi,0.0625);
    }
    v1/=n;
  } 
  return (0.439*v1 + 0.561*v); 
  //return v;
}
Double_t RooBs2PhiPhiFullPdf::Sin2Phi() const
{
  Double_t v=0;
  if(sigphi==0) v=sin(2*kai);
  else{
    double tl=kai-5*sigphi;
    double th=kai+5*sigphi;
    if(tl<0) tl=0;
    if(th<tl+5*sigphi) th=tl+5*sigphi;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v+=sin(2*xi)*reso_phi(kai-xi,sigphi);
    }
    v/=n;
  } 
 
  Double_t v1=0;
  if(sigphi==0) v=sin(2*kai);
  else{
    double tl=kai-5*0.0625;
    double th=kai+5*0.0625;
    if(tl<0) tl=0;
    if(th<tl+5*0.0625) th=tl+5*0.0625;
    int n=100;
    double dx=(th-tl)/(n-1);
    for (int i=0; i<n;i++) {
      double xi=tl+i*dx;
      v1+=sin(2*xi)*reso_phi(kai-xi,0.0625);
    }
    v1/=n;
  } 
  
  return (0.439*v1 + 0.561*v);
  //return v;
}

