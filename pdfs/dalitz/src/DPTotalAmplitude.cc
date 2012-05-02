
#include "DPTotalAmplitude.hh"
//#include "DPJpsiKSpinTwo.hh"
//#include "DPJpsiKSpinOne.hh"
//#include "DPJpsiKSpinZero.hh"
#include "DPJpsiKaon.hh"
#include "DPZplusK.hh"
#include "DPHelpers.hh"

#include "TComplex.h"
#include "TLorentzVector.h"

#include <iostream>

DPTotalAmplitude::DPTotalAmplitude()
{
  mJpsi=3.096916;
  // Construct all components we need
  DPComponent* tmp;
  // B0 --> J/psi K*
  tmp=new DPJpsiKaon(0, 1, 5.279, 0.89594, 0.0487, 0.493677,
                         0.13957018, 5.0, 1.5, 3.096916,1);
//  tmp->setHelicityAmplitudes(0.0, 1., 0., 0.0 , 0.0, 0.0); // 0,+,- mag,phase
  tmp->setHelicityAmplitudes(0.75, 0.42, 0.05, 0.0, 3.1, -1.45); // 0,+,- mag,phase
//  KpiComponents.push_back(tmp);
  // B0 --> J/psi K*(1410)
  tmp=new DPJpsiKaon(0, 1, 5.279, 1.414, 0.232, 0.493677,
                         0.13957018, 5.0, 1.5, 3.096916,1);
  tmp->setHelicityAmplitudes(0.6, 0.5, 0.62, 0.0, 3.1, -1.45); // 0,+,- mag,phase
//  KpiComponents.push_back(tmp);
  // B0 --> J/psi K*(1680)
  tmp=new DPJpsiKaon(0, 1, 5.279, 1.717, 0.322, 0.493677,
                         0.13957018, 5.0, 1.5, 3.096916,1);
  tmp->setHelicityAmplitudes(0.075, 0.042, 0.005, 0.0, 3.1, -1.45); // 0,+,- mag,phase
//  KpiComponents.push_back(tmp);
  // B0 --> J/psi K0(1430)
  tmp=new DPJpsiKaon(0, 0, 5.279, 1.425, 0.270, 0.493677,
                         0.13957018, 5.0, 1.5, 3.096916,0);
  tmp->setHelicityAmplitudes(0.075, .0, .0, 0.0, 0.0, 0.0); // 0,+,- mag,phase
//  KpiComponents.push_back(tmp);
  // B0 --> J/psi K2(1430)
  tmp=new DPJpsiKaon(0, 0, 5.279, 1.4324, 0.109, 0.493677,
                         0.13957018, 5.0, 1.5, 3.096916,2);
  tmp->setHelicityAmplitudes(0.075, 0.042, 0.005, 0.0, 3.1, -1.45); // 0,+,- mag,phase
//  KpiComponents.push_back(tmp);

  // Kpi s-wave using LASS
  tmp=new DPJpsiKaon(0, 0, 5.279, 1.425, 0.270, 0.493677,
                     0.13957018, 5.0, 1.5, 3.096916,0,
                     "LASS", 0.00415, 0.00136);
   KpiComponents.push_back(tmp);

  // B0 --> Z+ K-
  tmp=new DPZplusK(0,0,5.279,4.430,0.100,0.493677,
                     0.13957018, 5.0, 1.5, 3.096916,0);
  tmp->setHelicityAmplitudes(0.75, 0.42, 0.05, 0.0, 3.1, -1.45);
//  ZComponents.push_back(tmp);
}

DPTotalAmplitude::~DPTotalAmplitude()
{
  // destroy components
  for (unsigned int i=0;i<KpiComponents.size();++i)
  {
    delete KpiComponents[i];
    KpiComponents[i]=0;
  }
  for (unsigned int i=0;i<ZComponents.size();++i)
  {
    delete ZComponents[i];
    ZComponents[i]=0;
  }
}

double DPTotalAmplitude::matrixElement(double m23, double cosTheta1, double cosTheta2,
                         double phi)
{
  double result=0; 

  TComplex tmp(0,0);

  // Now sum over final state helicities (this is not general code, but
  // knows about internals of components

  for (int twoLambda=-2; twoLambda<=2; twoLambda+=4) // Sum over +-1
  {
    tmp=TComplex(0,0);
    for (int twoLambdaPsi=-2; twoLambdaPsi<=2; twoLambdaPsi+=2) // Sum over -1,0,+1
    {
      for (unsigned int i=0;i<KpiComponents.size();++i) // sum over all components
      {
        tmp+=KpiComponents[i]->amplitude(m23, cosTheta1, cosTheta2, phi,
                                         twoLambda, twoLambdaPsi);
      }
      // Now comes sum over Z+ components and lambdaPsiPrime
      for (unsigned int i=0; i<ZComponents.size();++i)
      {
        // Need angle between reference axis
        TLorentzVector pMuPlus;
        TLorentzVector pMuMinus;
        TLorentzVector pPi;
        TLorentzVector pK;
        DPHelpers::calculateFinalStateMomenta(5.279, m23, mJpsi,
              cosTheta1,  cosTheta2, phi, 0.105, 0.105, 0.13957018, 0.493677,
              pMuPlus, pMuMinus, pPi, pK);
        TLorentzVector pB(0,0,0,5.279);
        // Cos of the angle between psi reference axis
        double cosARefs=DPHelpers::referenceAxisCosAngle(pB, pMuPlus, pMuMinus, pPi, pK);
        // Sum over lambdaPsiPrime
        for (int twoLambdaPsiPrime=-2; twoLambdaPsiPrime<=2; twoLambdaPsiPrime+=2)
        {
          tmp+=wigner.function(cosARefs,twoLambdaPsiPrime/2,twoLambdaPsi/2)*
               ZComponents[i]->amplitude(m23, cosTheta1, cosTheta2, phi,
                                         twoLambda,twoLambdaPsiPrime);
        }
      }
    }
    result+=tmp.Rho2();
  }


  return result;
}

