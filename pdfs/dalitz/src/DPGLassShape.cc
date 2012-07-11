#include "DPGLassShape.hh"
#include "DPBarrierL0.hh"
#include "DPBarrierL1.hh"
#include "DPBarrierL2.hh"
#include "DPBarrierL3.hh"

#include <iostream>

#define DOUBLE_TOLERANCE 1E-6

DPGLassShape::DPGLassShape(double mR, double gammaR, int L,
                    double m1, double m2, double R, double a, double r):
   mR(mR)
  ,gammaR(gammaR)
  ,LR(L)
  ,m1(m1)
  ,m2(m2)
  ,a(a)
  ,r(r)
  ,fraction(0.5)
  ,phaseR(0)
  ,phaseB(0)
{
  switch (L)
  {
    case 0: barrier=new DPBarrierL0(R);
            break;
    default: std::cout<<"WARNING: Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             barrier=new DPBarrierL0(R);
             break;
  }

  pR0=daughterMomentum(mR);
}

DPGLassShape::~DPGLassShape()
{
  if (barrier)
  {
    delete barrier;
  }
}

TComplex DPGLassShape::massShape(double m)
{
// Calculate delta_R
  double tanDeltaR=mR*gamma(m)/(mR*mR-m*m);
  double deltaR=0;
  if ( (mR-m) < DOUBLE_TOLERANCE )
  {
    deltaR=TMath::Pi()/2.0;
  }
  else
  {
    deltaR=TMath::ATan(tanDeltaR);
  }

// Calculate delta_B
  double q=daughterMomentum(m);
  double cotDeltaB=1./(a*q)+r*q/2.;
  double deltaB=0;
  if ( q < DOUBLE_TOLERANCE )
  {
    if (a>0)
    {
      deltaB=0;
    }
    else
    {
      deltaB=TMath::Pi();
    }
  }
  else if (cotDeltaB < DOUBLE_TOLERANCE || TMath::IsNaN(1./cotDeltaB))
  {
    deltaB=TMath::Pi()/2.0;
  }
  else
  {
    deltaB=TMath::ATan(1.0/cotDeltaB);
  }

//  double sinDeltas=TMath::Sin(deltaR+deltaB);
//  double cosDeltas=TMath::Cos(deltaR+deltaB);
//  TComplex result(sinDeltas*cosDeltas,sinDeltas*sinDeltas);
  TComplex result(0,0);
  result+=fraction*TComplex::Sin(deltaR)*TComplex::Exp(TComplex::I()*(deltaR+phaseR))*
          TComplex::Exp(2.0*TComplex::I()*(deltaB+phaseB));
  result+=(1-fraction)*TComplex::Sin(deltaB+phaseB)*TComplex::Exp(TComplex::I()*(deltaB+phaseB));

  return result;
}

double DPGLassShape::gamma(double m)
{
  double pp=daughterMomentum(m);  // momentum of daughter at the actual mass
  double bb=barrier->barrier(pR0,pp);  // Barrier factor
  double gg=gammaR*mR/m*bb*bb*TMath::Power(pp/pR0,2*LR+1);

  return gg;
}

double DPGLassShape::daughterMomentum(double m)
{
  double momentum;

  momentum=(m*m-(m1+m2)*(m1+m2))*(m*m-(m1-m2)*(m1-m2));
  momentum=TMath::Sqrt(momentum);
  momentum/=2*m;

  return momentum;
}

void DPGLassShape::setParameters(double* pars)
{
  mR=pars[0];
  gammaR=pars[1];
  a=pars[2];
  r=pars[3];
  fraction=pars[4];
  phaseR=pars[5];  
  phaseB=pars[6];  
}

void DPGLassShape::setResonanceParameters(double a_lass, double r_lass)
{
	a = a_lass;
	r = r_lass;
	return;
}
