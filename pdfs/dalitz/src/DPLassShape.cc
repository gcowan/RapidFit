#include "DPLassShape.hh"
#include "DPBarrierL0.hh"
#include "DPBarrierL1.hh"
#include "DPBarrierL2.hh"
#include "DPBarrierL3.hh"

#include <iostream>

#define DOUBLE_TOLERANCE 1E-6

DPLassShape::DPLassShape(double mRR, double gammaRR, int L,
                         double mm1, double mm2, double RR, double aa, double rr):
   mR(mRR)
  ,gammaR(gammaRR)
  ,LR(L)
  ,m1(mm1)
  ,m2(mm2)
  ,R(RR)
  ,a(aa)
  ,r(rr)
{
  switch (LR)
  {
    case 0: barrier=new DPBarrierL0(R);
            break;
    default: std::cout<<"WARNING: Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             barrier=new DPBarrierL0(R);
             break;
  }

  pR0=daughterMomentum(mR);
}

DPLassShape::DPLassShape( const DPLassShape& other ) : DPMassShape( other ),
        mR( other.mR ), gammaR( other.gammaR ), LR(other.LR),
        m1(other.m1), m2(other.m2), a(other.a), r(other.r), R(other.R),
        barrier(NULL)
{
        if ( other.barrier != NULL )
        {
                switch (LR)
                {
                        case 0: barrier=new DPBarrierL0(R);
                        break;
                default: std::cout<<"WARNING: Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
                        barrier=new DPBarrierL0(R);
                        break;
                }
        }       
        pR0=daughterMomentum(mR);
}       

DPLassShape::~DPLassShape()
{
  if (barrier)
  {
    delete barrier;
  }
}

TComplex DPLassShape::massShape(double m)
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

  double sinDeltas=TMath::Sin(deltaR+deltaB);
  double cosDeltas=TMath::Cos(deltaR+deltaB);
  TComplex result(sinDeltas*cosDeltas,sinDeltas*sinDeltas);

  return result;
}

double DPLassShape::gamma(double m)
{
  double pp=daughterMomentum(m);  // momentum of daughter at the actual mass
  double bb=barrier->barrier(pR0,pp);  // Barrier factor
  double gg=gammaR*mR/m*bb*bb*TMath::Power(pp/pR0,2*LR+1);

  return gg;
}

double DPLassShape::daughterMomentum(double m)
{
  double momentum;

  momentum=(m*m-(m1+m2)*(m1+m2))*(m*m-(m1-m2)*(m1-m2));
  momentum=TMath::Sqrt(momentum);
  momentum/=2*m;

  return momentum;
}

void DPLassShape::setResonanceParameters(double a_lass, double r_lass)
{
	a = a_lass;
	r = r_lass;
	return;
}
