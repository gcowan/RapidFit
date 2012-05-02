#include "DPBWResonanceShape.hh"
#include "DPBarrierL0.hh"
#include "DPBarrierL1.hh"
#include "DPBarrierL2.hh"
#include "DPBarrierL3.hh"

#include <iostream>

DPBWResonanceShape::DPBWResonanceShape(double mR, double gammaR, int L,
                    double m1, double m2, double R):
   mR(mR)
  ,gammaR(gammaR)
  ,LR(L)
  ,m1(m1)
  ,m2(m2)
{
  switch (L)
  {
    case 0: barrier=new DPBarrierL0(R);
            break;
    case 1: barrier=new DPBarrierL1(R);
            break;
    case 2: barrier=new DPBarrierL2(R);
            break;
    case 3: barrier=new DPBarrierL3(R);
            break;
    default: std::cout<<"WARNING: Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             barrier=new DPBarrierL0(R);
             break;
  }

  pR0=daughterMomentum(mR);
}

DPBWResonanceShape::~DPBWResonanceShape()
{
  if (barrier)
  {
    delete barrier;
  }
}

TComplex DPBWResonanceShape::massShape(double m)
{
  TComplex result(1,0);
  TComplex denominator(mR*mR-m*m,-mR*gamma(m));

  result/=denominator;
  return result;
}

double DPBWResonanceShape::gamma(double m)
{
  double pp=daughterMomentum(m);  // momentum of daughter at the actual mass
  double bb=barrier->barrier(pR0,pp);  // Barrier factor
  double gg=gammaR*mR/m*bb*bb*TMath::Power(pp/pR0,2*LR+1);

  return gg;
}

double DPBWResonanceShape::daughterMomentum(double m)
{
  double momentum;

  momentum=(m*m-(m1+m2)*(m1+m2))*(m*m-(m1-m2)*(m1-m2));
  momentum=TMath::Sqrt(momentum);
  momentum/=2*m;

  return momentum;
}

void DPBWResonanceShape::setResonanceParameters( double mass, double sigma )
{
	//std::cout << "DPBWResonanceShape setting" << std::endl;
	mR = mass;
	gammaR = sigma;	
	//std::cout << "DPBWResonanceShape set" << std::endl;
}

