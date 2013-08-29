#include "DPBWResonanceShape.hh"
#include "DPBarrierL0.hh"
#include "DPBarrierL1.hh"
#include "DPBarrierL2.hh"
#include "DPBarrierL3.hh"
#include "DPBarrierL4.hh"
#include "DPBarrierL5.hh"

#include <iostream>

DPBWResonanceShape::DPBWResonanceShape(double mRR, double gammaRR, int L,
                    double mm1, double mm2, double RR):
   mR(mRR)
  ,gammaR(gammaRR)
  ,LR(L)
  ,m1(mm1)
  ,m2(mm2)
  ,R(RR)
{
  switch (LR)
  {
    case 0: barrier=new DPBarrierL0(R);
            break;
    case 1: barrier=new DPBarrierL1(R);
            break;
    case 2: barrier=new DPBarrierL2(R);
            break;
    case 3: barrier=new DPBarrierL3(R);
            break;
    case 4: barrier=new DPBarrierL4(R);
            break;
    case 5: barrier=new DPBarrierL5(R);
            break;
    default: std::cout<<"WARNING: Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             barrier=new DPBarrierL0(R);
             break;
  }

  pR0=daughterMomentum(mR);
}

DPBWResonanceShape::DPBWResonanceShape( const DPBWResonanceShape& other ) : DPMassShape(other),
   mR(other.mR), gammaR(other.gammaR), LR(other.LR), m1(other.m1), m2(other.m2), R(other.R), barrier(NULL)
{
    std::cout << "In DPBW copy const " << mR << " " << gammaR << " " << LR << " " << R << std::endl;
  	if ( other.barrier != NULL )
	{
  		switch (LR)
  		{
    			case 0: barrier=new DPBarrierL0(R);
            			break;
    			case 1: barrier=new DPBarrierL1(R);
            			break;
    			case 2: barrier=new DPBarrierL2(R);
            			break;
    			case 3: barrier=new DPBarrierL3(R);
            			break;
    			case 4: barrier=new DPBarrierL4(R);
            			break;
    			case 5: barrier=new DPBarrierL5(R);
            			break;
    			default: std::cout<<"WARNING: Do not know which barrier factor to use.  Using L=0 and you should check what are you doing.\n";
             			barrier=new DPBarrierL0(R);
             			break;
  		}
	}
  	pR0=daughterMomentum(mR);
}

DPBWResonanceShape::~DPBWResonanceShape()
{
  //std::cout << "DPBW destructor" << std::endl;
  if (barrier)
  {
    delete barrier;
  }
}

TComplex DPBWResonanceShape::massShape(double m)
{
  //std::cout << "DPBWResonanceShape m" << m << " " << mR << std::endl;
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
  //std::cout << "DPBWResonanceShape g " << R << " " << m << " " << mR << " " << gammaR << " " << pp << " " << pR0 << " " << bb <<  " " << m1 << " " << m2 <<std::endl;

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
	//std::cout << "DPBWResonanceShape setting " << mass << std::endl;
	mR = mass;
	gammaR = sigma;
	pR0 = daughterMomentum(mR);
    //std::cout << "DPBWResonanceShape set" << std::endl;
}

