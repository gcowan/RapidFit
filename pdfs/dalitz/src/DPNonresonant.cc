#include "DPNonresonant.hh"

#include <iostream>

DPNonresonant::DPNonresonant(double mRR, double gammaRR, int L,
                    double mm1, double mm2, double RR):
   mR(mRR)
  ,gammaR(gammaRR)
  ,LR(L)
  ,m1(mm1)
  ,m2(mm2)
{
}

DPNonresonant::DPNonresonant( const DPNonresonant& other ) : DPMassShape(other),
   mR(other.mR), gammaR(other.gammaR), LR(other.LR), m1(other.m1), m2(other.m2)
{
}

DPNonresonant::~DPNonresonant()
{
}

TComplex DPNonresonant::massShape(double m)
{
  TComplex result(1,0);

  return result;
}


