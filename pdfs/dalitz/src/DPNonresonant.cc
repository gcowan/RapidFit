#include "DPNonresonant.hh"

#include <iostream>

DPNonresonant::DPNonresonant()
{
}

DPNonresonant::DPNonresonant(const DPNonresonant& other) : DPMassShape(other)
{
}

DPNonresonant::~DPNonresonant()
{
}

TComplex DPNonresonant::massShape(double m)
{
  (void)m;
  TComplex result(1,0);
  return result;
}

