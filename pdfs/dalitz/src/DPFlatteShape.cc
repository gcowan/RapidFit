#include "DPFlatteShape.hh"
#include "TComplex.h"
DPFlatteShape::DPFlatteShape(double in_mean, double in_g0, double in_m0a, double in_m0b, double in_g1, double in_m1a, double in_m1b) : 
    mean(in_mean)
  , g0(in_g0) 
  , m0a(in_m0a)
  , m0b(in_m0b)
  , g1(in_g1)
  , m1a(in_m1a)
  , m1b(in_m1b)
{
}
DPFlatteShape::DPFlatteShape(const DPFlatteShape& other) : 
    DPMassShape(other)
  , mean(other.mean)
  , g0(other.g0)
  , m0a(other.m0a)
  , m0b(other.m0b)
  , g1(other.g1)
  , m1a(other.m1a)
  , m1b(other.m1b)
{
}
DPFlatteShape::~DPFlatteShape()
{
}
TComplex DPFlatteShape::massShape(double x)
{
  if (g0<0 || g1<0)
  {
  return TComplex(0,0);
  }
  double s = x*x;
  // Energy, centre of mass p^2 of first channel
  double E0a = 0.5 * (s + m0a*m0a - m0b*m0b) / x;
  double qSq0 = E0a*E0a - m0a*m0a; 
  // Energy, centre of mass p^2 of second channel
  double E1a = 0.5 * (s + m1a*m1a - m1b*m1b) / x;
  double qSq1 = E1a*E1a - m1a*m1a; 
  TComplex gamma0 = (qSq0 > 0) ? TComplex(g0*sqrt(qSq0),0) : TComplex(0, g0*sqrt(-qSq0));
  TComplex gamma1 = (qSq1 > 0) ? TComplex(g1*sqrt(qSq1),0) : TComplex(0, g1*sqrt(-qSq1));
  TComplex gamma = gamma0 + gamma1;
  TComplex partA(mean*mean - s, 0);
  TComplex partB = TComplex(0.0, 2*mean/x) * gamma;
  TComplex denom = partA - partB;
  TComplex T(1,0);
  T = T / denom;
  return T;
}

void DPFlatteShape::setParameters(double* pars)
{
  setResonanceParameters(pars[0],pars[1],pars[2]);
}

void DPFlatteShape::setResonanceParameters(double in_mean, double in_g0, double in_g1)
{
  mean=in_mean;
  g0=in_g0;
  g1=in_g1;
}

void DPFlatteShape::setResonanceParameters(double n, double o)
{
  (void)n;
  (void)o;
  return;
}

