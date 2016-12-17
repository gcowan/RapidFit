#ifndef DPFlatteShape_H
#define DPFlatteShape_H
#include "DPMassShape.hh"
#include "TComplex.h"
class DPFlatteShape: public virtual DPMassShape
{

  public:

    DPFlatteShape(double in_mean, double in_g0, double in_m0a, double in_m0b, double in_g1, double in_m1a, double in_m1b);
    DPFlatteShape( const DPFlatteShape& );
    ~DPFlatteShape();

    TComplex massShape(double x);

    void setParameters(double* pars);
  
  protected:
    double mean, g0, m0a, m0b, g1, m1a, m1b;
  private:
    void setResonanceParameters( double mass, double sigma ); // Do not use!!
    void setResonanceParameters(double in_mean, double in_g0, double in_g1);
};
#endif

