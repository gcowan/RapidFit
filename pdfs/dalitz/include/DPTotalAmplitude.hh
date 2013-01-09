#ifndef DP_TOTALAMPLITUDE_HH
#define DP_TOTALAMPLITUDE_HH

#include <vector>

#include "DPComponent.hh"
#include "DPWignerFunctionJ1over2.hh"

class DPTotalAmplitude
{
  public:

    DPTotalAmplitude();
    ~DPTotalAmplitude();

  double matrixElement(double m23, double cosTheta1, double cosTheta2, double phi, int pionID);

  private:

    std::vector<DPComponent*> KpiComponents;
    std::vector<DPComponent*> ZComponents;

    DPWignerFunctionJ1over2 wigner;

    double mJpsi;
};

#endif
