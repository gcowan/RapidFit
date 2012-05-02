#ifndef DP_TOTALAMPLITUDE_HH
#define DP_TOTALAMPLITUDE_HH

#include <vector>

#include "DPComponent.hh"
#include "DPWignerFunctionJ1.hh"

class DPTotalAmplitude
{
  public:

    DPTotalAmplitude();
    ~DPTotalAmplitude();

    double matrixElement(double m23, double cosTheta1, double cosTheta2, double phi);

  private:

    std::vector<DPComponent*> KpiComponents;
    std::vector<DPComponent*> ZComponents;

    DPWignerFunctionJ1 wigner;

    double mJpsi;
};

#endif
