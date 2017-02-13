#ifndef DP_MASS_SHAPE
#define DP_MASS_SHAPE
#include <complex>
#include <vector>
class DPMassShape
{
	public:
		DPMassShape() {};
		DPMassShape(const DPMassShape&) {};
		~DPMassShape() {};
		virtual std::complex<double> massShape(const double m) const = 0;
		virtual void setParameters(const std::vector<double>& pars) = 0;
};
#endif
