#ifndef DP_NONRESONANT_SHAPE
#define DP_NONRESONANT_SHAPE
#include "DPMassShape.hh"
#include <complex>
class DPNonresonant: public virtual DPMassShape
{
	public:
		DPNonresonant() {}
		DPNonresonant(const DPNonresonant& other) : DPMassShape(other) {}
		~DPNonresonant() {}
		std::complex<double> massShape(const double m) const;
		void setParameters(const std::vector<double>& pars) {(void)pars;};
};
#endif
