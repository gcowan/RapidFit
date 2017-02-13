#ifndef DPFlatteShape_H
#define DPFlatteShape_H
#include "DPMassShape.hh"
#include <complex>
class DPFlatteShape: public virtual DPMassShape
{
	public:
		DPFlatteShape(double in_mean, double in_g0, double in_m0a, double in_m0b, double in_g1, double in_m1a, double in_m1b);
		DPFlatteShape(const DPFlatteShape&);
		~DPFlatteShape() {}
		std::complex<double> massShape(const double x) const;
		void setParameters(const std::vector<double>& pars);
	protected:
		double mean, g0, m0a, m0b, g1, m1a, m1b;
	private:
		void setResonanceParameters(const double in_mean, const double in_g0, double in_g1);
};
#endif

