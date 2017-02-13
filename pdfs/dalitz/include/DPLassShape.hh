#ifndef DP_LASS_SHAPE
#define DP_LASS_SHAPE
#include "DPMassShape.hh"
#include "DPBarrierFactor.hh"
#include <complex>
class DPLassShape: public virtual DPMassShape
{
	public:
		DPLassShape(double mR, double gammaR, int L, double m1,double m2, double R, double a, double r);
		DPLassShape(const DPLassShape&);
		~DPLassShape() {}
		std::complex<double> massShape(const double m) const;
		void setResonanceParameters(const double a, const double r);
		void setParameters(const std::vector<double>& pars);
	private:
		void Init();
		double mR;
		double gammaR;
		int LR;
		double m1;
		double m2;
		double a;
		double r;
		DPBarrierFactor barrier;
		double pR0;  // Momentum of daughters at mR
		double gamma(const double m) const;
};
#endif

