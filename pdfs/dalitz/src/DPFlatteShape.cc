#include "DPFlatteShape.hh"
#include <complex>
DPFlatteShape::DPFlatteShape(const double in_mean, const double in_g0, const double in_m0a, const double in_m0b, const double in_g1, const double in_m1a, const double in_m1b) :
	  mean(in_mean)
	, g0(in_g0) 
	, m0a(in_m0a)
	, m0b(in_m0b)
	, g1(in_g1)
	, m1a(in_m1a)
	, m1b(in_m1b)
{
}

DPFlatteShape::DPFlatteShape(const DPFlatteShape& other) : DPMassShape(other)
	, mean(other.mean)
	, g0(other.g0)
	, m0a(other.m0a)
	, m0b(other.m0b)
	, g1(other.g1)
	, m1a(other.m1a)
	, m1b(other.m1b)
{
}

std::complex<double> DPFlatteShape::massShape(const double x) const
{
	if (g0<0 || g1<0)
		return std::complex<double>(0,0);
	double s = x*x;
	// Energy, centre of mass p^2 of first channel
	double E0a = 0.5 * (s + m0a*m0a - m0b*m0b) / x;
	double qSq0 = E0a*E0a - m0a*m0a; 
	// Energy, centre of mass p^2 of second channel
	double E1a = 0.5 * (s + m1a*m1a - m1b*m1b) / x;
	double qSq1 = E1a*E1a - m1a*m1a; 
	std::complex<double> gamma0 = (qSq0 > 0) ? std::complex<double>(g0*sqrt(qSq0),0) : std::complex<double>(0, g0*sqrt(-qSq0));
	std::complex<double> gamma1 = (qSq1 > 0) ? std::complex<double>(g1*sqrt(qSq1),0) : std::complex<double>(0, g1*sqrt(-qSq1));
	std::complex<double> gamma = gamma0 + gamma1;
	std::complex<double> partA(mean*mean - s, 0);
	std::complex<double> partB = std::complex<double>(0.0, 2*mean/x) * gamma;
	std::complex<double> denom = partA - partB;
	std::complex<double> T(1,0);
	T = T / denom;
	return T;
}

void DPFlatteShape::setParameters(const std::vector<double>& pars)
{
	setResonanceParameters(pars[0],pars[1],pars[2]);
}

void DPFlatteShape::setResonanceParameters(const double in_mean, const double in_g0, const double in_g1)
{
	mean=in_mean;
	g0=in_g0;
	g1=in_g1;
}
