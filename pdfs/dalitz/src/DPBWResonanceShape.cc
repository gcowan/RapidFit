#include "DPBWResonanceShape.hh"
#include "DPHelpers.hh"
#include <iostream>

DPBWResonanceShape::DPBWResonanceShape(double mRR, double gammaRR, int L, double mm1, double mm2, double RR):
	 mR(mRR)
	,gammaR(gammaRR)
	,LR(L)
	,m1(mm1)
	,m2(mm2)
{
	pR0=DPHelpers::daughterMomentum(mR,m1,m2);
	barrier = DPBarrierFactor((unsigned int)LR,RR,pR0);
}

DPBWResonanceShape::DPBWResonanceShape( const DPBWResonanceShape& other ) : DPMassShape(other)
	,mR(other.mR)
	,gammaR(other.gammaR)
	,LR(other.LR)
	,m1(other.m1)
	,m2(other.m2)
	,pR0(other.pR0)
	,barrier(other.barrier)
{
}

std::complex<double> DPBWResonanceShape::massShape(const double m) const
{
	double width = gamma(m);
	std::complex<double> result(1,0);
	std::complex<double> denominator(mR*mR-m*m,-mR*width);
	result/=denominator;
	return result;
}

double DPBWResonanceShape::gamma(const double m) const
{
	double pp=DPHelpers::daughterMomentum(m,m1,m2);;  // momentum of daughter at the actual mass
	double bb=barrier.barrier_sq(pp);  // Barrier factor
	double gg=gammaR*(mR/m)*bb*std::pow(pp/pR0,2*LR+1);
	if(std::isnan(pp)) std::cerr << "\t\tDaughter momentum is nan" << std::endl;
	if(std::isnan(bb)) std::cerr << "\t\tBarrier factor is nan" << std::endl;
	if(std::isnan(gg)) std::cerr << "\t\tMass-dependent width is nan" << std::endl;
	if(LR!=barrier.getspin()) std::cerr << "\t\tStored spin of resonance and barrier factor do not match" << std::endl;
	return gg;
}

void DPBWResonanceShape::setParameters(const std::vector<double>& pars)
{
	setResonanceParameters(pars[0],pars[1]);
	pR0=DPHelpers::daughterMomentum(mR,m1,m2);
	barrier.setparameters(pars[2],pR0);
}

void DPBWResonanceShape::setResonanceParameters(const double mass, const double sigma )
{
	mR = mass;
	gammaR = sigma;
}

