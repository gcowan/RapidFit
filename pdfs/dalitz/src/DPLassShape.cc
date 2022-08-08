#include "DPLassShape.hh"
#include "DPHelpers.hh"
#include <iostream>

#define DOUBLE_TOLERANCE 1E-6

DPLassShape::DPLassShape(double mRR, double gammaRR, int L, double mm1, double mm2, double RR, double aa, double rr):
	 mR(mRR)
	,gammaR(gammaRR)
	,LR(L)
	,m1(mm1)
	,m2(mm2)
	,a(aa)
	,r(rr)
{
	pR0=DPHelpers::daughterMomentum(mR,m1,m2);
	barrier = DPBarrierFactor((unsigned int)LR,RR,pR0);
}

DPLassShape::DPLassShape(const DPLassShape& other) : DPMassShape(other)
	,mR(other.mR)
	,gammaR(other.gammaR)
	,LR(other.LR)
	,m1(other.m1)
	,m2(other.m2)
	,a(other.a)
	,r(other.r)
	,barrier(other.barrier)
	,pR0(other.pR0)
{
}

std::complex<double> DPLassShape::massShape(const double m) const
{
// Calculate delta_R
	double tanDeltaR=mR*gamma(m)/(mR*mR-m*m);
	double deltaR=0;
	if ( (mR-m) < DOUBLE_TOLERANCE )
		deltaR=M_PI/2.0;
	else
		deltaR=std::atan(tanDeltaR);
// Calculate delta_B
	double q=DPHelpers::daughterMomentum(m,m1,m2);
	double cotDeltaB=1./(a*q)+r*q/2.;
	double deltaB=0;
	if ( q < DOUBLE_TOLERANCE )
	{
		if (a>0)
			deltaB=0;
		else
			deltaB=M_PI;
	}
	else if (cotDeltaB < DOUBLE_TOLERANCE || std::isnan(1./cotDeltaB))
		deltaB=M_PI/2.0;
	else
		deltaB=std::atan(1.0/cotDeltaB);
	double sinDeltas=std::sin(deltaR+deltaB);
	double cosDeltas=std::cos(deltaR+deltaB);
	std::complex<double> result(sinDeltas*cosDeltas,sinDeltas*sinDeltas);
	return result;
}

double DPLassShape::gamma(const double m) const
{
	double pp=DPHelpers::daughterMomentum(m,m1,m2);  // momentum of daughter at the actual mass
	double bb=barrier.barrier_sq(pp);  // Barrier factor
	double gg=gammaR*mR/m*bb*std::pow(pp/pR0,2*LR+1);
	return gg;
}

void DPLassShape::setParameters(const std::vector<double>& pars)
{
	setResonanceParameters(pars[0],pars[1]);
	return;
}
void DPLassShape::setResonanceParameters(const double a_lass, const double r_lass)
{
	a = a_lass;
	r = r_lass;
	return;
}

