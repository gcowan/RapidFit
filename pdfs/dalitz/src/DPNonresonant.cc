#include "DPNonresonant.hh"

#include <iostream>

std::complex<double> DPNonresonant::massShape(const double m) const
{
	(void)m;
	std::complex<double> result(1,0);
	return result;
}

