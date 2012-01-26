#ifndef Bs2PhiPhi_AccCorr_H
#define Bs2PhiPhi_AccCorr_H

class Bs2PhiPhi_AccCorr
{
	public:
	static int    gNVariables;
	static int    gNCoefficients;
	static double gDMean;
	static double gXMean[3];
	static double gXMin[3];
	static double gXMax[3];
	static double gCoefficient[35];
	static double gCoefficientRMS[35];
	static int    gPower[35*3];
	double MDF(double*);

};

#endif
