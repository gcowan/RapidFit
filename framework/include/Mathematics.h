/**
  @namespace Mathematics

  Namespace holding some common mathematical functions that are used
  by many PDFs. For example, gaussian's convolved with trigonometric
  functions.

  @author Greig A Cowan, Pete Clalrke greig.cowan@cern.ch
  @date 2009-12-22
 */
#ifndef RAPIDFIT_MATHEMATICS
#define RAPIDFIT_MATHEMATICS
#include "IDataSet.h"
#include "IPDF.h"
#include "RapidFitIntegrator.h"
#include "PhaseSpaceBoundary.h"
#include "TMath.h"

namespace Mathematics
{
	//	This comes up A LOT, save some CPU cycles and cache it!
	static const double sqrt_2 = sqrt(2.);
	static const double _over_sqrt_2 = 1./sqrt_2;
	static const double invpi=TMath::InvPi();
	static const double global_frac = 0.28125 * invpi;
	static const double third= 1./3.;
	static const double root_6 = sqrt(6.);
	static const double root_3 = sqrt(3.);
	static double rootpi= sqrt(atan2(0.,-1.));

	inline double SQRT_2(){ return sqrt_2; }//=sqrt(2.);
	inline double _Over_SQRT_2(){ return _over_sqrt_2; }// = 1./SQRT_2;
	inline double _Over_PI(){ return invpi; }// = TMath::InvPi();//1./TMath::Pi();
	inline double Global_Frac(){ return global_frac; }// = 0.28125 * _Over_PI;
	inline double Root_6(){ return root_6; }
	inline double Root_3(){ return root_3; }
	inline double Third(){ return third; }
	inline double Rootpi(){ return rootpi; }

	//--------------------------- exp and exp*sin and exp*cos time functions -------------------------
	// time functions for use in PDFs with resolution

	double Exp( const double t, const double gamma, const double resolution );
	double ExpInt( const double tlow, const double thigh, const double gamma, const double resolution  );

	//Added thes to include an upper time acceptance of form ( 1 - beta*t)
	double Exp_betaAcceptance( const double t, const double gamma, const double resolution, const double betaParameter );
	double ExpInt_betaAcceptance( const double tlow, const double thigh, const double gamma, const double resolution, const double betaParameter  );

	double ExpCosh( const double t, const double gamma, const double deltaGamma, const double resolution );
	double ExpCoshInt( const double tlow, const double thigh, const double gamma, const double deltaM, const double resolution );

	double ExpSinh( const double t, const double gamma, const double deltaGamma, const double resolution );
	double ExpSinhInt( const double tlow, const double thigh, const double gamma, const double deltaM, const double resolution );

	double ExpCos( const double t, const double gamma, const double deltaM, const double resolution );
	double ExpCosInt( const double tlow, const double thigh, const double gamma, const double deltaM, const double resolution );

	double ExpSin( const double t, const double gamma, const double deltaM, const double resolution );
	double ExpSinInt( const double tlow, const double thigh, const double gamma, const double deltaM, const double resolution );

	double expErfInt(const double tlimit, const double tau, const double sigma);
	void getBs2JpsiPhiAngularFunctions( double & f1, double & f2, double & f3, double & f4, double & f5, double & f6, const double cosTheta, const double phi, const double cosPsi);
	void getBs2JpsiPhiAngularFunctionsWithSwave( double & f1, double & f2, double & f3, double & f4, double & f5, double & f6, double & f7, double & f8, double & f9, double & f10, const double cosTheta, const double phi, const double cosPsi);
	void getBs2PhiPhiAngularFunctions( double & f1, double & f2, double & f3, double & f4, double & f5, double & f6, const double ct1, const double ct2, const double phi);

  	void calculateAcceptanceWeights( IDataSet * dataSet, IPDF * PDF );
	void calculateAcceptanceWeightsWithSwave( IDataSet * dataSet, IPDF * PDF );
}

#endif
