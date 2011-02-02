/**
  @namespace Mathematics

  Namespace holding some common mathematical functions that are used
  by many PDFs. For example, gaussian's convolved with trigonometric
  functions.

  @author Greig A Cowan, Pete Clalrke greig.cowan@cern.ch
  @date 2009-12-22
 */

#include "IDataSet.h"
#include "IPDF.h"
#include "RapidFitIntegrator.h"
#include "PhaseSpaceBoundary.h"

namespace Mathematics
{

	//--------------------------- exp and exp*sin and exp*cos time functions -------------------------
	// time functions for use in PDFs with resolution

	double Exp( double t, double gamma, double resolution );
	double ExpInt( double tlow, double thigh, double gamma, double resolution  );
	double ExpInt( double tlow, double thigh, double gamma, double resolution, double acceptanceParameter  );  //Overloaded to include aceptance function
	
	double ExpCosh( double t, double gamma, double deltaGamma, double resolution );
	double ExpCoshInt( double tlow, double thigh, double gamma, double deltaM, double resolution );
	
	double ExpSinh( double t, double gamma, double deltaGamma, double resolution );
	double ExpSinhInt( double tlow, double thigh, double gamma, double deltaM, double resolution );
	
	double ExpCos( double t, double gamma, double deltaM, double resolution );
	double ExpCosInt( double tlow, double thigh, double gamma, double deltaM, double resolution );
	
	double ExpSin( double t, double gamma, double deltaM, double resolution );
	double ExpSinInt( double tlow, double thigh, double gamma, double deltaM, double resolution );
	
	double expErfInt(double tlimit, double tau, double sigma);
	void getBs2JpsiPhiAngularFunctions( double & f1, double & f2, double & f3, double & f4, double & f5, double & f6
			, double cosTheta, double phi, double cosPsi);
	void getBs2JpsiPhiAngularFunctionsWithSwave( double & f1, double & f2, double & f3, double & f4, double & f5, double & f6,
			double & f7, double & f8, double & f9, double & f10, double cosTheta, double phi, double cosPsi);
	
  	void calculateAcceptanceWeights( IDataSet * dataSet, IPDF * PDF );
	void calculateAcceptanceWeightsWithSwave( IDataSet * dataSet, IPDF * PDF );
}

