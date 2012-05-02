// $Id: CrystalBall.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class CrystalBall CrystalBall.cpp
 *
 *  RapidFit PDF for Bs mass
 *
 *
 *  @author Pete
 *  @date 2011-07-30
 */

#include "CrystalBall.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"

PDF_CREATOR( CrystalBall );

//Constructor
CrystalBall::CrystalBall(PDFConfigurator* configurator) :
	// Physics parameters
	BasePDF()
	, m0Name	( configurator->getName("m0") )
	, sigmaName	( configurator->getName("sigma") )
	, alphaName	( configurator->getName("alpha") )
	, nName		( configurator->getName("n") )
	// Observables
	, recoMassName	( configurator->getName("mass") )
{
	MakePrototypes();
}

//Make the data point and parameter set
void CrystalBall::MakePrototypes()
{
	// Observables
	allObservables.push_back( recoMassName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( m0Name );
	parameterNames.push_back( sigmaName );
	parameterNames.push_back( alphaName );
	parameterNames.push_back( nName );
	allParameters = ParameterSet(parameterNames);
}

//Destructor
CrystalBall::~CrystalBall()
{
}


double CrystalBall::ApproxErf(double arg) const
{
  static const double erflim = 5.0;
  if( arg > erflim )
    return 1.0;
  if( arg < -erflim )
    return -1.0;

  return RooMath::erf(arg);
}

//Calculate the function value
double CrystalBall::Evaluate(DataPoint * measurement)
{
	// Get the physics parameters
	double m0  = allParameters.GetPhysicsParameter( m0Name )->GetValue();
	double sigma = allParameters.GetPhysicsParameter( sigmaName )->GetValue();
	double alpha = allParameters.GetPhysicsParameter( alphaName )->GetValue();
	double n = allParameters.GetPhysicsParameter( nName )->GetValue();

	// Get the observable
	double mass = measurement->GetObservable( recoMassName )->GetValue();

	double returnValue = 0;

	double t = (mass - m0)/sigma;
	if (alpha < 0) t = -t;

	double absAlpha = fabs((double)alpha);

	if (t >= -absAlpha) {
		returnValue = exp(-0.5*t*t);
	}
	else {
		double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
		double b = n/absAlpha - absAlpha;

		returnValue = a/TMath::Power(b - t, n);
	}

	return returnValue;
}

// Normalisation
double CrystalBall::Normalisation(PhaseSpaceBoundary * boundary)
{
	// Get the physics parameters
	double m0  = allParameters.GetPhysicsParameter( m0Name )->GetValue();
	double sigma = allParameters.GetPhysicsParameter( sigmaName )->GetValue();
	double alpha = allParameters.GetPhysicsParameter( alphaName )->GetValue();
	double n = allParameters.GetPhysicsParameter( nName )->GetValue();

	double result = 0.0;
        double mlow(0.), mhigh(0.);
	IConstraint * massBound = boundary->GetConstraint( recoMassName );
        if ( massBound->GetUnit() == "NameNotFoundError" )
        {
                cerr << "Bound on mass not provided" << endl;
                result = -1.;
        }
        else
        {
                mlow  = massBound->GetMinimum();
                mhigh = massBound->GetMaximum();
	}

	static const double sqrtPiOver2 = 1.2533141373;
	static const double sqrt2 = 1.4142135624;

	//assert(code==1);
	bool useLog = false;

	if( fabs(n-1.0) < 1.0e-05 )
		useLog = true;

	double sig = fabs((double)sigma);

	double tmin = (mlow - m0)/sig;
	double tmax = (mhigh - m0)/sig;

	if(alpha < 0) {
		double tmp = tmin;
		tmin = -tmax;
		tmax = -tmp;
	}

	double absAlpha = fabs((double)alpha);

	if( tmin >= -absAlpha ) {
		result += sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
				- ApproxErf(tmin/sqrt2) );
	}
	else if( tmax <= -absAlpha ) {
		double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
		double b = n/absAlpha - absAlpha;

		if(useLog) {
			result += a*sig*( log(b-tmin) - log(b-tmax) );
		}
		else {
			result += a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
					- 1.0/(TMath::Power(b-tmax,n-1.0)) );
		}
	}
	else {
		double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
		double b = n/absAlpha - absAlpha;

		double term1 = 0.0;
		if(useLog) {
			term1 = a*sig*(  log(b-tmin) - log(n/absAlpha));
		}
		else {
			term1 = a*sig/(1.0-n)*(   1.0/(TMath::Power(b-tmin,n-1.0))
					- 1.0/(TMath::Power(n/absAlpha,n-1.0)) );
		}

		double term2 = sig*sqrtPiOver2*(   ApproxErf(tmax/sqrt2)
				- ApproxErf(-absAlpha/sqrt2) );


		result += term1 + term2;
	}

	return result;
}

