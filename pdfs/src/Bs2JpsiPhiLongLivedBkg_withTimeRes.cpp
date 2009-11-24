// $Id: Bs2JpsiPhiLongLivedBkg_withTimeRes.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Bs2JpsiPhiLongLivedBkg_withTimeRes Bs2JpsiPhiLongLivedBkg_withTimeRes.cpp
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#include "Bs2JpsiPhiLongLivedBkg_withTimeRes.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"

//Constructor
Bs2JpsiPhiLongLivedBkg_withTimeRes::Bs2JpsiPhiLongLivedBkg_withTimeRes() : 
	// Observables
	  timeName	( "time" )

	// Physics parameters
	, tau1Name	( "tau_LL1" )
	, tau2Name	( "tau_LL2" )
	, f_LL1Name	( "f_LL1" )
	, sigmaLL1Name	( "sigmaLL1" )
	, sigmaLL2Name	( "sigmaLL2" )
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhiLongLivedBkg_withTimeRes::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );

        //Make the parameter set
        vector<string> parameterNames;
        parameterNames.push_back( tau1Name );
        parameterNames.push_back( tau2Name );
        parameterNames.push_back( f_LL1Name );
        parameterNames.push_back( sigmaLL1Name );
        parameterNames.push_back( sigmaLL2Name );
        allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2JpsiPhiLongLivedBkg_withTimeRes::~Bs2JpsiPhiLongLivedBkg_withTimeRes()
{
}

//Calculate the function value
double Bs2JpsiPhiLongLivedBkg_withTimeRes::Evaluate(DataPoint * measurement)
{
	// Observable
        double time = measurement->GetObservable( timeName )->GetValue();

	// Do we need the angular normalisation fraction or will this automatically
	// be calculated?
	double angularNorm = 1./( 2. * 2. * 2.*TMath::Pi() );

  	double tau1 = allParameters.GetPhysicsParameter( tau1Name )->GetValue();
  	double tau2 = allParameters.GetPhysicsParameter( tau2Name )->GetValue();
  	double f_LL1 = allParameters.GetPhysicsParameter( f_LL1Name )->GetValue();
  	double sigmaLL1 = allParameters.GetPhysicsParameter( sigmaLL1Name )->GetValue();
  	double sigmaLL2 = allParameters.GetPhysicsParameter( sigmaLL2Name )->GetValue();
  
	double R1 = erfc( time, tau1, sigmaLL1); 
	double R2 = erfc( time, tau2, sigmaLL2); 
	
	double val = f_LL1 * R1 + (1. - f_LL1) * R2;
  	return angularNorm * val;
}

// Maple integral of the PDF _without_ time resolution
// int(f_LL1 * exp( -time/tau1 ) + (1 - f_LL1)*exp( -time/tau2 ), time=tmin..tmax);
//                  tmin               tmin               tmin                           tmax               tmax               tmax
// f_LL1 tau1 exp(- ----) + tau2 exp(- ----) - tau2 exp(- ----) f_LL1 - f_LL1 tau1 exp(- ----) - tau2 exp(- ----) + tau2 exp(- ----) f_LL1
//                  tau1               tau2               tau2                           tau1               tau2               tau2
//

double Bs2JpsiPhiLongLivedBkg_withTimeRes::Normalisation(PhaseSpaceBoundary * boundary)
{
        double tmin = 0.;
        double tmax = 0.;
        IConstraint * timeBound = boundary->GetConstraint("time");
        if ( timeBound->GetUnit() == "NameNotFoundError" )
        {
                cerr << "Bound on time not provided" << endl;
                return -1.;
        }
        else
        {
                tmin = timeBound->GetMinimum();
                tmax = timeBound->GetMaximum();
        }
	
        double angularNorm = 1./( 2. * 2. * 2.*TMath::Pi() );

        double tau1 = allParameters.GetPhysicsParameter( tau1Name )->GetValue();
        double tau2 = allParameters.GetPhysicsParameter( tau2Name )->GetValue();
        double f_LL1 = allParameters.GetPhysicsParameter( f_LL1Name )->GetValue();
        double sigmaLL1 = allParameters.GetPhysicsParameter( sigmaLL1Name )->GetValue();
        double sigmaLL2 = allParameters.GetPhysicsParameter( sigmaLL2Name )->GetValue();

	double R1 = erfcInt( tmin, tau1, sigmaLL1);
	double R2 = erfcInt( tmax, tau1, sigmaLL1);
	double R3 = erfcInt( tmin, tau2, sigmaLL2);
	double R4 = erfcInt( tmax, tau2, sigmaLL2);

  	double val = f_LL1 * (R2 - R1) + (1.-f_LL1) * ( R4 - R3 );
 
	return angularNorm * val;
}

// Convolution of a gaussian resolution of width sigma with exp(-time*Gamma) gives:
// R = 1/2 * exp(-time*Gamma + sigma^2*Gamma^2/2)*erfc( -(time - sigma^2*Gamma )/(sqrt(2)*sigma))
double Bs2JpsiPhiLongLivedBkg_withTimeRes::erfc( double time, double tau, double sigma)
{
	double R = 0.5 * exp( -time/tau + sigma*sigma/(2.*tau*tau) )
                       * RooMath::erfc( -( time - sigma*sigma/tau )/(sqrt(2.)*sigma));
	return R;
}

// Mathematica integral of the exp * erf
//Integrate[(1*Exp[-(x/t) + s^2/(2*t^2)]* Erfc[-((x - s^2/t)/(Sqrt[2]*s))])/2, x] ==
//(t*(Erf[x/(Sqrt[2]*s)] - E^((s^2 - 2*t*x)/(2*t^2))* Erfc[(s^2 - t*x)/(Sqrt[2]*s*t)]))/2
double Bs2JpsiPhiLongLivedBkg_withTimeRes::erfcInt( double tlimit, double tau, double sigma)
{
	double val = 0.5 * (tau * ( RooMath::erf( tlimit/(sqrt(2.)*sigma) ) 
				  - exp( (sigma*sigma - 2*tau*tlimit)/(2.*tau*tau) )
				  * RooMath::erfc( (sigma*sigma - tau*tlimit)/(sqrt(2.)*sigma*tau) )
				  )
			   );
	return val;
}
