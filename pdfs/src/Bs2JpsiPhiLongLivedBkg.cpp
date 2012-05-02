// $Id: Bs2JpsiPhiLongLivedBkg.cpp,v 1.3 2009/11/11 18:30:10 bwynne Exp $
/** @class Bs2JpsiPhiLongLivedBkg Bs2JpsiPhiLongLivedBkg.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi long lived background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-07-30
 */

#include "Bs2JpsiPhiLongLivedBkg.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( Bs2JpsiPhiLongLivedBkg );

//Constructor
Bs2JpsiPhiLongLivedBkg::Bs2JpsiPhiLongLivedBkg( PDFConfigurator* configurator ) :
	// Physics parameters
	tau1Name	( configurator->getName("tau_LL1") )
	, tau2Name	( configurator->getName("tau_LL2") )
	, f_LL1Name	( configurator->getName("f_LL1") )

        // Observables
        , timeName      ( configurator->getName("time") )
	, constraint_timeName()
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhiLongLivedBkg::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	constraint_timeName=timeName;


        //Make the parameter set
        vector<string> parameterNames;
        parameterNames.push_back( tau1Name );
        parameterNames.push_back( tau2Name );
        parameterNames.push_back( f_LL1Name );
        allParameters = ParameterSet(parameterNames);
}

//Destructor
Bs2JpsiPhiLongLivedBkg::~Bs2JpsiPhiLongLivedBkg()
{
}

//Calculate the function value
double Bs2JpsiPhiLongLivedBkg::Evaluate(DataPoint * measurement)
{
	// Observable
        double time = measurement->GetObservable( timeName )->GetValue();

	if ( time < 0.0 ) return 0.0;
  	double tau1 = allParameters.GetPhysicsParameter( tau1Name )->GetValue();
  	double tau2 = allParameters.GetPhysicsParameter( tau2Name )->GetValue();
  	double f_LL1 = allParameters.GetPhysicsParameter( f_LL1Name )->GetValue();
  	double val = f_LL1 * exp( -time/tau1 ) + (1. - f_LL1)*exp( -time/tau2 );

  	return val;
}

// Maple integral
// int(f_LL1 * exp( -time/tau1 ) + (1 - f_LL1)*exp( -time/tau2 ), time=tmin..tmax);
//                  tmin               tmin               tmin                           tmax               tmax               tmax
// f_LL1 tau1 exp(- ----) + tau2 exp(- ----) - tau2 exp(- ----) f_LL1 - f_LL1 tau1 exp(- ----) - tau2 exp(- ----) + tau2 exp(- ----) f_LL1
//                  tau1               tau2               tau2                           tau1               tau2               tau2
double Bs2JpsiPhiLongLivedBkg::Normalisation(PhaseSpaceBoundary * boundary)
{
        double tmin = 0.;
        double tmax = 0.;
        IConstraint * timeBound = boundary->GetConstraint( constraint_timeName );
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

	if ( tmin < 0. ) tmin = 0.;
	if ( tmax < 0. ) tmax = 0.;

        double tau1 = allParameters.GetPhysicsParameter( tau1Name )->GetValue();
        double tau2 = allParameters.GetPhysicsParameter( tau2Name )->GetValue();
        double f_LL1 = allParameters.GetPhysicsParameter( f_LL1Name )->GetValue();

  	double val = f_LL1 * tau1 * ( exp(-tmin/tau1) - exp(-tmax/tau1) )
	      + (1.-f_LL1) * tau2 * ( exp(-tmin/tau2) - exp(-tmax/tau2) );

	return val;
}
