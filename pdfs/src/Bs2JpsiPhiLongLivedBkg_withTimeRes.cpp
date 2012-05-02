// $Id: Bs2JpsiPhiLongLivedBkg_withTimeRes.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Bs2JpsiPhiLongLivedBkg_withTimeRes Bs2JpsiPhiLongLivedBkg_withTimeRes.cpp
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#include "Bs2JpsiPhiLongLivedBkg_withTimeRes.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( Bs2JpsiPhiLongLivedBkg_withTimeRes );

//Constructor
Bs2JpsiPhiLongLivedBkg_withTimeRes::Bs2JpsiPhiLongLivedBkg_withTimeRes( PDFConfigurator* configurator ) :

	// Physics parameters
	  tauLL1Name	( configurator->getName("tau_LL1") )
	, tauLL2Name	( configurator->getName("tau_LL2") )
	, f_LL1Name     ( configurator->getName("f_LL1") )
	, sigmaLL1Name	( configurator->getName("sigma_LL1") )
	, sigmaLL2Name	( configurator->getName("sigma_LL2") )
        , timeResLL1FracName( configurator->getName("timeResLL1Frac") )
	// Observables
        , timeName      ( configurator->getName("time") )

	, tauLL1(), tauLL2(), f_LL1(), sigmaLL(), sigmaLL1(),
	sigmaLL2(), timeResLL1Frac(), tlow(), thigh(), time()
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
	parameterNames.push_back( f_LL1Name );
        parameterNames.push_back( tauLL1Name );
        parameterNames.push_back( tauLL2Name );
	parameterNames.push_back( timeResLL1FracName );
        parameterNames.push_back( sigmaLL1Name );
        parameterNames.push_back( sigmaLL2Name );
        allParameters = ParameterSet(parameterNames);
}

//Destructor
Bs2JpsiPhiLongLivedBkg_withTimeRes::~Bs2JpsiPhiLongLivedBkg_withTimeRes()
{
}

bool Bs2JpsiPhiLongLivedBkg_withTimeRes::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
        bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
        tauLL1      = allParameters.GetPhysicsParameter( tauLL1Name )->GetValue();
        tauLL2      = allParameters.GetPhysicsParameter( tauLL2Name )->GetValue();
        f_LL1       = allParameters.GetPhysicsParameter( f_LL1Name )->GetValue();
        sigmaLL1    = allParameters.GetPhysicsParameter( sigmaLL1Name )->GetValue();
        sigmaLL2    = allParameters.GetPhysicsParameter( sigmaLL2Name )->GetValue();
        timeResLL1Frac = allParameters.GetPhysicsParameter( timeResLL1FracName )->GetValue();
	return isOK;
}

//Calculate the function value
double Bs2JpsiPhiLongLivedBkg_withTimeRes::Evaluate(DataPoint * measurement)
{
	// Observable
        time = measurement->GetObservable( timeName )->GetValue();

	if( timeResLL1Frac >= 0.9999 )
        {
                // Set the member variable for time resolution to the first value and calculate
                sigmaLL = sigmaLL1;
                return buildPDFnumerator();
        }
        else
        {
                // Set the member variable for time resolution to the first value and calculate
                sigmaLL = sigmaLL1;
                double val1 = buildPDFnumerator();
                // Set the member variable for time resolution to the second value and calculate
                sigmaLL = sigmaLL2;
                double val2 = buildPDFnumerator();
                return timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2;
        }
}

double Bs2JpsiPhiLongLivedBkg_withTimeRes::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions

	if( f_LL1 >= 0.9999 ) {
		if( tauLL1 <= 0 ) {
			cout << " In Bs2JpsiPhiLongLivedBkg_withTimeRes() you gave a negative or zero lifetime for tauLL1 " << endl ;
			throw(10) ;
		}
		double val = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
		return val;
	}
	else {
		if( (tauLL1 <= 0) ||  (tauLL2 <= 0) ) {
			cout << " In Bs2JpsiPhiLongLivedBkg_withTimeRes() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
			throw(10) ;
		}
		double val1 = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
		double val2 = Mathematics::Exp(time, 1./tauLL2, sigmaLL);
		double val = f_LL1 * val1 + (1. - f_LL1) * val2;
		return val;
	}
}

double Bs2JpsiPhiLongLivedBkg_withTimeRes::Normalisation(PhaseSpaceBoundary * boundary)
{
	IConstraint * timeBound = boundary->GetConstraint( timeName );
        if ( timeBound->GetUnit() == "NameNotFoundError" )
        {
                cerr << "Bound on time not provided" << endl;
                return -1.;
        }
        else
        {
                tlow = timeBound->GetMinimum();
                thigh = timeBound->GetMaximum();
        }

	if( timeResLL1Frac >= 0.9999 )
        {
                // Set the member variable for time resolution to the first value and calculate
                sigmaLL = sigmaLL1;
                return buildPDFdenominator();
        }
        else
        {
                // Set the member variable for time resolution to the first value and calculate
                sigmaLL = sigmaLL1;
                double val1 = buildPDFdenominator();
                // Set the member variable for time resolution to the second value and calculate
                sigmaLL = sigmaLL2;
                double val2 = buildPDFdenominator();
                return timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2;
        }
}

double Bs2JpsiPhiLongLivedBkg_withTimeRes::buildPDFdenominator()
{
	// Sum of two exponentials, using the time resolution functions

	if( f_LL1 >= 0.9999 ) {
		if( tauLL1 <= 0 ) {
			cout << " In Bs2JpsiPhiLongLivedBkg_withTimeRes() you gave a negative or zero lifetime for tauLL1 " << endl ;
			throw(10) ;
		}
		double val = Mathematics::ExpInt(tlow, thigh, 1./tauLL1, sigmaLL);
		return val;
	}
	else {
		if( (tauLL1 <= 0) ||  (tauLL2 <= 0) ) {
			cout << " In Bs2JpsiPhiLongLivedBkg_withTimeRes() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
			throw(10) ;
		}
		double val1 = Mathematics::ExpInt(tlow, thigh, 1./tauLL1, sigmaLL);
		double val2 = Mathematics::ExpInt(tlow, thigh, 1./tauLL2, sigmaLL);
		double val = f_LL1 * val1 + (1. - f_LL1) * val2;
		return val;
	}

}


