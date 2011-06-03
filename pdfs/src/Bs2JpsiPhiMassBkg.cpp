// $Id: Bs2JpsiPhiMassBkg.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhiMassBkg Bs2JpsiPhiMassBkg.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi mass background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-07-30
 */

#include "Bs2JpsiPhiMassBkg.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

#define DOUBLE_TOLERANCE 1E-6

//Constructor
Bs2JpsiPhiMassBkg::Bs2JpsiPhiMassBkg() :
	// Physics parameters
	  alphaM_prName	( make_pair("alphaM_pr",-1) )
        // Observables
        , recoMassName  ( make_pair("mass",-1) )
	, constraint_recoMassName()
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhiMassBkg::MakePrototypes()
{
        allObservables.push_back( recoMassName.first );
	constraint_recoMassName = recoMassName;

        //Make the parameter set
        vector<string> parameterNames;
        parameterNames.push_back( alphaM_prName.first );
        allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2JpsiPhiMassBkg::~Bs2JpsiPhiMassBkg()
{
}


//Calculate the function value
double Bs2JpsiPhiMassBkg::Evaluate(DataPoint * measurement)
{
  	double alphaM_pr = allParameters.GetPhysicsParameter( &alphaM_prName )->GetValue();

	// Get the observable
        double mass = measurement->GetObservable( &recoMassName )->GetValue();

	double val = exp( -alphaM_pr * mass);
	
	//To take out a large scale variation  - this is arbitrary provided it is same in Evaluate() and Normalisation()
	double scaleFactor = exp( -alphaM_pr * 5366.0 );
	val /= scaleFactor ;
	
  	return val;
}


double Bs2JpsiPhiMassBkg::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	//	Stupid gcc
	(void)measurement;
	double mhigh, mlow ;

	IConstraint * massBound = boundary->GetConstraint( &constraint_recoMassName );
	if ( massBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on mass not provided in Bs2JpsiPhiMassBkg" << endl;
		return 1.0 ;
	}
	else
	{
		mlow = massBound->GetMinimum();
		mhigh = massBound->GetMaximum();
	}

	double alphaM_pr = allParameters.GetPhysicsParameter( &alphaM_prName )->GetValue();
	double integral ;

	if( fabs( alphaM_pr - 0. ) < DOUBLE_TOLERANCE ) {
		integral = mhigh-mlow ;   // this was added by PELC to catch a divide by zero Nov-2010
	}
	else {
		integral = (1.0/alphaM_pr)* (exp(-alphaM_pr*mlow) - exp(-alphaM_pr*mhigh)) ;

		//Code prepared to take out a large scale variation - this is arbitrary provided it is same in Evaluate() and Normalisation()
		double scaleFactor = exp( -alphaM_pr * 5366.0 );
		integral /= scaleFactor ;
		
	}
	return integral;

}
