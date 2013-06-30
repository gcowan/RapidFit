// $Id: Bs2JpsiPhiMassBkg.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhiMassBkg Bs2JpsiPhiMassBkg.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi mass background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-07-30
 */

#include "TMath.h"
#include <cmath>

#include "Bs2JpsiPhiMassBkg.h"
#include <iostream>
#include "math.h"
#include <float.h>

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

PDF_CREATOR( Bs2JpsiPhiMassBkg );

//Constructor
Bs2JpsiPhiMassBkg::Bs2JpsiPhiMassBkg(PDFConfigurator* configurator) :
	// Physics parameters
	  alphaM_prName	( configurator->getName("alphaM_pr" ))
        // Observables
        , recoMassName  (configurator->getName( "mass" ))
	, constraint_recoMassName( "mass" )
{
	cout << "Constructing Bs2JpsiPhiMassBkg" << endl;
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhiMassBkg::MakePrototypes()
{
        allObservables.push_back( recoMassName );
	constraint_recoMassName = recoMassName;

        //Make the parameter set
        vector<string> parameterNames;
        parameterNames.push_back( alphaM_prName );
        allParameters = ParameterSet(parameterNames);
}

//Destructor
Bs2JpsiPhiMassBkg::~Bs2JpsiPhiMassBkg()
{
}


//Calculate the function value
double Bs2JpsiPhiMassBkg::Evaluate(DataPoint * measurement)
{
  	double alphaM_pr = allParameters.GetPhysicsParameter( alphaM_prName )->GetValue();

	// Get the observable
        double mass = measurement->GetObservable( recoMassName )->GetValue();

	double val = exp( -alphaM_pr * mass);
	
	//To take out a large scale variation  - this is arbitrary provided it is same in Evaluate() and Normalisation()
	double scaleFactor = exp( -alphaM_pr * 5366.0 );
	val /= scaleFactor ;
	
  	return val;
}


double Bs2JpsiPhiMassBkg::Normalisation(PhaseSpaceBoundary * boundary)
{
	double mhigh, mlow ;

	IConstraint * massBound = boundary->GetConstraint( constraint_recoMassName );
	if ( massBound->GetUnit() == "NameNotFoundError" )
	{
		PDF_THREAD_LOCK
		cerr << "Bound on mass not provided in Bs2JpsiPhiMassBkg" << endl;
		PDF_THREAD_UNLOCK
		return 1.0 ;
	}
	else
	{
		mlow = massBound->GetMinimum();
		mhigh = massBound->GetMaximum();
	}

	double alphaM_pr = allParameters.GetPhysicsParameter( alphaM_prName )->GetValue();
	double integral =0;

	if( fabs( alphaM_pr - 0. ) < DOUBLE_TOLERANCE ) {
		integral = mhigh-mlow ;   // this was added by PELC to catch a divide by zero Nov-2010
	}
	else {
		integral = (1.0/alphaM_pr)* (exp(-alphaM_pr*mlow) - exp(-alphaM_pr*mhigh)) ;

		//Code prepared to take out a large scale variation - this is arbitrary provided it is same in Evaluate() and Normalisation()
		double scaleFactor = exp( -alphaM_pr * 5366.0 );
		integral /= scaleFactor ;
	}

	if( std::isnan(integral) )
	{
		PDF_THREAD_LOCK
		cout << "scale factor: " << exp( -alphaM_pr * 5366.0 ) << endl;
		cout << "alphaM_pr: " << alphaM_pr << endl;
		cout << integral << endl;
		boundary->Print();
		allParameters.Print();
		PDF_THREAD_UNLOCK
	}

	return integral;
}

