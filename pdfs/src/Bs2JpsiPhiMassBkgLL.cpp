// $Id: Bs2JpsiPhiMassBkgLL.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhiMassBkgLL Bs2JpsiPhiMassBkgLL.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi mass background
 *
 *  @author Pete Clarke - copy to give two different PDFs for LL and prompt
 *  @date 2010-01-24
 */

#include "Bs2JpsiPhiMassBkgLL.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

#include <float.h>

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

PDF_CREATOR( Bs2JpsiPhiMassBkgLL );

//Constructor
Bs2JpsiPhiMassBkgLL::Bs2JpsiPhiMassBkgLL( PDFConfigurator* configurator ) : 
	// Physics parameters
	  alphaM_llName	( configurator->getName("alphaM_ll") )
        // Observables
        , recoMassName  ( configurator->getName("mass") )
	, constraint_recoMassName( configurator->getName("mass") )
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhiMassBkgLL::MakePrototypes()
{
        allObservables.push_back( recoMassName );
	constraint_recoMassName = recoMassName;

        //Make the parameter set
        vector<string> parameterNames;
        parameterNames.push_back( alphaM_llName );
        allParameters = ParameterSet(parameterNames);
}

//Destructor
Bs2JpsiPhiMassBkgLL::~Bs2JpsiPhiMassBkgLL()
{
}


//Calculate the function value
double Bs2JpsiPhiMassBkgLL::Evaluate(DataPoint * measurement)
{
  	double alphaM_ll = allParameters.GetPhysicsParameter( alphaM_llName )->GetValue();
	
	// Get the observable
        double mass = measurement->GetObservable( recoMassName )->GetValue();
	
	double val = exp( -alphaM_ll * mass);
	
	//To take out a large scale variation  - this is arbitrary provided it is same in Evaluate() and Normalisation()
	double scaleFactor = exp( -alphaM_ll * 5366.0 );
	val /= scaleFactor ;
	
  	return val;
}


double Bs2JpsiPhiMassBkgLL::Normalisation(PhaseSpaceBoundary * boundary)
{
	double mhigh, mlow ;
	
	IConstraint * massBound = boundary->GetConstraint( constraint_recoMassName );
	if ( massBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on mass not provided in Bs2JpsiPhiMassBkgLL" << endl;
		return 1.0 ;
	}
	else
	{
		mlow = massBound->GetMinimum();
		mhigh = massBound->GetMaximum();
	}
	
	double alphaM_ll = allParameters.GetPhysicsParameter( alphaM_llName )->GetValue();
	double integral =0;
	
	if( fabs(alphaM_ll-0.)<DOUBLE_TOLERANCE  ) {
		integral = mhigh-mlow ;   // this was added by PELC to catch a divide by zero Nov-2010
	}
	else {
		integral = (1.0/alphaM_ll)* (exp(-alphaM_ll*mlow) - exp(-alphaM_ll*mhigh)) ;

		//Code prepared to take out a large scale variation - this is arbitrary provided it is same in Evaluate() and Normalisation()
		double scaleFactor = exp( -alphaM_ll * 5366.0 );
		integral /= scaleFactor ;
	}
	return integral;
}
