// $Id: RaPDF_Bs2JpsiPhiPromptBkg.cpp,v 1.2 2009/11/13 09:57:06 gcowan Exp $
/** @class RaPDF_Bs2JpsiPhiPromptBkg RaPDF_Bs2JpsiPhiPromptBkg.cpp
 *
 *  PDF for Bs2JpsiPhi prompt background
 *
 *  @author Greig A Cowan greig.cowan@cern.ch
 *  @date 2009-11-12
 */

#include "RaPDF_Bs2JpsiPhiPromptBkg.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
RaPDF_Bs2JpsiPhiPromptBkg::RaPDF_Bs2JpsiPhiPromptBkg() : 
	// Observables
	  timeName	( "time" )
	
	// Physics parameters
	, meanName	( "mean" )
	, sigmaName	( "sigma" )
{
	MakePrototypes();
}

//Make the data point and parameter set
void RaPDF_Bs2JpsiPhiPromptBkg::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );

        //Make the parameter set
        vector<string> parameterNames;
        parameterNames.push_back( meanName );
        parameterNames.push_back( sigmaName );
        allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
RaPDF_Bs2JpsiPhiPromptBkg::~RaPDF_Bs2JpsiPhiPromptBkg()
{
}


//Calculate the function value
double RaPDF_Bs2JpsiPhiPromptBkg::Evaluate(DataPoint * measurement)
{
	// Observable
        double time = measurement->GetObservable( timeName )->GetValue();

	if ( time < 0 ) return 0.0;
  	double tau = 1.;
  	double val = 1.;

	double mean= allParameters.GetPhysicsParameter( meanName )->GetValue();
        double sigma = allParameters.GetPhysicsParameter( sigmaName )->GetValue();

        return val;
}


double RaPDF_Bs2JpsiPhiPromptBkg::Normalisation(PhaseSpaceBoundary * boundary)
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

  	double val = 1.;

	return val;
}
