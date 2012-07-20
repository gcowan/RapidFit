// $Id: LogNormalDistribution.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class LogNormalDistribution LogNormalDistribution.cpp
 *
 *  RapidFit PDF for Bs mass
 *
 *
 *  @author Pete
 *  @date 2011-07-30
 */

#include "LogNormalDistribution.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( LogNormalDistribution );

//Constructor
LogNormalDistribution::LogNormalDistribution(PDFConfigurator* configurator) :
	// Physics parameters
	LNsigmaName	( configurator->getName("LNsigma" ) ),
	LNthetaName	( configurator->getName("LNtheta" ) ),
	LNmName	( configurator->getName("LNm" ) ),
	// Observables
	LNxName( configurator->getName("LNx") )
{
	
	std::cout << "Constructing PDF: LogNormalDistribution " << std::endl ;

	MakePrototypes();
}

//Make the data point and parameter set
void LogNormalDistribution::MakePrototypes()
{
	// Observables
	allObservables.push_back( LNxName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( LNsigmaName );
	parameterNames.push_back( LNthetaName );
	parameterNames.push_back( LNmName );
	allParameters = ParameterSet(parameterNames);
}

//Destructor
LogNormalDistribution::~LogNormalDistribution()
{
}


bool LogNormalDistribution::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  	sigma = allParameters.GetPhysicsParameter( LNsigmaName )->GetValue();
  	theta = allParameters.GetPhysicsParameter( LNthetaName )->GetValue();
  	m     = allParameters.GetPhysicsParameter( LNmName )->GetValue();
	
	return isOK;
}



//Calculate the function value
double LogNormalDistribution::Evaluate(DataPoint * measurement)
{
	// Get the observable
	double x = measurement->GetObservable( LNxName )->GetValue();
	
	double returnValue = 0.000000000001;

	if ( x > theta ) returnValue = TMath::LogNormal( x, sigma, theta, m );

  	return returnValue;
}


// Normalisation
double LogNormalDistribution::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void)boundary;
	// Assumes that the mass integration limits are +/- Infinity
	// So take sufficiently large mass window.
	return 1.0;
}

