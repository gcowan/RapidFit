// $Id: MistagDistribution.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class MistagDistribution MistagDistribution.cpp
 *
 *  RapidFit PDF for Bs mass
 *
 *
 *  @author Pete
 *  @date 2011-07-30
 */

#include "MistagDistribution.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( MistagDistribution );

//Constructor
MistagDistribution::MistagDistribution(PDFConfigurator* configurator) : BasePDF(),
	// Physics parameters
	GFgammaName	( configurator->getName("MistagGamma" ) ),
	GFmuName	( configurator->getName("MistagMu" ) ),
	GFbetaName	( configurator->getName("MistagBeta" ) ),
	// Observables
	GFxName( configurator->getName("mistag") )
	, gamma(0.), mu(0.), beta(0.)
{
	std::cout << "Constructing PDF: MistagDistribution " << std::endl ;

	MakePrototypes();
}


//Make the data point and parameter set
void MistagDistribution::MakePrototypes()
{
	// Observables
	allObservables.push_back( GFxName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( GFgammaName );
	parameterNames.push_back( GFmuName );
	parameterNames.push_back( GFbetaName );

	allParameters = ParameterSet(parameterNames);
}

//Destructor
MistagDistribution::~MistagDistribution()
{
}


bool MistagDistribution::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  	gamma = allParameters.GetPhysicsParameter( GFgammaName )->GetValue();
  	mu    = allParameters.GetPhysicsParameter( GFmuName )->GetValue();
  	beta  = allParameters.GetPhysicsParameter( GFbetaName )->GetValue();
	
	return isOK;
}



//Calculate the function value
double MistagDistribution::Evaluate(DataPoint * measurement)
{
	// Get the observable
	double x = 0.5 - measurement->GetObservable( GFxName )->GetValue();
	
	double returnValue = 0.00001;
	
	if( x > mu ) returnValue = TMath::GammaDist( x, gamma, mu, beta );

  	return returnValue ;
}


// Normalisation
double MistagDistribution::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void)boundary;
	return 1.0;
}

