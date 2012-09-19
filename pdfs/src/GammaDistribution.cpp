// $Id: GammaDistribution.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class GammaDistribution GammaDistribution.cpp
 *
 *  RapidFit PDF for Bs mass
 *
 *
 *  @author Pete
 *  @date 2011-07-30
 */

#include "GammaDistribution.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( GammaDistribution );

//Constructor
GammaDistribution::GammaDistribution(PDFConfigurator* configurator) :
	// Physics parameters
	GFgammaName	( configurator->getName("GFgamma" ) ),
	GFmuName	( configurator->getName("GFmu" ) ),
	GFbetaName	( configurator->getName("GFbeta" ) ),
	// Observables
	GFxName( configurator->getName("GFx") )
{
	
	std::cout << "Constructing PDF: GammaDistribution " << std::endl ;

	//PDF = new TF1("sigmatPDF","[3]*TMath::GammaDist(x, [0], [2], [1] )",0, 0.1); 
	//PDF->SetParameters(10.4,0.0027, 0.003,1);

	MakePrototypes();
}

//Make the data point and parameter set
void GammaDistribution::MakePrototypes()
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
GammaDistribution::~GammaDistribution()
{
}


bool GammaDistribution::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
  	gamma = allParameters.GetPhysicsParameter( GFgammaName )->GetValue();
  	mu    = allParameters.GetPhysicsParameter( GFmuName )->GetValue();
  	beta  = allParameters.GetPhysicsParameter( GFbetaName )->GetValue();
	
	return isOK;
}



//Calculate the function value
double GammaDistribution::Evaluate(DataPoint * measurement)
{
	//PDF->SetParameters(gamma,beta,mu,1.);
	//PDF->SetParameters(10.4,0.0027, 0.003,1);
	
	// Get the observable
	x = measurement->GetObservable( GFxName )->GetValue();
	
	double returnValue = 0.000000000001;

	if( x > mu ) returnValue = TMath::GammaDist( x, gamma, mu, beta );

  	return returnValue;
}


// Normalisation
double GammaDistribution::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void)boundary;
	// Assumes that the mass integration limits are +/- Infinity
	// So take sufficiently large mass window.
	return 1.0;
}

