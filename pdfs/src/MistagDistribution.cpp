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
	GFshoulderName	( configurator->getName("MistagShoulder" ) ),
	// Observables
	GFxName( configurator->getName("mistag") ),
	tagName( configurator->getName("tag") )
	, gamma(0.), mu(0.), beta(0.)
	, NumericallyIntegrate(false)
{
	NumericallyIntegrate = configurator->isTrue( "NumericallyIntegrate" );
	std::cout << "Constructing PDF: MistagDistribution " << std::endl ;

	MakePrototypes();
}


//Make the data point and parameter set
void MistagDistribution::MakePrototypes()
{
	// Observables
	allObservables.push_back( GFxName );
	allObservables.push_back( tagName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( GFgammaName );
	parameterNames.push_back( GFmuName );
	parameterNames.push_back( GFbetaName );
	parameterNames.push_back( GFshoulderName );

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
  	shoulder  = allParameters.GetPhysicsParameter( GFshoulderName )->GetValue();

	return isOK;
}



//Calculate the function value
double MistagDistribution::Evaluate(DataPoint * measurement)
{
	if( abs(measurement->GetObservable( tagName )->GetValue()) <= std::numeric_limits<double>::epsilon() ) return 1.;
	//cout << measurement->GetObservable( tagName )->GetValue();
	// Get the observable
	x = 0.5 - measurement->GetObservable( GFxName )->GetValue();

	double returnValue ;

	if( x > mu ) returnValue = TMath::GammaDist( x, gamma, mu, beta );
	else if( shoulder > 0.0) returnValue = shoulder ;
	else returnValue = 1E-9 ;

  	return returnValue;
}


// Normalisation
double MistagDistribution::Normalisation( DataPoint* input, PhaseSpaceBoundary * boundary)
{
	if( NumericallyIntegrate ) return -1.;

	(void)boundary;
	if( boundary->GetConstraint( GFxName )->IsDiscrete() )
	{
		return this->Evaluate( input );
	}
	if( abs(input->GetObservable( tagName )->GetValue()) <= std::numeric_limits<double>::epsilon() ) return 1.;
	//return -1.;
	return 1.0 +  shoulder*mu;
}

