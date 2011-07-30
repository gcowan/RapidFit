// $Id: BsMass.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class BsMass BsMass.cpp
 *
 *  RapidFit PDF for Bs mass
 *
 *
 *  @author Pete
 *  @date 2011-07-30
 */

#include "BsMass.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
BsMass::BsMass(PDFConfigurator configurator) :
	// Physics parameters
	  f_sig_m1Name	( configurator.getName("f_sig_m1") )
	, sigma_m1Name	( configurator.getName("sigma_m1") )
	, ratio_21Name	( configurator.getName("ratio_21") )
	, m_BsName		( configurator.getName("m_Bs") )
	// Observables
	, recoMassName	( configurator.getName("mass") )
{
	MakePrototypes();
}

//Make the data point and parameter set
void BsMass::MakePrototypes()
{
	// Observables
	allObservables.push_back( recoMassName );

        //Make the parameter set
        vector<string> parameterNames;
        parameterNames.push_back( f_sig_m1Name );
        parameterNames.push_back( sigma_m1Name );
        parameterNames.push_back( ratio_21Name );
        parameterNames.push_back( m_BsName );
        allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
BsMass::~BsMass()
{
}


//Calculate the function value
double BsMass::Evaluate(DataPoint * measurement)
{
	// Get the physics parameters
  	double f_sig_m1  = allParameters.GetPhysicsParameter( f_sig_m1Name )->GetValue();
  	double sigma_m1 = allParameters.GetPhysicsParameter( sigma_m1Name )->GetValue();
  	double ratio_21 = allParameters.GetPhysicsParameter( ratio_21Name )->GetValue();
  	double m_Bs = allParameters.GetPhysicsParameter( m_BsName )->GetValue();

	// Get the observable
	double mass = measurement->GetObservable( recoMassName )->GetValue();
	
	//Construct second width from ratio
	double sigma_m2 = sigma_m1 * ratio_21 ;

	double returnValue = 0;
	
	if( f_sig_m1 > 0.9999 ) {
		double factor1 = 1./(sigma_m1*sqrt(2.*TMath::Pi()));
		double deltaMsq = ( mass - m_Bs )*( mass - m_Bs );		
		double exp1 = exp( -deltaMsq / ( 2. * sigma_m1 * sigma_m1 ) );
		returnValue = factor1 * exp1 ;		
	}
	else {
		double factor1 = 1./(sigma_m1*sqrt(2.*TMath::Pi()));
		double factor2 = 1./(sigma_m2*sqrt(2.*TMath::Pi()));
		double deltaMsq = ( mass - m_Bs )*( mass - m_Bs );
		double exp1 = exp( -deltaMsq / ( 2. * sigma_m1 * sigma_m1 ) );
		double exp2 = exp( -deltaMsq / ( 2. * sigma_m2 * sigma_m2 ) );
		returnValue = f_sig_m1 * factor1 * exp1 + (1. - f_sig_m1) * factor2 * exp2;
	}

  	return returnValue ;
}


// Normalisation
double BsMass::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	//	Stupid gcc
	(void)measurement;
	(void)boundary;
	// Assumes that the mass integration limits are +/- Infinity
	// So take sufficiently large mass window.
	return 1.0;
}

