/** @class Bs2DsPiMassSignal Bs2DsPiMassSignal.cpp
 *
 *  RapidFit PDF for Bs2DsPi mass signal
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-07-30
 */

#include "Bs2DsPiMassSignal.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
Bs2DsPiMassSignal::Bs2DsPiMassSignal() : 
	// Physics parameters
	  f_sig_m1Name	( "f_sig_m1" )
	, sigma_m1Name	( "sigma_m1" ) 
	, sigma_m2Name	( "sigma_m2" ) 
	, m_BsName	( "m_Bs" )
	// Observables
	, recoMassName	( "mass")
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2DsPiMassSignal::MakePrototypes()
{
	// Observables
	allObservables.push_back( recoMassName );

        //Make the parameter set
        vector<string> parameterNames;
        parameterNames.push_back( f_sig_m1Name );
        parameterNames.push_back( sigma_m1Name );
        parameterNames.push_back( sigma_m2Name );
        parameterNames.push_back( m_BsName );
        allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2DsPiMassSignal::~Bs2DsPiMassSignal()
{
}


//Calculate the function value
double Bs2DsPiMassSignal::Evaluate(DataPoint * measurement)
{
	// Get the physics parameters
  	double f_sig_m1  = allParameters.GetPhysicsParameter( f_sig_m1Name )->GetValue();
  	double sigma_m1 = allParameters.GetPhysicsParameter( sigma_m1Name )->GetValue();
  	double sigma_m2 = allParameters.GetPhysicsParameter( sigma_m2Name )->GetValue();
  	double m_Bs = allParameters.GetPhysicsParameter( m_BsName )->GetValue();
  	
	// Get the observable
        double mass = measurement->GetObservable( recoMassName )->GetValue();

	double factor1 = 1./(sigma_m1*sqrt(2.*TMath::Pi()));
	double factor2 = 1./(sigma_m2*sqrt(2.*TMath::Pi()));
	double deltaMsq = ( mass - m_Bs )*( mass - m_Bs );

	double exp1 = exp( -deltaMsq / ( 2. * sigma_m1 * sigma_m1 ) );
	double exp2 = exp( -deltaMsq / ( 2. * sigma_m2 * sigma_m2 ) );

	double val = f_sig_m1 * factor1 * exp1 + (1. - f_sig_m1) * factor2 * exp2;
	
  	return val;
}

double Bs2DsPiMassSignal::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	return 1.0;
}
