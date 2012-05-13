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

PDF_CREATOR( Bs2DsPiMassSignal );

//Constructor
Bs2DsPiMassSignal::Bs2DsPiMassSignal( PDFConfigurator* configurator ) : 
	// Physics parameters
	  f_sig_m1Name	( configurator->getName("f_sig_m1") )
	, sigma_m1Name	( configurator->getName("sigma_m1") ) 
	, alpha_m1Name	( configurator->getName("alpha_m1") ) 
	, eta_m1Name	( configurator->getName("eta_m1") ) 
	, sigma_m2Name	( configurator->getName("sigma_m2") ) 
	, alpha_m2Name	( configurator->getName("alpha_m2") ) 
	, eta_m2Name	( configurator->getName("eta_m2") ) 
	, m_BsName	( configurator->getName("m_Bs") )
	// Observables
	, recoMassName	( configurator->getName("mass") )
	, constraint_recoMassName( configurator->getName("mass") )

{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2DsPiMassSignal::MakePrototypes()
{
	// Observables
	allObservables.push_back( recoMassName );
	constraint_recoMassName = recoMassName;

        //Make the parameter set
        vector<string> parameterNames;
        parameterNames.push_back( f_sig_m1Name );
        parameterNames.push_back( sigma_m1Name );
		parameterNames.push_back( alpha_m1Name );
		parameterNames.push_back( eta_m1Name );
        parameterNames.push_back( sigma_m2Name );
		parameterNames.push_back( alpha_m2Name );
		parameterNames.push_back( eta_m2Name );	
        parameterNames.push_back( m_BsName );
        allParameters = ParameterSet(parameterNames);
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
  	double alpha_m1 = allParameters.GetPhysicsParameter( alpha_m1Name )->GetValue();
  	double eta_m1 = allParameters.GetPhysicsParameter( eta_m1Name )->GetValue();  	
	double sigma_m2 = allParameters.GetPhysicsParameter( sigma_m2Name )->GetValue();
  	double alpha_m2 = allParameters.GetPhysicsParameter( alpha_m2Name )->GetValue();
  	double eta_m2 = allParameters.GetPhysicsParameter( eta_m2Name )->GetValue();  	
  	double m_Bs = allParameters.GetPhysicsParameter( m_BsName )->GetValue();
  	
	// Get the observable
        double mass = measurement->GetObservable( recoMassName )->GetValue();


	double deltaMs = ( mass - m_Bs );
	double condition1 = deltaMs / sigma_m1;
	double condition2 = deltaMs / sigma_m2;

	double CB_1;
	double CB_2;
	
	if (alpha_m1 < 0) condition1 = -condition1; 
	if (alpha_m2 < 0) condition2 = -condition2; 
	
	Double_t absAlpha_m1 = fabs((Double_t)alpha_m1); 
	Double_t absAlpha_m2 = fabs((Double_t)alpha_m2); 
	
	if (condition1 >= -absAlpha_m1) { 
		CB_1 = exp(-0.5*condition1*condition1); 
	} 
	else { 
		Double_t a = TMath::Power(eta_m1/absAlpha_m1,eta_m1)*exp(-0.5*absAlpha_m1*absAlpha_m1); 
		Double_t b= eta_m1/absAlpha_m1 - absAlpha_m1; 
		
		CB_1 = a/TMath::Power(b - condition1, eta_m1); 
	} 

	if (condition2 >= -absAlpha_m2) { 
		CB_2 = exp(-0.5*condition2*condition2); 
	} 
	else { 
		Double_t a2 = TMath::Power(eta_m2/absAlpha_m2,eta_m2)*exp(-0.5*absAlpha_m2*absAlpha_m2); 
		Double_t b2= eta_m2/absAlpha_m2 - absAlpha_m2; 
		
		CB_2 = a2/TMath::Power(b2 - condition2, eta_m2); 
	} 	
	
	
	double val = f_sig_m1 *  CB_1 + (1. - f_sig_m1) * CB_2;

	return val;
	
 }

double Bs2DsPiMassSignal::Normalisation(PhaseSpaceBoundary * boundary)
{
	(void)boundary;
/*
	double mhigh, mlow ;

	IConstraint * massBound = boundary->GetConstraint( constraint_recoMassName );
	if ( massBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on mass not provided in Bs2DsPiMassSignal" << endl;
		return 1.0 ;
	}
	else
	{
		mlow = massBound->GetMinimum();
		mhigh = massBound->GetMaximum();
	}
	
	// Get the physics parameters
  	double f_sig_m1  = allParameters.GetPhysicsParameter( f_sig_m1Name )->GetValue();
  	double sigma_m1 = allParameters.GetPhysicsParameter( sigma_m1Name )->GetValue();
  	double alpha_m1 = allParameters.GetPhysicsParameter( alpha_m1Name )->GetValue();
  	double eta_m1 = allParameters.GetPhysicsParameter( eta_m1Name )->GetValue();  	
	double sigma_m2 = allParameters.GetPhysicsParameter( sigma_m2Name )->GetValue();
  	double alpha_m2 = allParameters.GetPhysicsParameter( alpha_m2Name )->GetValue();
  	double eta_m2 = allParameters.GetPhysicsParameter( eta_m2Name )->GetValue();	
  	double m_Bs = allParameters.GetPhysicsParameter( m_BsName )->GetValue();
	
	//Reality checks
	if( sigma_m1 <= 0. ) {cout << "Bs2DsPiMassSignal::Normalisation() : sigma_m1 < 0 : " << sigma_m1 << endl ;exit(1); }
	if( sigma_m2 <= 0. ) {cout << "Bs2DsPiMassSignal::Normalisation() : sigma_m2 < 0 : " << sigma_m2 << endl ; exit(1); }	
	
*/
/*	
	double integral1 = 0.;
	double integral2 = 0.;
	
	double m1 = eta_m1 - 1.; //minus?!
	double m2 = eta_m2 - 1.;

	double c1 = m_Bs - alpha_m1*sigma_m1;
	double c2 = m_Bs - alpha_m2*sigma_m2;
	 
	   if( mlow > c1 ) {
	     // Integral from Gaussian part
		   integral1 = 0.5 * ( erf( mhigh/(sqrt(2.)*sigma_m1) ) - erf( mlow/(sqrt(2.)*sigma_m1 )) );
	   } else if( mhigh < c1 ) {
	     // Integral from powerlaw part
	     double k = eta_m1*sigma_m1*exp(-0.5*alpha_m1*alpha_m1)/m1*alpha_m1;
	     if( eta_m1 == 1. ) {
	       double B = eta_m1/alpha_m1 - alpha_m1;
	       integral1 = k/pow(1-alpha_m1*alpha_m1/eta_m1,m1)*log( (B+(m_Bs-mlow)/sigma_m1)/(B+(m_Bs-mhigh)/sigma_m1) );
	     } else {
	       integral1 = k*( pow( 1 + alpha_m1/eta_m1*( (m_Bs-mhigh)/sigma_m1 - alpha_m1 ),-m1 )
	                -pow( 1 + alpha_m1/eta_m1*( (m_Bs-mlow)/sigma_m1 - alpha_m1 ),-m1 ) );
	     }
	   } else {
	 
	     // Integral from both parts
	     double k = eta_m1*sigma_m1*exp(-0.5*alpha_m1*alpha_m1)/m1*alpha_m1;
	     integral1 = sigma_m1*sqrt(M_PI/2.)*( erf(alpha_m1/sqrt(2)) - erf((m_Bs-mhigh)/sqrt(2)*sigma_m1) );
	     if( eta_m1 == 1. ) {
	       double B = eta_m1/alpha_m1 - alpha_m1;
	       integral1 += k/pow(1-alpha_m1*alpha_m1/eta_m1,m1)*log( (B+(m_Bs-mlow)/sigma_m1)/(B+alpha_m1) );
	     } else {
	       integral1 += k*( 1. - pow( 1 + alpha_m1/eta_m1*( (m_Bs-mlow)/sigma_m1 - alpha_m1 ),-m1 ) );
	     }
	   }
*/	
	
	/*
	if( mlow > c2 ) {
		// Integral from Gaussian part
		integral1 = 0.5 * ( erf( mhigh/(sqrt(2.)*sigma_m1) ) - erf( mlow/(sqrt(2.)*sigma_m1 )) );
		//integral2 = sigma_m2*sqrt(M_PI/2.)*( erf((m_Bs-mlow)/sqrt(2)*sigma_m2) - erf((m_Bs-mhigh)/sqrt(2)*sigma_m2) );
	} else if( mhigh < c2 ) {
		// Integral from powerlaw part
		double k = eta_m2*sigma_m2*exp(-0.5*alpha_m2*alpha_m2)/m2*alpha_m2;
		if( eta_m2 == 1. ) {
			double B = eta_m2/alpha_m2 - alpha_m2;
			integral2 = k/pow(1-alpha_m2*alpha_m2/eta_m2,m2)*log( (B+(m_Bs-mlow)/sigma_m2)/(B+(m_Bs-mhigh)/sigma_m2) );
		} else {
			integral2 = k*( pow( 1 + alpha_m2/eta_m2*( (m_Bs-mhigh)/sigma_m2 - alpha_m2 ),-m2 )
	                -pow( 1 + alpha_m2/eta_m2*( (m_Bs-mlow)/sigma_m2 - alpha_m2 ),-m2 ) );
		}
	} else {
		// Integral from both parts
		double k = eta_m2*sigma_m2*exp(-0.5*alpha_m2*alpha_m2)/m2*alpha_m2;
		integral2 = sigma_m2*sqrt(M_PI/2.)*( erf(alpha_m2/sqrt(2)) - erf((m_Bs-mhigh)/sqrt(2)*sigma_m2) );
		if( eta_m2 == 1. ) {
			double B = eta_m2/alpha_m2 - alpha_m2;
			integral2 += k/pow(1-alpha_m2*alpha_m2/eta_m2,m2)*log( (B+(m_Bs-mlow)/sigma_m2)/(B+alpha_m2) );
		} else {
			integral2 += k*( 1. - pow( 1 + alpha_m2/eta_m2*( (m_Bs-mlow)/sigma_m2 - alpha_m2 ),-m2 ) );
		}
	}*/
	//double Integral = f_sig_m1 *  integral1 + (1. - f_sig_m1) * integral2;	
	return -1;
}

