// $Id: Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist.cpp
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution
 *
 *  @author Pete Clarke
 *  @date 2010-12-12
 */

#include "Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist );

//Constructor
Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist::Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist( PDFConfigurator* configurator ) :

	// Physics parameters
	  tauLL1Name	( configurator->getName("tau_LL1") )
	, tauLL2Name	( configurator->getName("tau_LL2") )
//	, f_JpsiName	( configurator->getName("f_Jpsi") )
	, f_LL1Name      ( configurator->getName("f_LL1") )
        , sigmaLL1Name  ( configurator->getName("sigma_LL1") )
        , sigmaLL2Name  ( configurator->getName("sigma_LL2") )
        , timeResLL1FracName( configurator->getName("timeResLL1Frac") )
	, f_JpsiName    ( configurator->getName("f_Jpsi") )
	, f_NoJpsiName  ( configurator->getName("f_NoJpsi") )

// Observables
        , timeName      ( configurator->getName("time") )
        , cosThetaName  ( configurator->getName("cosTheta") )
        , phiName       ( configurator->getName("phi") )
        , cosPsiName    ( configurator->getName("cosPsi") )

	, tauLL1(), tauLL2(), f_LL1(), sigmaLL(), sigmaLL1(),
	sigmaLL2(), timeResLL1Frac(), f_Jpsi(), f_NoJpsi(),
	tlow(), thigh(), time(), cosTheta(), phi(), cosPsi()

{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( f_LL1Name );
	parameterNames.push_back( tauLL1Name );
	parameterNames.push_back( tauLL2Name );
	parameterNames.push_back( timeResLL1FracName );
	parameterNames.push_back( sigmaLL1Name );
	parameterNames.push_back( sigmaLL2Name );
//	parameterNames.push_back( f_JpsiName );
	parameterNames.push_back( f_NoJpsiName );

	allParameters = ParameterSet(parameterNames);
}

//Destructor
Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist::~Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist()
{
}

bool Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
        bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
        tauLL1      = allParameters.GetPhysicsParameter( tauLL1Name )->GetValue();
        tauLL2      = allParameters.GetPhysicsParameter( tauLL2Name )->GetValue();
        f_LL1       = allParameters.GetPhysicsParameter( f_LL1Name )->GetValue();
        sigmaLL1    = allParameters.GetPhysicsParameter( sigmaLL1Name )->GetValue();
        sigmaLL2    = allParameters.GetPhysicsParameter( sigmaLL2Name )->GetValue();
        timeResLL1Frac = allParameters.GetPhysicsParameter( timeResLL1FracName )->GetValue();
//	f_Jpsi =  allParameters.GetPhysicsParameter( f_JpsiName )->GetValue();
	f_NoJpsi = allParameters.GetPhysicsParameter( f_NoJpsiName )->GetValue();

	return isOK;
}

//Calculate the function value
double Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist::Evaluate(DataPoint * measurement)
{
	// Observable
	time = measurement->GetObservable( timeName )->GetValue();
	cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
	phi      = measurement->GetObservable( phiName )->GetValue();
	cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();

	if( timeResLL1Frac >= 0.9999 )
        {
                // Set the member variable for time resolution to the first value and calculate
                sigmaLL = sigmaLL1;
                return buildPDFnumerator()*angularFactor();
        }
        else
        {
                // Set the member variable for time resolution to the first value and calculate
                sigmaLL = sigmaLL1;
                double val1 = buildPDFnumerator();
                // Set the member variable for time resolution to the second value and calculate
                sigmaLL = sigmaLL2;
                double val2 = buildPDFnumerator();
                return (timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2) * angularFactor();
        }
}

double Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist::buildPDFnumerator()
{
	// Sum of two exponentials, using the time resolution functions

	if( f_LL1 >= 0.9999 ) {
		if( tauLL1 <= 0 ) {
			cout << " In Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist() you gave a negative or zero lifetime for tauLL1 " << endl ;
			exit(1) ;
		}
		double val = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
		return val;
	}
	else {
		if( (tauLL1 <= 0) ||  (tauLL2 <= 0) ) {
			cout << " In Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
			exit(1) ;
		}
		double val1 = Mathematics::Exp(time, 1./tauLL1, sigmaLL);
		double val2 = Mathematics::Exp(time, 1./tauLL2, sigmaLL);
		double val = f_LL1 * val1 + (1. - f_LL1) * val2;
		return val;
	}
}

double Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist::Normalisation(PhaseSpaceBoundary * boundary)
{
        IConstraint * timeBound = boundary->GetConstraint( timeName );
        if ( timeBound->GetUnit() == "NameNotFoundError" )
        {
                cerr << "Bound on time not provided" << endl;
                return -1.;
        }
        else
        {
                tlow = timeBound->GetMinimum();
                thigh = timeBound->GetMaximum();
        }

	if( timeResLL1Frac >= 0.9999 )
        {
                // Set the member variable for time resolution to the first value and calculate
                sigmaLL = sigmaLL1;
                return buildPDFdenominator();
        }
        else
        {
                // Set the member variable for time resolution to the first value and calculate
                sigmaLL = sigmaLL1;
                double val1 = buildPDFdenominator();
                // Set the member variable for time resolution to the second value and calculate
                sigmaLL = sigmaLL2;
                double val2 = buildPDFdenominator();
                return timeResLL1Frac*val1 + (1. - timeResLL1Frac)*val2;
        }
}

double Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist::buildPDFdenominator()
{
	// Sum of two exponentials, using the time resolution functions

	if( f_LL1 >= 0.9999 ) {
		if( tauLL1 <= 0 ) {
			cout << " In Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist() you gave a negative or zero lifetime for tauLL1 " << endl ;
			exit(1) ;
		}
		double val = Mathematics::ExpInt(tlow, thigh, 1./tauLL1, sigmaLL);
		return val;
	}
	else {
		if( (tauLL1 <= 0) ||  (tauLL2 <= 0) ) {
			cout << " In Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist() you gave a negative or zero lifetime for tauLL1/2 " << endl ;
			exit(1) ;
		}
		double val1 = Mathematics::ExpInt(tlow, thigh, 1./tauLL1, sigmaLL);
		double val2 = Mathematics::ExpInt(tlow, thigh, 1./tauLL2, sigmaLL);
		double val = f_LL1 * val1 + (1. - f_LL1) * val2;
		return val;
	}

}


double Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist::angularFactor( )
{

	// Distribution for cosPsi : f_cp

	double p0_cp = 1.0;
	double p2_cp = 0.373 ;
	double p4_cp = -0.216 ;
	double cpsq = cosPsi*cosPsi ;
	double cp4  = cpsq*cpsq ;

	double f_cp		 = p0_cp + (p2_cp*cpsq) + (p4_cp*cp4) ;
	double f_cp_norm = 2.*(p0_cp + (p2_cp/3) + (p4_cp/5) ) ;


	//Second Distribution for cosPsi  : f_cpB
	double p0_cpB = 1.0;
	double p2_cpB = 0.216; //double check

	double f_cpB = p0_cpB + (p2_cpB*cpsq);
	double f_cpB_norm = 2.*(p0_cpB + (p2_cpB/3)); //double check

	// Distribution for phi : f_p

	double p0_p = 1.0 ;
	double p2_p = 0.287 ;
	double p4_p = 0.180 ;
	double pi = TMath::Pi() ;
	double pi3 = pi*pi*pi ;
	double pi5 = pi3*pi*pi ;

	double v = ( TMath::Abs(phi) -(pi/2.) ) ;
	double vsq = v*v ;
	double v4 = vsq*vsq ;
	double f_p		= p0_p + (p2_p*vsq) + (p4_p*v4) ;
	double f_p_norm = 2.*( (p0_p*pi) + (p2_p*pi3/12.) + (p4_p*pi5/80.) ) ;

	//Second Distribution for phi : f_pB

	double p0_pB = 1.0;

	double f_pB = p0_pB;
	double f_pB_norm = 2 * pi; //delete?

	// Distribution for cosTheta : f_ct

	double f1_ct = 0.064 ;
	double f2_ct = 0.624 ;
	double f0_ct = 1 - f1_ct - f2_ct ;
	double sig1 = 0.23 ;
	double sig2 = 1.66 ;
	double ctsq = cosTheta*cosTheta ;

	double bw1 = 1. / ( ctsq + (0.25*sig1*sig1) ) ;
	double bw2 = 1. / ( ctsq + (0.25*sig2*sig2) ) ;
	double bw1_norm = (4.*atan(2./sig1)/sig1) ;
	double bw2_norm = (4.*atan(2./sig2)/sig2) ;


	//double f_ct = f0_ct + f1_ct*bw1 + f2_ct*bw2 ;
	double f_ct = f0_ct/2. + f1_ct*bw1/bw1_norm + f2_ct*bw2/bw2_norm ;
	//double f_ct_norm = 1. ;
	//double f_ct_norm = f0_ct*2. + f1_ct*bw1_norm + f2_ct*bw2_norm ;
	//double f_ct = 1. ;
	double f_ct_norm = 1. ;


	//Second Distribution for cos Theta : f_ctB

	double f_ctB = 1.0;
	double f_ctB_norm = 2;

	// Put it together
	//double factor = (f_cp * f_p * f_ct ) / ( f_cp_norm * f_p_norm * f_ct_norm ) ;

	//double N_a = 0.996;
	//double N_b = 0.004;



	//New factor
	double factor = f_NoJpsi * ((f_cp * f_p * f_ct) / (f_cp_norm * f_p_norm * f_ct_norm) )
			+ (1-f_NoJpsi) * ((f_cpB * f_pB * f_ctB) /( f_cpB_norm * f_pB_norm * f_ctB_norm ));


	// For testing
	//double factor = ( 1. * 1. * 1. ) / ( 2. * 2*pi * 2. ) ;


	return factor ;

}
