// $Id: Bs2JpsiPhi.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhi Bs2JpsiPhi.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-07-30
 */

#include "Bs2JpsiPhi.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
Bs2JpsiPhi::Bs2JpsiPhi() : 
	  gammaName     ( "gamma" )
	, deltaGammaName( "deltaGamma" )
	, deltaMName    ( "deltaM")
	, Phi_sName     ( "Phi_s")
	, Azero_sqName  ( "Azero_sq" )
	//, Apara_sqName  ( "Apara_sq")
	, Aperp_sqName  ( "Aperp_sq" )
	, delta_zeroName( "delta_zero" )
	, delta_paraName( "delta_para" )
	, delta_perpName( "delta_perp" )
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhi::MakePrototypes()
{
	//Make the DataPoint prototype
	timeName     = "time";
	cosThetaName = "cosTheta";
	phiName      = "phi";
	cosPsiName   = "cosPsi";
	tagName      = "tag";
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );
	allObservables.push_back( tagName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( deltaGammaName );
	parameterNames.push_back( Aperp_sqName );
	parameterNames.push_back( Azero_sqName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_zeroName );
	parameterNames.push_back( deltaMName );
	parameterNames.push_back( Phi_sName );
	//parameterNames.push_back( Apara_sqName );
	parameterNames.push_back( "tagFraction" );
	parameterNames.push_back( "resolution" );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2JpsiPhi::~Bs2JpsiPhi()
{
}

//Calculate the function value
double Bs2JpsiPhi::Evaluate(DataPoint * measurement)
{
	// The angular functions f1->f6 as defined in roadmap Table 1.
	double f1, f2, f3, f4, f5, f6;
	getAngularFunctions( f1, f2, f3, f4, f5, f6, measurement );	

	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59
	double AzeroAzero, AparaApara, AperpAperp;
	double ImAparaAperp, ReAzeroApara, ImAzeroAperp;
	getTimeDependentAmplitudes( AzeroAzero, AparaApara, AperpAperp
				  , ImAparaAperp, ReAzeroApara, ImAzeroAperp
				  , measurement);

	// What about overall normalisation?
	double eval = f1 * AzeroAzero
		    + f2 * AparaApara
		    + f3 * AperpAperp
		    + f4 * ImAparaAperp
		    + f5 * ReAzeroApara
		    + f6 * ImAzeroAperp;
	return 0.5*eval;
}

void Bs2JpsiPhi::getAngularFunctions( double & f1
					, double & f2
					, double & f3
					, double & f4
					, double & f5
					, double & f6
					, DataPoint * measurement)
{
        // Observables (the stuff your experiment measures)
	double cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
        double phi      = measurement->GetObservable( phiName )->GetValue();
        double cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();

	double sinTheta  = 1. - cosTheta*cosTheta;
	double sinPsi    = 1. - cosPsi*cosPsi;

	double cosPhi    = cos(phi);
	double sinPhi	 = sin(phi);

	double sin2Theta = 2.*sinTheta*cosTheta;
	double sin2Psi	 = 2.*sinPsi*cosPsi;
	double sin2Phi	 = 2.*sinPhi*cosPhi;
	 
	double norm = 9./32./TMath::Pi();
	f1 = 2. * cosPsi*cosPsi * ( 1. - sinTheta*sinTheta * cosPhi*cosPhi ) / norm;
	f2 =      sinPsi*sinPsi * ( 1. - sinTheta*sinTheta * sinPhi*sinPhi ) / norm;
	f3 =      sinPsi*sinPsi * sinTheta*sinTheta / norm;
	f4 = -1.* sinPsi*sinPsi * sin2Theta * sinPhi / norm;
	f5 = sin2Psi * sinTheta*sinTheta * sin2Phi/sqrt(2.) / norm;
	f6 = sin2Psi * sin2Theta * cosPhi/sqrt(2.) / norm;
	return;
}

void Bs2JpsiPhi::getTimeDependentAmplitudes(  double & AzeroAzero
						, double & AparaApara
						, double & AperpAperp
						, double & ImAparaAperp
						, double & ReAzeroApara
						, double & ImAzeroAperp
						, DataPoint * measurement )
{
	// Observable
	double time = measurement->GetObservable( timeName )->GetValue();
	int tag  = (int)measurement->GetObservable( tagName )->GetValue();
        
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	double gamma, deltaGamma, deltaM, Phi_s;
        double Azero_sq, Apara_sq, Aperp_sq;
        double delta_zero, delta_para, delta_perp;
        getPhysicsParameters( gamma, deltaGamma, deltaM, Phi_s, Azero_sq, Apara_sq, Aperp_sq, delta_zero, delta_para, delta_perp);

	// Calculate the quantities that will be used in the returned amplitudes
	double Azero = sqrt( Azero_sq );
	double Apara = sqrt( Apara_sq );
	double Aperp = sqrt( Aperp_sq );
	
	double sinDeltaPerpPara = sin( delta_perp - delta_para );
	double cosDeltaPerpPara = cos( delta_perp - delta_para );
	double sinDeltaPerp     = sin( delta_perp );
	double cosDeltaPerp     = cos( delta_perp );
	double cosDeltaPara	= cos( delta_para );
	
	double expGT = exp( -gamma*time );
	
	double coshDeltaGammaT = cosh( deltaGamma*time/2.);
	double sinhDeltaGammaT = sinh( deltaGamma*time/2.);
	
	double sinPhis = sin( Phi_s );
	double cosPhis = cos( Phi_s );
	
	double sinDeltaMT = sin( deltaM*time );
	double cosDeltaMT = cos( deltaM*time );
	
	// Now calculate the amplitudes
	AzeroAzero = Azero_sq * expGT * ( coshDeltaGammaT - cosPhis * sinhDeltaGammaT + tag * sinPhis * sinDeltaMT ); 
	AparaApara = Apara_sq * expGT * ( coshDeltaGammaT - cosPhis * sinhDeltaGammaT + tag * sinPhis * sinDeltaMT ); 
	AperpAperp = Aperp_sq * expGT * ( coshDeltaGammaT + cosPhis * sinhDeltaGammaT - tag * sinPhis * sinDeltaMT ); 
	
	ImAparaAperp = Apara*Aperp * expGT * ( - cosDeltaPerpPara * sinPhis * sinhDeltaGammaT 
					       + tag * sinDeltaPerpPara * cosDeltaMT
					       - tag * cosDeltaPerpPara * cosPhis * sinDeltaMT );

	ReAzeroApara = Azero*Apara * expGT * cosDeltaPara * ( coshDeltaGammaT - cosPhis * sinhDeltaGammaT
							    + tag * sinPhis * sinDeltaMT );

	ImAzeroAperp = Azero*Aperp * expGT * ( - cosDeltaPerp * sinPhis * sinhDeltaGammaT
                                               + tag * sinDeltaPerp * cosDeltaMT 
                                               - tag * cosDeltaPerp * cosPhis * sinDeltaMT );
	return;
}

// Start calculating the integrals of the time dependent parts of the amplitudes
void Bs2JpsiPhi::getTimeAmplitudeIntegrals( double & AzeroAzeroInt
                                                 , double & AparaAparaInt
                                                 , double & AperpAperpInt
                                                 , double & ImAparaAperpInt
                                                 , double & ReAzeroAparaInt
                                                 , double & ImAzeroAperpInt
				       	         , DataPoint * measurement
				       	         , PhaseSpaceBoundary * boundary)
{
        // Observable
        int tag = (int)measurement->GetObservable( tagName )->GetValue();	

	// Get the latest values of the physics parameters that will be used in the evaluation of the integral
	double gamma, deltaGamma, deltaM, Phi_s;
	double Azero_sq, Apara_sq, Aperp_sq;
	double delta_zero, delta_para, delta_perp;
	getPhysicsParameters( gamma, deltaGamma, deltaM, Phi_s, Azero_sq, Apara_sq, Aperp_sq, delta_zero, delta_para, delta_perp);

	// Need these when working out the time integral of the amplitudes
	double Azero = sqrt( Azero_sq );
	double Apara = sqrt( Apara_sq );
	double Aperp = sqrt( Aperp_sq );

	// Get the bounds of the integration
	double timeMin = boundary->GetConstraint( timeName )->GetMinimum();
	double timeMax = boundary->GetConstraint( timeName )->GetMaximum();

	AzeroAzeroInt = getAzeroAzeroTimeInt( Azero_sq, gamma, deltaGamma, deltaM, Phi_s, timeMin, timeMax, tag );	
	AparaAparaInt = getAparaAparaTimeInt( Apara_sq, gamma, deltaGamma, deltaM, Phi_s, timeMin, timeMax, tag );	
	AperpAperpInt = getAperpAperpTimeInt( Aperp_sq, gamma, deltaGamma, deltaM, Phi_s, timeMin, timeMax, tag );	
	
	ImAparaAperpInt = getAparaAperpTimeInt( Apara, Aperp, gamma, deltaGamma, deltaM, Phi_s, delta_para, delta_perp, timeMin, timeMax, tag );
        ReAzeroAparaInt = getAzeroAparaTimeInt( Azero, Apara, gamma, deltaGamma, deltaM, Phi_s, delta_para, timeMin, timeMax, tag );
        ImAzeroAperpInt = getAzeroAperpTimeInt( Azero, Aperp, gamma, deltaGamma, deltaM, Phi_s, delta_perp, timeMin, timeMax, tag );
	return;
}

double Bs2JpsiPhi::getEvenTimeComponentInt( double gamma
                                        , double deltaGamma
                                        , double deltaM
                                        , double Phi_s
                                        , double tlo
                                        , double thi
                                        , int tag ) const
{
	// Calculate some things that we need for the integral;
	double cosPhis = cos(Phi_s);
	double sinPhis = sin(Phi_s);
	double gammaL  = gamma + deltaGamma/2.;
	double gammaH  = gamma - deltaGamma/2.;
	double expLint = getExpInt( gammaL, tlo, thi );
	double expHint = getExpInt( gammaH, tlo, thi );
	double expSinInt = getExpSinInt( gamma, deltaM, tlo, thi );	

	if( tlo < 0. ) tlo = 0. ;

      	double result = (1. + cosPhis)*expLint + (1. - cosPhis)*expHint + tag * (2. * sinPhis)*expSinInt;
      	return result;
}

double Bs2JpsiPhi::getOddTimeComponentInt( double gamma
                                        , double deltaGamma
                                        , double deltaM
                                        , double Phi_s
                                        , double tlo
                                        , double thi
                                        , int tag ) const
{
        // Calculate some things that we need for the integral;
        double cosPhis = cos(Phi_s);
        double sinPhis = sin(Phi_s);
        double gammaL  = gamma + deltaGamma/2.;
        double gammaH  = gamma - deltaGamma/2.;
        double expLint = getExpInt( gammaL, tlo, thi );
        double expHint = getExpInt( gammaH, tlo, thi );
        double expSinInt = getExpSinInt( gamma, deltaM, tlo, thi );

        if( tlo < 0. ) tlo = 0. ;

        double result = (1. - cosPhis)*expLint + (1. + cosPhis)*expHint - tag * (2. * sinPhis)*expSinInt;
        return result;
}

double Bs2JpsiPhi::getAzeroAzeroTimeInt( double Azero_sq
                                        , double gamma
                                        , double deltaGamma
                                        , double deltaM
                                        , double Phi_s
                                        , double timeMin
                                        , double timeMax
					, int tag )  const
{

        double evenComponentInt = getEvenTimeComponentInt(gamma, deltaGamma, deltaM, Phi_s, timeMin, timeMax, tag );
        return Azero_sq*evenComponentInt;
}

double Bs2JpsiPhi::getAparaAparaTimeInt( double Apara_sq
                                        , double gamma
                                        , double deltaGamma
                                        , double deltaM
                                        , double Phi_s
                                        , double timeMin
                                        , double timeMax
					, int tag )  const
{
        double evenComponentInt = getEvenTimeComponentInt(gamma, deltaGamma, deltaM, Phi_s, timeMin, timeMax, tag);
	return Apara_sq*evenComponentInt;
}

double Bs2JpsiPhi::getAperpAperpTimeInt( double Aperp_sq
                                        , double gamma
                                        , double deltaGamma
                                        , double deltaM
                                        , double Phi_s
                                        , double timeMin
                                        , double timeMax
					, int tag )  const
{
        double oddComponentInt = getOddTimeComponentInt(gamma, deltaGamma, deltaM, Phi_s, timeMin, timeMax, tag);
        return Aperp_sq*oddComponentInt;
}

double Bs2JpsiPhi::getAparaAperpTimeInt( double Apara
					, double Aperp
                                        , double gamma
                                        , double deltaGamma
                                        , double deltaM
                                        , double Phi_s
                                        , double delta_para
                                        , double delta_perp
                                        , double timeMin
                                        , double timeMax
                                        , int tag ) const
{
	double sinDelta = sin(delta_para - delta_perp);
	double cosDelta = cos(delta_para - delta_perp);
	double sinPhis  = sin(Phi_s);
	double cosPhis  = cos(Phi_s);
        double gammaL   = gamma + deltaGamma/2.;
        double gammaH   = gamma - deltaGamma/2.;

	double expLint   = getExpInt(gammaL, timeMin, timeMax);
	double expHint   = getExpInt(gammaH, timeMin, timeMax);
	double expSinInt = getExpSinInt(gamma, deltaM, timeMin, timeMax);
	double expCosInt = getExpCosInt(gamma, deltaM, timeMin, timeMax);
	      	
	double result = 2.*( sinDelta*expCosInt - cosDelta*cosPhis*expSinInt )*tag
        		-  ( expHint - expLint )*cosDelta*sinPhis;
        return Apara*Aperp*result;
}

double Bs2JpsiPhi::getAzeroAparaTimeInt( double Azero
					, double Apara
                                        , double gamma
                                        , double deltaGamma
                                        , double deltaM
                                        , double Phi_s
                                        , double delta_para
                                        , double timeMin
                                        , double timeMax
					, int tag )  const
{
	double result = cos(delta_para)*getEvenTimeComponentInt(gamma, deltaGamma, deltaM, Phi_s, timeMin, timeMax, tag);
        return Azero*Apara*result;
}

double Bs2JpsiPhi::getAzeroAperpTimeInt( double Azero
					, double Aperp
                                        , double gamma
                                        , double deltaGamma
                                        , double deltaM
                                        , double Phi_s
                                        , double delta_perp
                                        , double timeMin
                                        , double timeMax
                                        , int tag ) const
{
        double sinDelta = sin(delta_perp);
        double cosDelta = cos(delta_perp);
        double sinPhis  = sin(Phi_s);
        double cosPhis  = cos(Phi_s);
        double gammaL   = gamma + deltaGamma/2.;
        double gammaH   = gamma - deltaGamma/2.;

        double expLint   = getExpInt(gammaL, timeMin, timeMax);
        double expHint   = getExpInt(gammaH, timeMin, timeMax);
        double expSinInt = getExpSinInt(gamma, deltaM, timeMin, timeMax);
        double expCosInt = getExpCosInt(gamma, deltaM, timeMin, timeMax);

        double result = 2.*( sinDelta*expCosInt - cosDelta*cosPhis*expSinInt )*tag
                        - ( expHint - expLint )*cosDelta*sinPhis;

        return Azero*Aperp*result;
}

// The following are some helper functions that are useful
// when doing the PDF integration.

// Integral of exp( -gamma*t ) from t1 to t2
double Bs2JpsiPhi::getExpInt( double G, double t1, double t2 ) const
{
        return (1./G) * ( exp(-G*t1) - exp(-G*t2) ) ;
}

// Integral of exp( -gamma*t ) * cos( dm*t )  from t1 to t2
double Bs2JpsiPhi::getExpCosInt( double G, double dm, double t1, double t2 ) const
{
        return (1./( G*G + dm*dm )) * (
                                ( exp(-G*t1)* (G*cos(dm*t1) - dm*sin(dm*t1)))
                               -( exp(-G*t2)* (G*cos(dm*t2) - dm*sin(dm*t2)))
                                );
}

// Integral of exp( -gamma*t ) * sin( dm*t )  from t1 to t2
double Bs2JpsiPhi::getExpSinInt( double G, double dm, double t1, double t2 ) const
{
        return (1./(G*G + dm*dm)) * (
                                ( exp(-G*t1)* (G*sin(dm*t1) + dm*cos(dm*t1)))
                               -( exp(-G*t2)* (G*sin(dm*t2) + dm*cos(dm*t2)))
                                );
}
// OK, we have now done the time integration part

// Now pull all of the integrals together...
double Bs2JpsiPhi::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
        // The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
        double AzeroAzeroInt, AparaAparaInt, AperpAperpInt;
        double ImAparaAperpInt, ReAzeroAparaInt, ImAzeroAperpInt;
        getTimeAmplitudeIntegrals( AzeroAzeroInt, AparaAparaInt, AperpAperpInt
                                , ImAparaAperpInt, ReAzeroAparaInt, ImAzeroAperpInt
				, measurement
                                , boundary);
   
	return 0.5*( AzeroAzeroInt + AparaAparaInt + AperpAperpInt ); // Angle factors normalised to 1
}

void Bs2JpsiPhi::getPhysicsParameters( double & gamma
					, double & deltaGamma
					, double & deltaM
					, double & Phi_s
				    	, double & Azero_sq
				    	, double & Apara_sq
				    	, double & Aperp_sq
				    	, double & delta_zero 
				    	, double & delta_para
				    	, double & delta_perp)
{
        // Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
    	deltaGamma = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
        deltaM     = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
        Phi_s      = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
        Azero_sq   = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
        //Apara_sq   = allParameters.GetPhysicsParameter( Apara_sqName )->GetValue();
        Aperp_sq   = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
        delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
        delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
        delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();

	Apara_sq = 1 - Azero_sq - Aperp_sq;

        return;
}
