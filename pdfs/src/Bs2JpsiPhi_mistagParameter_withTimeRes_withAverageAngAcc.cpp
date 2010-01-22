// $Id: Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi with average angular acceptance input as a paramter
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2010-01-12
 */

#include "Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc() : 
	// Physics parameters
	gammaName     ( "gamma" )
	, deltaGammaName( "deltaGamma" )
	, deltaMName    ( "deltaM")
	, Phi_sName     ( "Phi_s")
	, Azero_sqName  ( "Azero_sq" )
	, Aperp_sqName  ( "Aperp_sq" )
	, delta_zeroName( "delta_zero" )
	, delta_paraName( "delta_para" )
	, delta_perpName( "delta_perp" )
	, angAccI1Name	( "angAccI1" )
	, angAccI2Name	( "angAccI2" )
	, angAccI3Name	( "angAccI3" )
	, angAccI4Name	( "angAccI4" )
	, angAccI5Name	( "angAccI5" )
	, angAccI6Name	( "angAccI6" )
	, mistagName	( "mistag" )
	, timeResName	( "timeRes" )

	// Observables
	, timeName	( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	( "phi" )
	, cosPsiName	( "cosPsi" )
	, tagName	( "tag" )
	//, timeres	( "resolution" )
	, normalisationCacheValid(false)
, evaluationCacheValid(false)
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );
	allObservables.push_back( tagName );

	// Need to think about additional parameters like
	// event-by-event propertime resolution and acceptance.
	// This will require event-by-event PDF normalisation,
	// but we are already doing this for tagging.

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
	parameterNames.push_back( mistagName );
	parameterNames.push_back( timeResName );
	parameterNames.push_back( angAccI1Name );
	parameterNames.push_back( angAccI2Name );
	parameterNames.push_back( angAccI3Name );
	parameterNames.push_back( angAccI4Name );
	parameterNames.push_back( angAccI5Name );
	parameterNames.push_back( angAccI6Name );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::~Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
        // Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
        gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
        deltaGamma = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
        deltaMs    = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
        Phi_s      = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
        Azero_sq   = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
        Aperp_sq   = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
        delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
        delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
        delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
        omega = allParameters.GetPhysicsParameter( mistagName )->GetValue();
        timeRes = allParameters.GetPhysicsParameter( timeResName )->GetValue();
        angAccI1 = allParameters.GetPhysicsParameter( angAccI1Name )->GetValue();
        angAccI2 = allParameters.GetPhysicsParameter( angAccI2Name )->GetValue();
        angAccI3 = allParameters.GetPhysicsParameter( angAccI3Name )->GetValue();
        angAccI4 = allParameters.GetPhysicsParameter( angAccI4Name )->GetValue();
        angAccI5 = allParameters.GetPhysicsParameter( angAccI5Name )->GetValue();
        angAccI6 = allParameters.GetPhysicsParameter( angAccI6Name )->GetValue();

        Apara_sq = 1 - Azero_sq - Aperp_sq;
	AparaAperp = sqrt(Apara_sq)*sqrt(Aperp_sq);
        AzeroApara = sqrt(Azero_sq)*sqrt(Apara_sq);
        AzeroAperp = sqrt(Azero_sq)*sqrt(Aperp_sq);

	return isOK;
}

//Return a list of parameters not to be integrated
vector<string> Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::GetDoNotIntegrateList()
{
	vector<string> doNotIntList;
	doNotIntList.push_back(mistagName);
	return doNotIntList;
}

//Calculate the function value
double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::Evaluate(DataPoint * measurement)
{
	// Does not make sense to evaluate this PDF for time < 0
	double time = measurement->GetObservable( timeName )->GetValue();	
        double cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
        double phi      = measurement->GetObservable( phiName )->GetValue();
        double cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();
	int q = (int)measurement->GetObservable( tagName )->GetValue();

	if ( time < 0. ) return 0.;

	// The angular functions f1->f6 as defined in roadmap Table 1.
	double f1, f2, f3, f4, f5, f6;
	Mathematics::getBs2JpsiPhiAngularFunctions( f1, f2, f3, f4, f5, f6, cosTheta, phi, cosPsi );	

	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59
	// First for the B
	double AzeroAzeroB, AparaAparaB, AperpAperpB;
	double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB;
	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB
			, ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
			, measurement, 1);

	// Now for the Bbar
	double AzeroAzeroBbar, AparaAparaBbar, AperpAperpBbar;
	double ImAparaAperpBbar, ReAzeroAparaBbar, ImAzeroAperpBbar;
	getTimeDependentAmplitudes( AzeroAzeroBbar, AparaAparaBbar, AperpAperpBbar
			, ImAparaAperpBbar, ReAzeroAparaBbar, ImAzeroAperpBbar
			, measurement, -1);

	// Flavour tagging
	double epsilon[3];
	epsilon[0] = omega;
	epsilon[1] = 0.5;
	epsilon[2] = (1.0 - omega);
	double w1  = epsilon[q + 1];
	double w2  = 1.0 - w1;

	//W+
	double v1 = f1 * AzeroAzeroB
		+ f2 * AparaAparaB
		+ f3 * AperpAperpB
		+ f4 * ImAparaAperpB
		+ f5 * ReAzeroAparaB
		+ f6 * ImAzeroAperpB; 

	//W-
	double v2 = f1 * AzeroAzeroBbar
		+ f2 * AparaAparaBbar
		+ f3 * AperpAperpBbar
		+ f4 * ImAparaAperpBbar
		+ f5 * ReAzeroAparaBbar
		+ f6 * ImAzeroAperpBbar;  

	return ( w1*v1 + w2*v2 );
}


double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	// Now need to know the tag and the mistag
	int q = (int)measurement->GetObservable( tagName )->GetValue();

	double epsilon[3];
	epsilon[0] = omega;
	epsilon[1] = 0.5;
	epsilon[2] = (1.0 - omega);
	double w1  = epsilon[q + 1];
	double w2  = 1.0 - w1;

	if (!normalisationCacheValid)
	{
		// The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
		getTimeAmplitudeIntegrals(  cachedAzeroAzeroIntB
				, cachedAparaAparaIntB
				, cachedAperpAperpIntB
				, cachedAparaAperpIntB
				, cachedAzeroAparaIntB
				, cachedAzeroAperpIntB
				, boundary
				, 1);

		getTimeAmplitudeIntegrals(  cachedAzeroAzeroIntBbar
				, cachedAparaAparaIntBbar
				, cachedAperpAperpIntBbar
				, cachedAparaAperpIntBbar
				, cachedAzeroAparaIntBbar
				, cachedAzeroAperpIntBbar
				, boundary
				, -1);

		normalisationCacheValid = true;
	}

	double v1 = cachedAzeroAzeroIntB * angAccI1
		+ cachedAparaAparaIntB * angAccI2
		+ cachedAperpAperpIntB * angAccI3 
		+ cachedAparaAperpIntB * angAccI4 
		+ cachedAzeroAparaIntB * angAccI5 
		+ cachedAzeroAperpIntB * angAccI6; 

	double v2 = cachedAzeroAzeroIntBbar * angAccI1
		+ cachedAparaAparaIntBbar * angAccI2
		+ cachedAperpAperpIntBbar * angAccI3 
		+ cachedAparaAperpIntBbar * angAccI4 
		+ cachedAzeroAparaIntBbar * angAccI5 
		+ cachedAzeroAperpIntBbar * angAccI6; 

	return w1*v1 + w2*v2;
}

void Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::getTimeDependentAmplitudes(  double & AzeroAzero
		, double & AparaApara
		, double & AperpAperp
		, double & ImAparaAperp
		, double & ReAzeroApara
		, double & ImAzeroAperp
		, DataPoint * measurement
		, int Btype )
{
	// Observable
	double time = measurement->GetObservable( timeName )->GetValue();

	// Quantities depending only on physics parameters can be cached
	if ( !evaluationCacheValid )
	{
		cachedAzero = sqrt( Azero_sq );
		cachedApara = sqrt( Apara_sq );
		cachedAperp = sqrt( Aperp_sq );

		cachedsinDeltaPerpPara	= sin( delta_perp - delta_para );
		cachedcosDeltaPerpPara	= cos( delta_perp - delta_para );
		cachedsinDeltaPerp	= sin( delta_perp );
		cachedcosDeltaPerp	= cos( delta_perp );
		cachedcosDeltaPara	= cos( delta_para );

		cachedsinPhis = sin( Phi_s );
		cachedcosPhis = cos( Phi_s );
		evaluationCacheValid = true;
	}

	// Quantities depending on time cannot be cached
	double expGT = exp( -gamma*time );
	double coshDeltaGammaT = cosh( deltaGamma*time/2.);
	double sinhDeltaGammaT = sinh( deltaGamma*time/2.);
	double sinDeltaMsT = sin( deltaMs*time );
	double cosDeltaMsT = cos( deltaMs*time );

	// Now calculate the amplitudes
	AzeroAzero = Azero_sq * expGT * ( coshDeltaGammaT - cachedcosPhis * sinhDeltaGammaT + Btype * cachedsinPhis * sinDeltaMsT ); 
	AparaApara = Apara_sq * expGT * ( coshDeltaGammaT - cachedcosPhis * sinhDeltaGammaT + Btype * cachedsinPhis * sinDeltaMsT ); 
	AperpAperp = Aperp_sq * expGT * ( coshDeltaGammaT + cachedcosPhis * sinhDeltaGammaT - Btype * cachedsinPhis * sinDeltaMsT ); 

	ImAparaAperp = cachedApara*cachedAperp * expGT * ( - cachedcosDeltaPerpPara * cachedsinPhis * sinhDeltaGammaT 
			+ Btype * cachedsinDeltaPerpPara * cosDeltaMsT
			- Btype * cachedcosDeltaPerpPara * cachedcosPhis * sinDeltaMsT );

	ReAzeroApara = cachedAzero*cachedApara * expGT * cachedcosDeltaPara * ( coshDeltaGammaT - cachedcosPhis * sinhDeltaGammaT
			+ Btype * cachedsinPhis * sinDeltaMsT );

	ImAzeroAperp = cachedAzero*cachedAperp * expGT * ( - cachedcosDeltaPerp * cachedsinPhis * sinhDeltaGammaT
			+ Btype * cachedsinDeltaPerp * cosDeltaMsT 
			- Btype * cachedcosDeltaPerp * cachedcosPhis * sinDeltaMsT );
	return;
}

void Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::getTimeAmplitudeIntegrals( double & AzeroAzeroInt
		, double & AparaAparaInt
		, double & AperpAperpInt
		, double & AparaAperpInt
		, double & AzeroAparaInt
		, double & AzeroAperpInt
		, PhaseSpaceBoundary * boundary
		, int Btype)
{
	double tlow = 0.;
	double thigh = 0.;
	IConstraint * timeBound = boundary->GetConstraint("time");
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		AzeroAzeroInt = -999.;
		AparaAparaInt = -999.;
		AperpAperpInt = -999.;
		return;
	}
	else
	{
		tlow = timeBound->GetMinimum();
		if (tlow < 0.) tlow = 0.;
		thigh = timeBound->GetMaximum();
	}

	double G_H    = gamma - 0.5*deltaGamma;
	double G_L    = gamma + 0.5*deltaGamma;

	double tauH   = (1.0 / G_H);
	double tauL   = (1.0 / G_L);
	double tauBar = (1.0 / gamma);

	// Need to introduce some caching of these values (cosPhis) here

	AzeroAzeroInt = getAzeroAzeroInt( tlow, thigh, tauL, tauH, tauBar, Btype);
	AparaAparaInt = getAparaAparaInt( tlow, thigh, tauL, tauH, tauBar, Btype);
	AperpAperpInt = getAperpAperpInt( tlow, thigh, tauL, tauH, tauBar, Btype);
	AparaAperpInt = getAparaAperpInt( tlow, thigh, Btype);
	AzeroAparaInt = getAzeroAparaInt( tlow, thigh, Btype);
	AzeroAperpInt = getAzeroAperpInt( tlow, thigh, Btype);
	// The interference terms above are expressed in terms of gamma, deltaGamma, not tauL, H. 
	// I should eventually get round to have everything in terms of gamma etc. 
	return;
}

/*
AzeroAzero[t_] := Exp[-gamma*t] * (   Cosh[deltaGamma*t/2]              
                                    - cosPhis * Sinh[deltaGamma*t/2]      
                                    + Btype * sinPhis * Sin[deltaMs*t] );
AzeroAzeroInt = CForm[FullSimplify[ Integrate[ AzeroAzero[t], {t, tmin, tmax}]]]
*/
inline double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::getAzeroAzeroInt(double tmin, double tmax, 
		double tauL, double tauH, double tauBar, int Btype)
{

	double gammadeltaMs = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (deltaMs*deltaMs)) );
	double cosPhi_s = cos(Phi_s);
	double sinPhi_s = sin(Phi_s);

	double valB = 0.5*Azero_sq*((1.0+cosPhi_s)*(tauL)*exp(-(1.0/tauL)*tmax)
			+(1.0-cosPhi_s)*(tauH)*exp(-(1.0/tauH)*tmax)
			+ Btype*2.0*gammadeltaMs*sinPhi_s*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(deltaMs*tmax)
					+deltaMs*cos(deltaMs*tmax))));

	double valA = 0.5*Azero_sq*((1.0+cosPhi_s)*(tauL)*exp(-(1.0/tauL)*tmin)
			+(1.0-cosPhi_s)*(tauH)*exp(-(1.0/tauH)*tmin)
			+ Btype*2.0*gammadeltaMs*sinPhi_s*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(deltaMs*tmin)
					+deltaMs*cos(deltaMs*tmin))));

	return (valA-valB);
}

/*
AparaApara[t_] := Exp[-gamma*t] * (   Cosh[deltaGamma*t/2]              
                                    - cosPhis * Sinh[deltaGamma*t/2]      
                                    + Btype * sinPhis * Sin[deltaMs*t] );
AparaAparaInt = CForm[FullSimplify[ Integrate[ AparaApara[t], {t, tmin, tmax}]]]
*/
inline double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::getAparaAparaInt(double tmin, double tmax,
		double tauL, double tauH, double tauBar, int Btype)
{

	double gammadeltaMs = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (deltaMs*deltaMs)) );
	double cosPhi_s = cos(Phi_s);
	double sinPhi_s = sin(Phi_s);

	double valB = 0.5*Apara_sq*((1.0+cosPhi_s)*(tauL)*exp(-(1.0/tauL)*tmax)
			+(1.0-cosPhi_s)*(tauH)*exp(-(1.0/tauH)*tmax)
			+ Btype*2.0*gammadeltaMs*sinPhi_s*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(deltaMs*tmax)
					+deltaMs*cos(deltaMs*tmax))));

	double valA = 0.5*Apara_sq*((1.0+cosPhi_s)*(tauL)*exp(-(1.0/tauL)*tmin)
			+(1.0-cosPhi_s)*(tauH)*exp(-(1.0/tauH)*tmin)
			+ Btype*2.0*gammadeltaMs*sinPhi_s*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(deltaMs*tmin)
					+deltaMs*cos(deltaMs*tmin))));
	
	return (valA-valB);	
}

/*
AperpAperp[t_] := Exp[-gamma*t] * (   Cosh[deltaGamma*t/2]              
                                    + cosPhis * Sinh[deltaGamma*t/2]      
                                    - Btype * sinPhis * Sin[deltaMs*t] );
AperpAperpInt = CForm[FullSimplify[ Integrate[ AperpAperp[t], {t, tmin, tmax}]]]
*/
inline double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::getAperpAperpInt(double tmin, double tmax,
		double tauL, double tauH, double tauBar, int Btype)
{

	double gammadeltaMs = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (deltaMs*deltaMs)) );
	double cosPhi_s = cos(Phi_s);
	double sinPhi_s = sin(Phi_s);

	double valB =  0.5*Aperp_sq*((1.0-cosPhi_s)*(tauL)*exp(-(1.0/tauL)*tmax)
			+(1.0+cosPhi_s)*(tauH)*exp(-(1.0/tauH)*tmax)
			- Btype*2.0*gammadeltaMs*sinPhi_s*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(deltaMs*tmax)
					+deltaMs*cos(deltaMs*tmax))));

	double valA =  0.5*Aperp_sq*((1.0-cosPhi_s)*(tauL)*exp(-(1.0/tauL)*tmin)
			+(1.0+cosPhi_s)*(tauH)*exp(-(1.0/tauH)*tmin)
			- Btype*2.0*gammadeltaMs*sinPhi_s*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(deltaMs*tmin)
					+deltaMs*cos(deltaMs*tmin))));

	return (valA-valB);
}

/*
AparaAperp[t_] := Exp[-gamma*t] * (-cosDeltaPerpMinusPara * sinPhis * Sinh[deltaGamma*t/2]              
                                    + Btype * sinDeltaPerpMinusPara * Cos[deltaMs*t]      
                                    - Btype * cosDeltaPerpMinusPara * cosPhis * Sinh[deltaGamma*t/2] );
AparaAperpInt = CForm[FullSimplify[ Integrate[ AparaAperp[t], {t, tmin, tmax}]]]
*/
inline double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::getAparaAperpInt(double tmin, double tmax, int Btype) 
{
	double cosPhis = cos(Phi_s);
	double sinPhis = sin(Phi_s);
	double cosDeltaPerpMinusPara = cos(delta_perp - delta_para);
	double sinDeltaPerpMinusPara = sin(delta_perp - delta_para);

	double v = (-2. * cosDeltaPerpMinusPara * deltaGamma * (Btype * cosPhis + sinPhis) * cosh((deltaGamma*tmax)/2.))/
     (exp(gamma * tmax) * (deltaGamma * deltaGamma - 4. * gamma * gamma)) + 
    (-2. * cosDeltaPerpMinusPara * deltaGamma * exp(gamma * tmax)*
        (deltaMs * deltaMs + gamma * gamma)*(Btype * cosPhis + sinPhis)*
        cosh((deltaGamma * tmin)/2.) + 
       Btype * (deltaGamma * deltaGamma - 4. * gamma * gamma) * sinDeltaPerpMinusPara*
        (exp(gamma * tmin) * (gamma * cos(deltaMs * tmax) - 
             deltaMs * sin(deltaMs * tmax)) + 
          exp(gamma * tmax) * (-(gamma * cos(deltaMs * tmin)) + 
             deltaMs * sin(deltaMs * tmin))) + 
       4. * cosDeltaPerpMinusPara * gamma * (deltaMs * deltaMs + gamma * gamma)*
        (Btype * cosPhis + sinPhis) *
        (exp(gamma*tmin) * sinh((deltaGamma*tmax)/2.) - 
          exp(gamma*tmax) * sinh((deltaGamma*tmin)/2.)))/
     (exp(gamma*(tmax + tmin)) * (deltaMs * deltaMs + gamma * gamma)*
       (-deltaGamma * deltaGamma + 4. * gamma * gamma));

	return AparaAperp * v;

	/*
	double valB = 0.5*k0*cos(delta_perp - delta_para)* ( tauH*exp(-(1.0/tauH)*tmax)
                                               + tauL*exp(-(1.0/tauL)*tmax)
                                               - cosPhi_s* ( tauH*exp(-(1.0/tauH)*tmax) 
                                                            - tauL*exp(-(1.0/tauL)*tmax) )
                                               + 2.0 * sinPhi_s * exp (-(1.0/tauBar)*tmax) * gammadeltaMs 
                                               * ( deltaMs * cos(deltaMs*tmax) + (1.0/tauBar)*sin(deltaMs*tmax)) );
  
	double valA = 0.5*k0*cos(delta_perp - delta_para)* ( tauH*exp(-(1.0/tauH)*tmin)
                                               + tauL*exp(-(1.0/tauL)*tmin)
                                               - cosPhi_s* ( tauH*exp(-(1.0/tauH)*tmin) 
                                                            - tauL*exp(-(1.0/tauL)*tmin) )
                                               + 2.0 * sinPhi_s * exp (-(1.0/tauBar)*tmin) * gammadeltaMs 
                                               * ( deltaMs * cos(deltaMs*tmin) + (1.0/tauBar)*sin(deltaMs*tmin)) );
  
  	return (valA-valB);
	*/
}


/*
AzeroApara[t_] := Exp[-gamma*t] * cosDeltaPara * (Cosh[deltaGamma*t/2]              
                                    - cosPhis * sinh[deltaGamma*t/2]      
                                    + Btype * sinPhis * Sin[deltaMs*t] );
AzeroAparaInt = CForm[FullSimplify[ Integrate[ AzeroApara[t], {t, tmin, tmax}]]]
*/
inline double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::getAzeroAparaInt(double tmin, double tmax, int Btype) 
{
	double cosPhis = cos(Phi_s);
	double sinPhis = sin(Phi_s);
	double cosDeltaPara = cos(delta_para);

	double v = 
		cosDeltaPara * (
			(2. * (-(cosPhis * deltaGamma) + 2. * gamma) * cosh((deltaGamma * tmax)/2.))/
       				(exp(gamma * tmax) * (deltaGamma * deltaGamma - 4. * gamma * gamma)) 
			+ (2. * exp(gamma * tmax) * (-(cosPhis*deltaGamma) + 2. * gamma) * (deltaMs * deltaMs + gamma * gamma) * cosh((deltaGamma * tmin)/2.)
			+ Btype * (deltaGamma * deltaGamma - 4. * gamma * gamma) * sinPhis * (exp(gamma * tmin) * (deltaMs * cos(deltaMs * tmax)
			+ gamma * sin(deltaMs*tmax)) - exp(gamma * tmax) * (deltaMs * cos(deltaMs * tmin) + gamma * sin(deltaMs * tmin)))
			- 2. * (deltaGamma - 2. * cosPhis * gamma) * (deltaMs * deltaMs + gamma * gamma) * (exp(gamma * tmin)* sinh((deltaGamma * tmax)/2.)
			- exp(gamma * tmax) * sinh((deltaGamma * tmin)/2.)))/(exp(gamma * (tmax + tmin)) * (deltaMs * deltaMs + gamma * gamma)*(-deltaGamma * deltaGamma + 4. * gamma * gamma))
			);
	return AzeroApara * v;
	/*
	double valB = 0.5*k0*((1.0+cosPhi_s)*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0-cosPhi_s)*(tauH)*exp(-(1.0/tauH)*tmax)
                        + Btype*2.0*gammadeltaMs*sinPhi_s*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(deltaMs*tmax)
                                        +deltaMs*cos(deltaMs*tmax))));

        double valA = 0.5*k0*((1.0+cosPhi_s)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cosPhi_s)*(tauH)*exp(-(1.0/tauH)*tmin)
                        + Btype*2.0*gammadeltaMs*sinPhi_s*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(deltaMs*tmin)
                                        +deltaMs*cos(deltaMs*tmin))));

        return (valA-valB);
	*/
}

/*
AzeroAperp[t_] := Exp[-gamma*t] * ( - cosDeltaPerp * sinPhis * Sinh[deltaGamma*t/2]              
                                    + Btype * sinDeltaPerp * Cos[deltaMs*t]      
                                    - Btype * cosDeltaPerp * cosPhis * Sin[deltaMs*t] );
AzeroAperpInt = CForm[FullSimplify[ Integrate[ AzeroAperp[t], {t, tmin, tmax}]]]
*/
inline double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::getAzeroAperpInt(double tmin, double tmax, int Btype)
{
	double cosPhis = cos(Phi_s);
	double sinPhis = sin(Phi_s);
	double cosDeltaPerp = cos(delta_perp);
	double sinDeltaPerp = sin(delta_perp);

        double v = 
	  (Btype * sinDeltaPerp * ((-(gamma * cos(deltaMs * tmax)) + deltaMs * sin(deltaMs * tmax))/exp(gamma * tmax)
	+ (gamma   * cos(deltaMs*tmin) - deltaMs * sin(deltaMs * tmin))/exp(gamma * tmin)))/(deltaMs * deltaMs + gamma * gamma)
	- (Btype * cosDeltaPerp * cosPhis * (-((deltaMs * cos(deltaMs * tmax) + gamma * sin(deltaMs * tmax))/exp(gamma * tmax)) 
        + (deltaMs * cos(deltaMs*tmin) + gamma   * sin(deltaMs * tmin))/exp(gamma * tmin)))/(deltaMs * deltaMs + gamma * gamma)
	- (2. * cosDeltaPerp * sinPhis * ((deltaGamma * cosh((deltaGamma * tmax)/2.) 
		+ 2. * gamma * sinh((deltaGamma * tmax)/2.))/exp(gamma * tmax) - (deltaGamma * cosh(( deltaGamma * tmin)/2.)
		+ 2. * gamma * sinh((deltaGamma * tmin)/2.))/exp(gamma * tmin)))/(deltaGamma * deltaGamma - 4.*gamma*gamma);

	return AzeroAperp * v;
}
