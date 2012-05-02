// $Id: Bs2JpsiPhi_mistagObservable.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhi_mistagObservable Bs2JpsiPhi_mistagObservable.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-07-30
 */

#include "Bs2JpsiPhi_mistagObservable.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
Bs2JpsiPhi_mistagObservable::Bs2JpsiPhi_mistagObservable() : 
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

	// Observables
	, timeName	( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	( "phi" )
	, cosPsiName	( "cosPsi" )
	, tagName	( "tag" )
	, mistagName	( "mistag" )
	//, timeres	( "resolution" )
	, normalisationCacheValid(false)
	, evaluationCacheValid(false)
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhi_mistagObservable::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );
	allObservables.push_back( tagName );
	allObservables.push_back( mistagName );

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
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2JpsiPhi_mistagObservable::~Bs2JpsiPhi_mistagObservable()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2JpsiPhi_mistagObservable::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	return allParameters.SetPhysicsParameters(NewParameterSet);
}

//Return a list of parameters not to be integrated
vector<string> Bs2JpsiPhi_mistagObservable::GetDoNotIntegrateList()
{
	vector<string> mistagList;
	mistagList.push_back(mistagName);
	return mistagList;
}

//Calculate the function value
double Bs2JpsiPhi_mistagObservable::Evaluate(DataPoint * measurement)
{
	// Does not make sense to evaluate this PDF for time < 0
	double time = measurement->GetObservable( timeName )->GetValue();	
	if ( time < 0. ) return 0.;

	// The angular functions f1->f6 as defined in roadmap Table 1.
	double f1, f2, f3, f4, f5, f6;
	getAngularFunctions( f1, f2, f3, f4, f5, f6, measurement );	

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

	// Now need to know the tag and the mistag
	int q = (int)measurement->GetObservable( tagName )->GetValue();
	double omega = measurement->GetObservable( mistagName )->GetValue();	

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


double Bs2JpsiPhi_mistagObservable::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	// Now need to know the tag and the mistag
	int q = (int)measurement->GetObservable( tagName )->GetValue();
	double omega = measurement->GetObservable( mistagName )->GetValue();	

	double epsilon[3];
	epsilon[0] = omega;
	epsilon[1] = 0.5;
	epsilon[2] = (1.0 - omega);
	double w1  = epsilon[q + 1];
	double w2  = 1.0 - w1;

	if (!normalisationCacheValid)
	{
		// The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
		double AzeroAzeroIntB, AparaAparaIntB, AperpAperpIntB;
		getTimeAmplitudeIntegrals(  AzeroAzeroIntB
				, AparaAparaIntB
				, AperpAperpIntB
				, boundary
				, 1);

		double AzeroAzeroIntBbar, AparaAparaIntBbar, AperpAperpIntBbar;
		getTimeAmplitudeIntegrals(  AzeroAzeroIntBbar
				, AparaAparaIntBbar
				, AperpAperpIntBbar
				, boundary
				, -1);

		cachedv1 = AzeroAzeroIntB + AparaAparaIntB + AperpAperpIntB;
		cachedv2 = AzeroAzeroIntBbar + AparaAparaIntBbar + AperpAperpIntBbar;
		normalisationCacheValid = true;
	}

	return ( w1*cachedv1 + w2*cachedv2 );
}

void Bs2JpsiPhi_mistagObservable::getAngularFunctions( double & f1
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

	double sinTheta  = sqrt(1. - cosTheta*cosTheta);
	double sinPsi    = sqrt(1. - cosPsi*cosPsi);

	double cosPhi    = cos(phi);
	double sinPhi	 = sin(phi);

	double sin2Theta = 2.*sinTheta*cosTheta;
	double sin2Psi	 = 2.*sinPsi*cosPsi;
	double sin2Phi	 = 2.*sinPhi*cosPhi;

	double norm = 9./32./TMath::Pi(); 
	//this is the factor that drops out when you integrate over the angles
	// i.e., int(2*cospsi*cospsi*(1-(1-costh*costh)*cos(phi)*cos(phi)),cospsi=-1..1, costh=-1..1,phi=-Pi..Pi);
	// same factor for f1, f2, f3. The remaining terms f4, f5, f6 give 0
	// int(-(1-cospsi*cospsi)*2*sqrt(1-costh*costh)*costh*sin(phi),cospsi=-1..1, costh=-1..1,phi=-Pi..Pi); 
	f1 =  2.* cosPsi*cosPsi * ( 1. - sinTheta*sinTheta * cosPhi*cosPhi ) * norm;
	f2 =      sinPsi*sinPsi * ( 1. - sinTheta*sinTheta * sinPhi*sinPhi ) * norm;
	f3 =      sinPsi*sinPsi * sinTheta*sinTheta * norm;
	f4 = -1.* sinPsi*sinPsi * sin2Theta * sinPhi * norm;
	f5 = sin2Psi * sinTheta*sinTheta * sin2Phi/sqrt(2.) * norm;
	f6 = sin2Psi * sin2Theta * cosPhi/sqrt(2.) * norm;
	return;
}

void Bs2JpsiPhi_mistagObservable::getTimeDependentAmplitudes(  double & AzeroAzero
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

	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	double gamma, deltaGamma, deltaMs, Phi_s;
	double Azero_sq, Apara_sq, Aperp_sq;
	double delta_zero, delta_para, delta_perp;
	getPhysicsParameters( gamma, deltaGamma, deltaMs, Phi_s, Azero_sq, Apara_sq, Aperp_sq, delta_zero, delta_para, delta_perp);

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

void Bs2JpsiPhi_mistagObservable::getTimeAmplitudeIntegrals( double & AzeroAzeroInt
		, double & AparaAparaInt
		, double & AperpAperpInt
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

	// Physics parameters
	double gamma, deltaGamma, deltaMs, Phi_s;
	double Azero_sq, Apara_sq, Aperp_sq;
	double delta_zero, delta_para, delta_perp;
	getPhysicsParameters( gamma, deltaGamma, deltaMs, Phi_s, Azero_sq, Apara_sq, Aperp_sq, delta_zero, delta_para, delta_perp);

	double G_H    = gamma - 0.5*deltaGamma;
	double G_L    = gamma + 0.5*deltaGamma;

	double tauH   = (1.0 / G_H);
	double tauL   = (1.0 / G_L);
	double tauBar = (1.0 / gamma);

	AzeroAzeroInt = getAzeroAzeroInt( tlow, thigh, Azero_sq, tauL, tauH, tauBar, deltaMs, Phi_s, Btype);
	AparaAparaInt = getAparaAparaInt( tlow, thigh, Apara_sq, tauL, tauH, tauBar, deltaMs, Phi_s, Btype);
	AperpAperpInt = getAperpAperpInt( tlow, thigh, Aperp_sq, tauL, tauH, tauBar, deltaMs, Phi_s, Btype);
	// No contribution from interference terms here since they drop out when the angular integration is done.
	return;
}


inline double Bs2JpsiPhi_mistagObservable::getAzeroAzeroInt(double tmin, double tmax, 
		double k0, double tauL, double tauH, double tauBar, double Dms, double phis, int Btype)
{

	double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
	double cosphis = cos(phis);
	double sinphis = sin(phis);

	double valB = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmax)
			+(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmax)
			+ Btype*2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
					+Dms*cos(Dms*tmax))));

	double valA = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
			+(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
			+ Btype*2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
					+Dms*cos(Dms*tmin))));

	return (valA-valB);
}


inline double Bs2JpsiPhi_mistagObservable::getAparaAparaInt(double tmin, double tmax,
		double k0, double tauL, double tauH, double tauBar,double Dms, double phis, int Btype)
{

	double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
	double cosphis = cos(phis);
	double sinphis = sin(phis);

	double valB = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmax)
			+(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmax)
			+ Btype*2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
					+Dms*cos(Dms*tmax))));

	double valA = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
			+(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
			+ Btype*2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
					+Dms*cos(Dms*tmin))));

	return (valA-valB);	
}

inline double Bs2JpsiPhi_mistagObservable::getAperpAperpInt(double tmin, double tmax,
		double k0, double tauL, double tauH, double tauBar,double Dms, double phis, int Btype)
{

	double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
	double cosphis = cos(phis);
	double sinphis = sin(phis);

	double valB =  0.5*k0*((1.0-cosphis)*(tauL)*exp(-(1.0/tauL)*tmax)
			+(1.0+cosphis)*(tauH)*exp(-(1.0/tauH)*tmax)
			- Btype*2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
					+Dms*cos(Dms*tmax))));

	double valA =  0.5*k0*((1.0-cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
			+(1.0+cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
			- Btype*2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
					+Dms*cos(Dms*tmin))));

	return (valA-valB);
}

// Need to work out what interference terms these correspond to and which strong phase is which.
// These are only needed (I think) for calculating the projections.
inline double Bs2JpsiPhi_mistagObservable::A4def(double tmin, double tmax,
		double k0, double tauL, double tauH, double tauBar, 
		double Dms, double phis, double tphase1, double tphase2, int Btype)
{

	double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
	double cosphis = cos(phis);
	double sinphis = sin(phis);

	double valB = 0.5*k0*cos(tphase2-tphase1)* ( tauH*exp(-(1.0/tauH)*tmax)
			+ tauL*exp(-(1.0/tauL)*tmax)
			- cosphis* ( tauH*exp(-(1.0/tauH)*tmax) 
				- tauL*exp(-(1.0/tauL)*tmax) )
			+ Btype * 2.0 * sinphis * exp (-(1.0/tauBar)*tmax) * gammaDms 
			* ( Dms* cos(Dms*tmax) + (1.0/tauBar)*sin(Dms*tmax)) );

	double valA = 0.5*k0*cos(tphase2-tphase1)* ( tauH*exp(-(1.0/tauH)*tmin)
			+ tauL*exp(-(1.0/tauL)*tmin)
			- cosphis* ( tauH*exp(-(1.0/tauH)*tmin) 
				- tauL*exp(-(1.0/tauL)*tmin) )
			+ Btype * 2.0 * sinphis * exp (-(1.0/tauBar)*tmin) * gammaDms 
			* ( Dms* cos(Dms*tmin) + (1.0/tauBar)*sin(Dms*tmin)) );

	return (valA-valB);

}

inline double Bs2JpsiPhi_mistagObservable::A5def(double tmin, double tmax,
		double k0, double tauL, double tauH, double tauBar, 
		double Dms, double phis, double tphase1, double tphase2, int Btype)
{

	double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
	double cosphis = cos(phis);
	double sinphis = sin(phis);

	double valB = k0 * ( gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*cos(Dms*tmax) 
				- Dms*sin(Dms*tmax))
			- gammaDms*cos(tphase1)*cosphis*exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)
				*sin(Dms*tmax)
				+ Dms*cos(Dms*tmax))
			- Btype * 0.5*((tauH)*exp(-(1.0/tauH)*tmax)-(tauL)*exp(-(1.0/tauL)*tmax))
			*cos(tphase1)*sinphis);

	double valA = k0 * ( gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*cos(Dms*tmin)
				- Dms*sin(Dms*tmin))
			- gammaDms*cos(tphase1)*cosphis*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
				+ Dms*cos(Dms*tmin))
			- Btype * 0.5*((tauH)*exp(-(1.0/tauH)*tmin)-(tauL)*exp(-(1.0/tauL)*tmin))*cos(tphase1)*sinphis);

	return (valA-valB);

}



void Bs2JpsiPhi_mistagObservable::getPhysicsParameters( double & gamma
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
