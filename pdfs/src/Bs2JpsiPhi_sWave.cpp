/** @class Bs2JpsiPhi_sWave Bs2JpsiPhi_sWave.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi with sWave contributions
 *
 *  @author Conor Fitzpatrick
 *  @date 2009-10-27
 */

#include "Bs2JpsiPhi_sWave.h"
#include <iostream>
#include "math.h"
#include "TMath.h"


using std::cout;
using std::endl;
//Constructor
Bs2JpsiPhi_sWave::Bs2JpsiPhi_sWave() :
	// Physics parameters
	gammaName	 ( "gamma" )
	, deltaGammaName( "deltaGamma" )
	, deltaMName	( "deltaM")
	, Phi_sName	 ( "Phi_s")
	, Azero_sqName  ( "Azero_sq" )
	, Aperp_sqName  ( "Aperp_sq" )
	, As_sqName	("As_sq")
	, delta_zeroName( "delta_zero" )
	, delta_paraName( "delta_para" )
	, delta_perpName( "delta_perp" )
	, delta_sName	("delta_s")

//Note: arXiv:0908.3627v3 uses the following parameters instead:
// Ap_sq = Azero_sq + Aperp_sq + Apara_sq
// Rpara = Apara_sq/Ap_sq
// Rperp = Aperp_sq/Ap_sq
// Rs = As_sq / (As_sq + Ap_sq)
// To get from these to our params, note that:
// Azero_sq + Aperp_sq + Apara_sq = 1 = Ap_sq
// so:
// Rpara = Apara_sq
// Rperp = Aperp_sq
// Rs = As_sq / (As_sq + 1) => 	As_sq = (1 - Rs)/Rs


	// Observables
	, timeName	( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	( "phi" )
	, cosPsiName	( "cosPsi" )
	, tagName	( "tag" )
	, mistagName	( "mistag" )
	, normalisationCacheValid(false)
	, evaluationCacheValid(false)
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhi_sWave::MakePrototypes()
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
	parameterNames.push_back( As_sqName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_zeroName );
	parameterNames.push_back( delta_sName );
	parameterNames.push_back( deltaMName );
	parameterNames.push_back( Phi_sName );
	allParameters = *( new ParameterSet(parameterNames) );
	valid = true;
}

//Destructor
Bs2JpsiPhi_sWave::~Bs2JpsiPhi_sWave()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2JpsiPhi_sWave::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	bool result = allParameters.SetPhysicsParameters(NewParameterSet);
	//cout << "SetPhysicsParameters: " << result << endl;
	return result;
}

//Calculate the function value
double Bs2JpsiPhi_sWave::Evaluate(DataPoint * measurement)
{
	// The angular functions f1->f6 as defined in roadmap Table 1.
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
	getAngularFunctions( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, measurement );

	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59
	// First for the B
	double AzeroAzeroB, AparaAparaB, AperpAperpB, AsAsB, ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB, ReAsAparaB, ImAsAperpB, ReAsAzeroB;
	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB, ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB, AsAsB, ReAsAparaB, ImAsAperpB, ReAsAzeroB, measurement, 1);

	// Now for the Bbar
	double AzeroAzeroBbar, AparaAparaBbar, AperpAperpBbar, AsAsBbar, ImAparaAperpBbar, ReAzeroAparaBbar, ImAzeroAperpBbar, ReAsAparaBbar, ImAsAperpBbar, ReAsAzeroBbar;
	getTimeDependentAmplitudes( AzeroAzeroBbar, AparaAparaBbar, AperpAperpBbar, ImAparaAperpBbar, ReAzeroAparaBbar, ImAzeroAperpBbar, AsAsBbar, ReAsAparaBbar, ImAsAperpBbar, ReAsAzeroBbar, measurement, -1);

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
		+ f6 * ImAzeroAperpB
		+ f7 * AsAsB
		+ f8 * ReAsAparaB
		+ f9 * ImAsAperpB
		+ f10 * ReAsAzeroB
		;

	//W-
	double v2 = f1 * AzeroAzeroBbar
		+ f2 * AparaAparaBbar
		+ f3 * AperpAperpBbar
		+ f4 * ImAparaAperpBbar
		+ f5 * ReAzeroAparaBbar
		+ f6 * ImAzeroAperpBbar
		+ f7 * AsAsBbar
		+ f8 * ReAsAparaBbar
		+ f9 * ImAsAperpBbar
		+ f10 * ReAsAzeroBbar
		;
	double eval = ( w1*v1 + w2*v2 );
	//cout << "Evaluate: " << eval << endl;
	return eval;
}


double Bs2JpsiPhi_sWave::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
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
		double AzeroAzeroIntB, AparaAparaIntB, AperpAperpIntB, AsAsIntB;
		getTimeAmplitudeIntegrals(  AzeroAzeroIntB
				, AparaAparaIntB
				, AperpAperpIntB
				, AsAsIntB
				, boundary
				, 1);

		double AzeroAzeroIntBbar, AparaAparaIntBbar, AperpAperpIntBbar, AsAsIntBbar;
		getTimeAmplitudeIntegrals(  AzeroAzeroIntBbar
				, AparaAparaIntBbar
				, AperpAperpIntBbar
				, AsAsIntBbar
				, boundary
				, -1);

		cachedv1 = AzeroAzeroIntB + AparaAparaIntB + AperpAperpIntB + AsAsIntB;
		cachedv2 = AzeroAzeroIntBbar + AparaAparaIntBbar + AperpAperpIntBbar + AsAsIntBbar; //FIXME should As terms be here?
		normalisationCacheValid = true;
	}

	double norm = ( w1*cachedv1 + w2*cachedv2 );
	//cout << "Normalisation: " << norm << endl;
	return norm;
}

void Bs2JpsiPhi_sWave::getAngularFunctions( double & f1
		, double & f2
		, double & f3
		, double & f4
		, double & f5
		, double & f6
		, double & f7
		, double & f8
		, double & f9
		, double & f10
		, DataPoint * measurement)
{
	// Observables (the stuff your experiment measures)

	double cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
	double sinTheta  = sqrt(1. - cosTheta*cosTheta);
	double cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();
	double sinPsi	= sqrt(1. - cosPsi*cosPsi);
	double phi	  = measurement->GetObservable( phiName )->GetValue();
	double cosPhi	= cos(phi);
	double sinPhi	 = sin(phi);
	double sin2Theta = 2.*sinTheta*cosTheta;
	double sin2Psi	 = 2.*sinPsi*cosPsi;
	double sin2Phi	 = 2.*sinPhi*cosPhi;
	double cosTheta_sq = cosTheta*cosTheta;
	double sinTheta_sq  = sinTheta*sinTheta;
	double cosPsi_sq   = cosPsi*cosPsi;
	double sinPsi_sq = sinPsi*sinPsi;
	double cosPhi_sq =  cosPhi*cosPhi;
	double sinPhi_sq =  sinPhi*sinPhi;

	double norm = 9./32./TMath::Pi();
	//this is the factor that drops out when you integrate over the angles
	// same factor for f1, f2, f3. The remaining terms f4, f5, f6 give 0.
	//TRANSVERSITY BASIS
	f1 =	2. * cosPsi_sq * ( 1. - sinTheta_sq * cosPhi_sq ) * norm; 		//Check: norm = (9./64./TMath::Pi())
	f2 =	1. * sinPsi_sq * ( 1. - sinTheta_sq * sinPhi_sq ) * norm;		//Check: norm = (9./64./TMath::Pi())
	f3 =	1. * sinPsi_sq * sinTheta_sq * norm;					//Check: norm = (9./64./TMath::Pi())
	f4 =	-1. * sinPsi_sq * sin2Theta * sinPhi * norm;				//Check: norm = 0	//was factor of -1 different from roadmap
	f5 = 	sinTheta_sq * sin2Phi * sin2Psi * norm/sqrt(2.);			//Check: norm = 0	//was factor of -1 different from roadmap
	f6 = 	sin2Theta * cosPhi * sin2Psi * norm/sqrt(2.);				//Check: norm = 0
	f7 = 	(2./3.) * (1. - sinTheta_sq * cosPhi_sq) * norm;			//Check: norm = (9./64./TMath::Pi())
	f8 = 	(1./3.)* sqrt(6.) * sinTheta_sq * sin2Phi * sinPsi * norm;		//Check: norm = 0	//was factor of -1 different from roadmap
	f9 = 	-(1./3.) * sqrt(6.) * sin2Theta * cosPhi* sinPsi * norm;		//Check: norm = 0 	//was factor of -1 different from roadmap
	f10=	-(4./3.) * sqrt(3.) * (1. - sinTheta_sq * cosPhi_sq) * cosPsi * norm; 	//Check: norm = 0	//was factor of -1 different from roadmap

	return;
}

void Bs2JpsiPhi_sWave::getTimeDependentAmplitudes(
		double & AzeroAzero
		, double & AparaApara
		, double & AperpAperp
		, double & ImAparaAperp
		, double & ReAzeroApara
		, double & ImAzeroAperp
		, double & AsAs
		, double & ReAsApara
		, double & ImAsAperp
		, double & ReAsAzero
		, DataPoint * measurement
		, int Btype )
{
	// Observable
	double time = measurement->GetObservable( timeName )->GetValue();

	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	double gamma, deltaGamma, deltaMs, Phi_s;
	double Azero_sq, Apara_sq, Aperp_sq, As_sq;
	double delta_zero, delta_para, delta_perp, delta_s;
	getPhysicsParameters( gamma, deltaGamma, deltaMs, Phi_s, Azero_sq, Apara_sq, Aperp_sq, As_sq, delta_zero, delta_para, delta_perp, delta_s);

	// Quantities depending only on physics parameters can be cached
	if ( !evaluationCacheValid )
	{
		cachedAzero = sqrt( Azero_sq );
		cachedApara = sqrt( Apara_sq );
		cachedAperp = sqrt( Aperp_sq );
		cachedAs = sqrt( As_sq );
		cachedsinDeltaPerpPara	= sin( delta_perp - delta_para );
		cachedcosDeltaPerpPara	= cos( delta_perp - delta_para );
		cachedsinDeltaPerpZero	= sin( delta_perp - delta_zero );
		cachedcosDeltaPerpZero	= cos( delta_perp - delta_zero );
		cachedcosDeltaParaZero	= cos( delta_para - delta_zero );
		cachedsinDeltaPerpS	= sin(delta_perp - delta_s);
		cachedsinDeltaParaS	= sin(delta_para - delta_s);
		cachedsinDeltaZeroS	= sin(delta_zero - delta_s);
		cachedcosDeltaZeroS	= cos(delta_zero - delta_s);
		cachedsinPhis = sin( Phi_s );
		cachedcosPhis = cos( Phi_s );
		evaluationCacheValid = true;
	}

	// Quantities depending on time cannot be cached
	double expGT = exp( -gamma*time );
	double coshDeltaGammaT = cosh( deltaGamma*time/2.);
	double sinhDeltaGammaT = sinh( deltaGamma*time/2.);
	double sinDeltaMsT = Btype * sin( deltaMs*time );
	double cosDeltaMsT = Btype * cos( deltaMs*time );
	double Term_AsAs_AperpAperp, Term_AzeroAzero_AparaApara;
	// Now calculate the amplitudes
	Term_AsAs_AperpAperp = expGT * ( coshDeltaGammaT + cachedcosPhis * sinhDeltaGammaT - cachedsinPhis * sinDeltaMsT );
	AperpAperp = 	Aperp_sq * Term_AsAs_AperpAperp;
	AsAs = 		As_sq *	Term_AsAs_AperpAperp;

	Term_AzeroAzero_AparaApara = expGT * ( coshDeltaGammaT - cachedcosPhis * sinhDeltaGammaT + cachedsinPhis * sinDeltaMsT );
	AzeroAzero = 	Azero_sq * Term_AzeroAzero_AparaApara;
	AparaApara = 	Apara_sq * Term_AzeroAzero_AparaApara;

	ImAparaAperp = cachedApara*cachedAperp * expGT * ( -1.0*cachedcosDeltaPerpPara * cachedsinPhis * sinhDeltaGammaT + cachedsinDeltaPerpPara * cosDeltaMsT - cachedcosDeltaPerpPara * cachedcosPhis * sinDeltaMsT );
	ReAzeroApara = cachedAzero*cachedApara * expGT * cachedcosDeltaParaZero * ( coshDeltaGammaT - cachedcosPhis * sinhDeltaGammaT + cachedsinPhis * sinDeltaMsT );
	ImAzeroAperp = cachedAzero*cachedAperp * expGT * ( -1.0*cachedcosDeltaPerpZero * cachedsinPhis * sinhDeltaGammaT + cachedsinDeltaPerpZero * cosDeltaMsT - cachedcosDeltaPerpZero * cachedcosPhis * sinDeltaMsT );

	ReAsApara = cachedAs*cachedApara * expGT   * ( -1.0*cachedsinDeltaParaS   * cachedsinPhis  * sinhDeltaGammaT + cachedcosDeltaParaS   * cosDeltaMsT - cachedsinDeltaParaS * cachedcosPhis * sinDeltaMsT );
	ImAsAperp = cachedAs*cachedAperp * expGT   * cachedsinDeltaPerpS * (coshDeltaGammaT + cachedcosPhis * sinhDeltaGammaT - cachedsinPhis * sinDeltaMsT );
	ReAsAzero = cachedAs*cachedAzero * expGT   * (-1.0*cachedsinDeltaZeroS * cachedsinPhis *  sinhDeltaGammaT + cachedcosDeltaZeroS * cosDeltaMsT - cachedsinDeltaZeroS * cachedcosPhis * sinDeltaMsT );
	return;
}

void Bs2JpsiPhi_sWave::getTimeAmplitudeIntegrals( double & AzeroAzeroInt
		, double & AparaAparaInt
		, double & AperpAperpInt
		, double & AsAsInt
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
		AsAsInt	= -999.;
		return;
	}
	else
	{
		tlow = timeBound->GetMinimum();
		thigh = timeBound->GetMaximum();
	}

	// Physics parameters
	double gamma, deltaGamma, deltaMs, Phi_s;
	double Azero_sq, Apara_sq, Aperp_sq, As_sq;
	double delta_zero, delta_para, delta_perp, delta_s;
	getPhysicsParameters( gamma, deltaGamma, deltaMs, Phi_s, Azero_sq, Apara_sq, Aperp_sq, As_sq, delta_zero, delta_para, delta_perp, delta_s);

	double G_H	= gamma - 0.5*deltaGamma;
	double G_L	= gamma + 0.5*deltaGamma;

	double tauH   = (1.0 / G_H);
	double tauL   = (1.0 / G_L);
	double tauBar = (1.0 / gamma);

	AzeroAzeroInt = getAzeroAzeroInt( tlow, thigh, Azero_sq, tauL, tauH, tauBar, deltaMs, Phi_s, Btype);
	AparaAparaInt = getAparaAparaInt( tlow, thigh, Apara_sq, tauL, tauH, tauBar, deltaMs, Phi_s, Btype);
	AperpAperpInt = getAperpAperpInt( tlow, thigh, Aperp_sq, tauL, tauH, tauBar, deltaMs, Phi_s, Btype);
	AsAsInt	   = getAsAsInt( tlow, thigh, As_sq,	tauL, tauH, tauBar, deltaMs, Phi_s, Btype);
	// No contribution from interference terms here since they drop out when the angular integration is done.
	return;
}


inline double Bs2JpsiPhi_sWave::getAzeroAzeroInt(double tmin, double tmax,
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


inline double Bs2JpsiPhi_sWave::getAparaAparaInt(double tmin, double tmax,
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


inline double Bs2JpsiPhi_sWave::getAperpAperpInt(double tmin, double tmax,
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

inline double Bs2JpsiPhi_sWave::getAsAsInt(double tmin, double tmax,
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

void Bs2JpsiPhi_sWave::getPhysicsParameters( double & gamma
		, double & deltaGamma
		, double & deltaM
		, double & Phi_s
		, double & Azero_sq
		, double & Apara_sq
		, double & Aperp_sq
		, double & As_sq
		, double & delta_zero
		, double & delta_para
		, double & delta_perp
		, double & delta_s
		)
{
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma	  = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	deltaGamma = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	deltaM	 = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
	Phi_s	  = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();

	Azero_sq   = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
	//Apara_sq   = allParameters.GetPhysicsParameter( Apara_sqName )->GetValue();
	Aperp_sq   = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
	As_sq	   = allParameters.GetPhysicsParameter( As_sqName )->GetValue();

	delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	delta_s = allParameters.GetPhysicsParameter( delta_sName )->GetValue();

	Apara_sq = 1 - Azero_sq - Aperp_sq; //Is this right?

	return;
}
