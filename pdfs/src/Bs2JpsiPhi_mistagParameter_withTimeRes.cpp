// $Id: Bs2JpsiPhi_mistagParameter_withTimeRes.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhi_mistagParameter_withTimeRes Bs2JpsiPhi_mistagParameter_withTimeRes.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi with mistag as a physics parameter
 *  and double gaussian time resolution included.
 *
 *  @author Greig Cowan greig.cowan@cern.ch
 *  @date 2010-01-19
 */

#include "Bs2JpsiPhi_mistagParameter_withTimeRes.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
Bs2JpsiPhi_mistagParameter_withTimeRes::Bs2JpsiPhi_mistagParameter_withTimeRes() : 
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
    	, mistagName	( "mistag" )
        // Time resolution parameters (will be fitted for)
        , timeResName ( "timeRes" )

	// Observables
	, timeName	( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	( "phi" )
	, cosPsiName	( "cosPsi" )
	, tagName	( "tag" )
	, normalisationCacheValid(false)
	, evaluationCacheValid(false)
{
	MakePrototypes();
	
	std::cout << "Constructing J/PsiPhi PDF with mistag as parameter and time resolution." << std::endl;

}

//Make the data point and parameter set
void Bs2JpsiPhi_mistagParameter_withTimeRes::MakePrototypes()
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
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2JpsiPhi_mistagParameter_withTimeRes::~Bs2JpsiPhi_mistagParameter_withTimeRes()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2JpsiPhi_mistagParameter_withTimeRes::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	return allParameters.SetPhysicsParameters(NewParameterSet);
}

//Return a list of parameters not to be integrated
vector<string> Bs2JpsiPhi_mistagParameter_withTimeRes::GetDoNotIntegrateList()
{
	vector<string> list;
	return list;
}

//Calculate the function value
double Bs2JpsiPhi_mistagParameter_withTimeRes::Evaluate(DataPoint * measurement)
{
	double time = measurement->GetObservable( timeName )->GetValue();

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
  	
	double gamma, deltaGamma, deltaMs, Phi_s;
	double Azero_sq, Apara_sq, Aperp_sq;
	double delta_zero, delta_para, delta_perp;
	double omega, timeRes;
	getPhysicsParameters( gamma, deltaGamma, deltaMs, Phi_s, Azero_sq, Apara_sq, Aperp_sq, delta_zero, delta_para, delta_perp, omega, timeRes);
	
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


double Bs2JpsiPhi_mistagParameter_withTimeRes::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	// Now need to know the tag and the mistag
    int q = (int)measurement->GetObservable( tagName )->GetValue();

	double gamma, deltaGamma, deltaMs, Phi_s;
	double Azero_sq, Apara_sq, Aperp_sq;
	double delta_zero, delta_para, delta_perp;
	double omega, timeRes;
	getPhysicsParameters( gamma, deltaGamma, deltaMs, Phi_s, Azero_sq, Apara_sq, Aperp_sq, delta_zero, delta_para, delta_perp, omega, timeRes);
		
	double epsilon[3];
  	epsilon[0] = omega;
  	epsilon[1] = 0.5;
  	epsilon[2] = (1.0 - omega);
  	double w1  = epsilon[q + 1];
  	double w2  = 1.0 - w1;

	if ( !normalisationCacheValid )
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

void Bs2JpsiPhi_mistagParameter_withTimeRes::getAngularFunctions( double & f1
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

void Bs2JpsiPhi_mistagParameter_withTimeRes::getTimeDependentAmplitudes(  double & AzeroAzero
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
	double omega, timeRes;
	getPhysicsParameters( gamma, deltaGamma, deltaMs, Phi_s, Azero_sq, Apara_sq, Aperp_sq, delta_zero, delta_para, delta_perp, omega, timeRes);
	
	// Quantities depending only on physics parameters can be cached
	if ( !evaluationCacheValid )
	{
		cachedAzero = sqrt( Azero_sq );
		cachedApara = sqrt( Apara_sq );
		cachedAperp = sqrt( Aperp_sq );
	
		cachedSinDeltaPerpPara	= sin( delta_perp - delta_para );
		cachedCosDeltaPerpPara	= cos( delta_perp - delta_para );
		cachedSinDeltaPerp	= sin( delta_perp );
		cachedCosDeltaPerp	= cos( delta_perp );
		cachedCosDeltaPara	= cos( delta_para );
	
		cachedSinPhis = sin( Phi_s );
		cachedCosPhis = cos( Phi_s );

		evaluationCacheValid = true;
	}
		
        // We always calculate things for the B first, and these are the same for the Bbar
        if ( Btype == 1 )
        {
		cachedExpCosh = Mathematics::ExpCosh( time, gamma, deltaGamma, timeRes );
		cachedExpSinh = Mathematics::ExpSinh( time, gamma, deltaGamma, timeRes);
		cachedExpCos  = Mathematics::ExpCos(  time, gamma, deltaMs, timeRes );
		cachedExpSin  = Mathematics::ExpSin(  time, gamma, deltaMs, timeRes );
	}

	AzeroAzero = Azero_sq * ( cachedExpCosh - cachedCosPhis * cachedExpSinh + Btype * cachedSinPhis * cachedExpSin ); 
	AparaApara = Apara_sq * ( cachedExpCosh - cachedCosPhis * cachedExpSinh + Btype * cachedSinPhis * cachedExpSin ); 
	AperpAperp = Aperp_sq * ( cachedExpCosh + cachedCosPhis * cachedExpSinh - Btype * cachedSinPhis * cachedExpSin ); 

	ImAparaAperp = cachedApara*cachedAperp * ( - cachedCosDeltaPerpPara * cachedSinPhis * cachedExpSinh 
					       + Btype * cachedSinDeltaPerpPara * cachedExpCos
					       - Btype * cachedCosDeltaPerpPara * cachedCosPhis * cachedExpSin );

	ReAzeroApara = cachedAzero*cachedApara * cachedCosDeltaPara * ( cachedExpCosh - cachedCosPhis * cachedExpSinh
							    + Btype * cachedSinPhis * cachedExpSin );

	ImAzeroAperp = cachedAzero*cachedAperp * ( - cachedCosDeltaPerp * cachedSinPhis * cachedExpSinh
                                               + Btype * cachedSinDeltaPerp * cachedExpCos 
                                               - Btype * cachedCosDeltaPerp * cachedCosPhis * cachedExpSin );
	return;
}

void Bs2JpsiPhi_mistagParameter_withTimeRes::getTimeAmplitudeIntegrals( double & AzeroAzeroInt
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
		if (tlow < 0 ) tlow = 0;  // ADDED BY PELC
		thigh = timeBound->GetMaximum();
	}
	
	// Physics parameters
	double gamma, deltaGamma, deltaMs, Phi_s;
	double Azero_sq, Apara_sq, Aperp_sq;
	double delta_zero, delta_para, delta_perp;
	double omega, timeRes;
	getPhysicsParameters( gamma, deltaGamma, deltaMs, Phi_s, Azero_sq, Apara_sq, Aperp_sq, delta_zero, delta_para, delta_perp, omega, timeRes);
	
	double cosPhis = cos(Phi_s);
	double sinPhis = sin(Phi_s);
	double expCoshInt = Mathematics::ExpCoshInt( tlow, thigh, gamma, deltaGamma, timeRes );
	double expSinhInt = Mathematics::ExpSinhInt( tlow, thigh, gamma, deltaGamma, timeRes );
	double expSinInt  = Mathematics::ExpSinInt(  tlow, thigh, gamma, deltaMs, timeRes );
	
        AzeroAzeroInt = Azero_sq * ( expCoshInt - cosPhis * expSinhInt + Btype * sinPhis * expSinInt );
        AparaAparaInt = Apara_sq * ( expCoshInt - cosPhis * expSinhInt + Btype * sinPhis * expSinInt );
        AperpAperpInt = Aperp_sq * ( expCoshInt + cosPhis * expSinhInt - Btype * sinPhis * expSinInt );
	// No contribution from interference terms here since they drop out when the angular integration is done.
	return;
}

void Bs2JpsiPhi_mistagParameter_withTimeRes::getPhysicsParameters( double & gamma
					, double & deltaGamma
					, double & deltaM
					, double & Phi_s
					, double & Azero_sq
					, double & Apara_sq
					, double & Aperp_sq
					, double & delta_zero 
					, double & delta_para
					, double & delta_perp
					, double & omega  
					, double & timeRes  )
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
	omega = allParameters.GetPhysicsParameter( mistagName )->GetValue();
	timeRes = allParameters.GetPhysicsParameter( timeResName )->GetValue();

	Apara_sq = 1 - Azero_sq - Aperp_sq;

	return;
}
