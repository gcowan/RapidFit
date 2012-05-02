// $Id: Bs2JpsiPhi_mistagObservable_withTimeRes.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhi_mistagObservable_withTimeRes Bs2JpsiPhi_mistagObservable_withTimeRes.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi including analytic time resolution
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-02
 */

#include "Bs2JpsiPhi_mistagObservable_withTimeRes.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"

//Constructor
Bs2JpsiPhi_mistagObservable_withTimeRes::Bs2JpsiPhi_mistagObservable_withTimeRes() : 
	// Physics parameters (will be fitted for)
	  gammaName     ( "gamma" )
	, deltaGammaName( "deltaGamma" )
	, deltaMName    ( "deltaM")
	, Phi_sName     ( "Phi_s")
	, Azero_sqName  ( "Azero_sq" )
	, Aperp_sqName  ( "Aperp_sq" )
	, delta_zeroName( "delta_zero" )
	, delta_paraName( "delta_para" )
	, delta_perpName( "delta_perp" )
        // Time resolution parameters (will be fitted for)
	, mean_time_res1Name ( "mean_time_res1" )
	, mean_time_res2Name ( "mean_time_res2" )
	, sigma_time_res1Name( "sigma_time_res1" )
	, sigma_time_res2Name( "sigma_time_res2" )
	, frac_time_res1Name ( "frac_time_res1" )
	// Observables (input to the fit)
	, timeName	( "time" )
	, timeErrName	( "timeErr" )
	, cosThetaName	( "cosTheta" )
	, phiName	( "phi" )
	, cosPsiName	( "cosPsi" )
	// Tagging (input to the fit, these could also be fit for)
	, tagName	( "tag" )
	, mistagName	( "mistag" )
	
	, normalisationCacheValid(false)
	, evaluationCacheValid(false)
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhi_mistagObservable_withTimeRes::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( timeErrName );
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
	parameterNames.push_back( mean_time_res1Name );
	parameterNames.push_back( mean_time_res2Name );
	parameterNames.push_back( sigma_time_res1Name );
	parameterNames.push_back( sigma_time_res2Name );
	parameterNames.push_back( frac_time_res1Name );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2JpsiPhi_mistagObservable_withTimeRes::~Bs2JpsiPhi_mistagObservable_withTimeRes()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2JpsiPhi_mistagObservable_withTimeRes::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	return allParameters.SetPhysicsParameters(NewParameterSet);
}

//Return a list of parameters not to be integrated
vector<string> Bs2JpsiPhi_mistagObservable_withTimeRes::GetDoNotIntegrateList()
{
	vector<string> mistagList;
	mistagList.push_back(mistagName);
	mistagList.push_back(timeErrName);
	return mistagList;
}

//Calculate the function value
double Bs2JpsiPhi_mistagObservable_withTimeRes::Evaluate(DataPoint * measurement)
{
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

	//cout << "angles " << f1 << " " << f2 << " " << f3 << " " << f4 << " " << f5 << " " << f6 << endl;
	
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
 
	//cout << "Evaluate " << w1*v1 + w2*v2 << endl; 
  	return ( w1*v1 + w2*v2 );
}


double Bs2JpsiPhi_mistagObservable_withTimeRes::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	// Now need to know the tag and the mistag
    	int q = (int)measurement->GetObservable( tagName )->GetValue();
    	double omega = measurement->GetObservable( mistagName )->GetValue();	
    	double timeErr = measurement->GetObservable( timeErrName )->GetValue();	

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
		double AparaAperpIntB, AzeroAparaIntB, AzeroAperpIntB;
		getTimeAmplitudeIntegrals(  AzeroAzeroIntB
							  , AparaAparaIntB
							  , AperpAperpIntB
							  , AparaAperpIntB
							  , AzeroAparaIntB
							  , AzeroAperpIntB
							  , boundary
							  , 1
							  , timeErr);

		double AzeroAzeroIntBbar, AparaAparaIntBbar, AperpAperpIntBbar;
		double AparaAperpIntBbar, AzeroAparaIntBbar, AzeroAperpIntBbar;
		getTimeAmplitudeIntegrals(  AzeroAzeroIntBbar
							  , AparaAparaIntBbar
							  , AperpAperpIntBbar
							  , AparaAperpIntBbar
							  , AzeroAparaIntBbar
							  , AzeroAperpIntBbar
							  , boundary
							  , -1
							  , timeErr);

		//cout << AzeroAzeroIntB << " " << AparaAparaIntB << " " << AperpAperpIntB << endl;	
		//cout << AzeroAzeroIntBbar << " " << AparaAparaIntBbar << " " << AperpAperpIntBbar << endl;	
		cachedv1 = AzeroAzeroIntB + AparaAparaIntB + AperpAperpIntB 
			 + 0.*(AparaAperpIntB + AzeroAparaIntB + AzeroAperpIntB);
		cachedv2 = AzeroAzeroIntBbar + AparaAparaIntBbar + AperpAperpIntBbar;
			 + 0.*(AparaAperpIntBbar + AzeroAparaIntBbar + AzeroAperpIntBbar);

		normalisationCacheValid = true;
	}
	//cout << "Normalisation\t" << w1*cachedv1 + w2*cachedv2 << endl;	
	return ( w1*cachedv1 + w2*cachedv2 );
}

void Bs2JpsiPhi_mistagObservable_withTimeRes::getAngularFunctions( double & f1
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
	//this is  1/factor that drops out when you integrate over the angles
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

void Bs2JpsiPhi_mistagObservable_withTimeRes::getTimeDependentAmplitudes(  double & AzeroAzero
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
	double timeErr = measurement->GetObservable( timeErrName )->GetValue();
        
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	double gamma, deltaGamma, deltaMs, Phi_s;
	double Azero_sq, Apara_sq, Aperp_sq;
	double delta_zero, delta_para, delta_perp;
	double mean_time_res1, sigma_time_res1;
	double mean_time_res2, sigma_time_res2;
	double frac_time_res1;
	getPhysicsParameters( gamma, deltaGamma, deltaMs, Phi_s
			    , Azero_sq, Apara_sq, Aperp_sq
			    , delta_zero, delta_para, delta_perp
			    , mean_time_res1, sigma_time_res1
			    , mean_time_res2, sigma_time_res2
			    , frac_time_res1);
	
	//cout << "inputs; gamma " << gamma << " deltaMs " << deltaMs << " " << mean_time_res1 << " " << mean_time_res2 << " " << sigma_time_res1 << " " << sigma_time_res2 << " frac_time_res " << frac_time_res1<< " " <<  time << " " << timeErr << endl;

	double gammaL =  gamma + deltaGamma/2.;
  	double gammaH =  gamma - deltaGamma/2.;

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

	if (isnan(deltaMs)){
		cout << "inputs; gamma " << gamma << " deltaGamma " << deltaGamma << " deltaMs " << deltaMs << " " << mean_time_res1 << " " << mean_time_res2 << " " << sigma_time_res1 << " " << sigma_time_res2 << " frac_time_res " << frac_time_res1<< " " <<  time << " " << timeErr << endl;
	}

	// We always calculate things for the B first, so this should be OK
	if ( Btype == 1 )
	{
		// These quantities don't change depending on the B type
        	cachedExpCosh = ExpCosh( gammaL, gammaH
				, mean_time_res1*timeErr, sigma_time_res1*timeErr
                            	, mean_time_res2*timeErr, sigma_time_res2*timeErr
                            	, frac_time_res1, time );
        	cachedExpSinh = ExpSinh( gammaL, gammaH
				, mean_time_res1*timeErr, sigma_time_res1*timeErr
                            	, mean_time_res2*timeErr, sigma_time_res2*timeErr
                            	, frac_time_res1, time );
        	cachedExpCos = ExpCos( gamma, deltaMs
                              , mean_time_res1*timeErr, sigma_time_res1*timeErr
                              , mean_time_res2*timeErr, sigma_time_res2*timeErr
                              , frac_time_res1, time );
        	cachedExpSin = ExpSin( gamma, deltaMs
                              , mean_time_res1*timeErr, sigma_time_res1*timeErr
                              , mean_time_res2*timeErr, sigma_time_res2*timeErr
                              , frac_time_res1, time ); 
		if (isnan(cachedExpCosh) || isnan(cachedExpSinh) || isnan(cachedExpCos) || isnan(cachedExpSin)){
		cout << "inputs; gamma " << gamma << " deltaGamma " << deltaGamma << " deltaMs " << deltaMs << " " << mean_time_res1 << " " << mean_time_res2 << " " << sigma_time_res1 << " " << sigma_time_res2 << " frac_time_res " << frac_time_res1<< " " <<  time << " " << timeErr << endl;
		cout << "ExpTrig " << cachedExpCosh << " " << cachedExpSinh << " " << cachedExpCos << " " << cachedExpSin << endl;
		exit(1);
		}
	}

	AzeroAzero = Azero_sq * ( cachedExpCosh - cachedcosPhis * cachedExpSinh + Btype * cachedsinPhis * cachedExpSin ); 
	AparaApara = Apara_sq * ( cachedExpCosh - cachedcosPhis * cachedExpSinh + Btype * cachedsinPhis * cachedExpSin ); 
	AperpAperp = Aperp_sq * ( cachedExpCosh + cachedcosPhis * cachedExpSinh - Btype * cachedsinPhis * cachedExpSin ); 
	        
	ImAparaAperp = cachedApara*cachedAperp * ( - cachedcosDeltaPerpPara * cachedsinPhis * cachedExpSinh
                                               + Btype * cachedsinDeltaPerpPara * cachedExpCos
                                               - Btype * cachedcosDeltaPerpPara * cachedcosPhis * cachedExpSin );
	       	
	ReAzeroApara = cachedAzero*cachedApara * cachedcosDeltaPara * ( cachedExpCosh - cachedcosPhis * cachedExpSinh
                                               + Btype * cachedsinPhis * cachedExpSin );
        	
	ImAzeroAperp = cachedAzero*cachedAperp * ( - cachedcosDeltaPerp * cachedsinPhis * cachedExpSinh
                                               + Btype * cachedsinDeltaPerp * cachedExpCos
                                               - Btype * cachedcosDeltaPerp * cachedcosPhis * cachedExpSin );
	//cout << "Amps " << AzeroAzero << " " << AparaApara << " " << AperpAperp << " " << ImAparaAperp << " " << ReAzeroApara << " " << ImAzeroAperp << endl;
	return;
}

void Bs2JpsiPhi_mistagObservable_withTimeRes::getTimeAmplitudeIntegrals( double & AzeroAzeroInt
                                                 , double & AparaAparaInt
                                                 , double & AperpAperpInt
                                                 , double & AparaAperpInt
                                                 , double & AzeroAparaInt
                                                 , double & AzeroAperpInt
						 , PhaseSpaceBoundary * boundary
						 , int Btype
						 , double timeErr)
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
		thigh = timeBound->GetMaximum();
	}
	
	// Physics parameters
        double gamma, deltaGamma, deltaMs, Phi_s;
        double Azero_sq, Apara_sq, Aperp_sq;
        double delta_zero, delta_para, delta_perp;
        double mean_time_res1, sigma_time_res1;
        double mean_time_res2, sigma_time_res2;
        double frac_time_res1;
        getPhysicsParameters( gamma, deltaGamma, deltaMs, Phi_s
                            , Azero_sq, Apara_sq, Aperp_sq
                            , delta_zero, delta_para, delta_perp
                            , mean_time_res1, sigma_time_res1
                            , mean_time_res2, sigma_time_res2
        		    , frac_time_res1);
	
	double AparaAperp = sqrt(Apara_sq)*sqrt(Aperp_sq);
	double AzeroApara = sqrt(Azero_sq)*sqrt(Apara_sq);
	double AzeroAperp = sqrt(Azero_sq)*sqrt(Aperp_sq);

	double gammaH = gamma - 0.5*deltaGamma;
	double gammaL = gamma + 0.5*deltaGamma;
	
	double tauH   = 1. / gammaH;
	double tauL   = 1. / gammaL;
	///double tauBar = 1. / gamma;

	//double timeErr = 1.0; // just a placeholder for now
	double expCosh_integral, expSinh_integral;
	double expCos_integral, expSin_integral;
	getExpTrigIntegrals( expCosh_integral, expSinh_integral
			   , expCos_integral, expSin_integral
			   , tlow, thigh
			   , tauL, tauH, deltaMs
			   , mean_time_res1 * timeErr
			   , mean_time_res2 * timeErr
			   , sigma_time_res1 * timeErr
			   , sigma_time_res2 * timeErr
			   , frac_time_res1);

	// Now need to calculate these terms with propertime resolution included. See Yuehong's code.
	AzeroAzeroInt = getAzeroAzeroInt( Azero_sq, expCosh_integral, expSinh_integral, expSin_integral, Phi_s, Btype);
	AparaAparaInt = getAparaAparaInt( Apara_sq, expCosh_integral, expSinh_integral, expSin_integral, Phi_s, Btype);
	AperpAperpInt = getAperpAperpInt( Aperp_sq, expCosh_integral, expSinh_integral, expSin_integral, Phi_s, Btype);
	//AparaAperpInt = getAparaAperpInt( AparaAperp, expSinh_integral, expCos_integral, expSin_integral, Phi_s, Btype, delta_perp, delta_para);
	//AzeroAparaInt = getAzeroAparaInt( AzeroApara, expCosh_integral, expSinh_integral, expSin_integral, Phi_s, Btype, delta_para);
	//AzeroAperpInt = getAzeroAperpInt( AzeroAperp, expSinh_integral, expCos_integral, expSin_integral, Phi_s, Btype, delta_perp);
	AparaAperpInt = 0.0; // since these terms have no contribution since angular terms int to 0
	AzeroAparaInt = 0.0;
	AzeroAperpInt = 0.0;
	return;
}

inline double Bs2JpsiPhi_mistagObservable_withTimeRes::getAzeroAzeroInt( double AzeroAzero
							    , double expCosh_integral
							    , double expSinh_integral
							    , double expSin_integral
							    , double phis
							    , int Btype)
{
        double v = expCosh_integral - cos(phis) * expSinh_integral 
		 + sin(phis) * expSin_integral * Btype;

        return AzeroAzero * v;
}

inline double Bs2JpsiPhi_mistagObservable_withTimeRes::getAparaAparaInt( double AparaApara
                                                            , double expCosh_integral
                                                            , double expSinh_integral
                                                            , double expSin_integral
                                                            , double phis
                                                            , int Btype)
{
        double v = expCosh_integral - cos(phis) * expSinh_integral
        	 + sin(phis)*expSin_integral * Btype;

        return AparaApara * v;
}

inline double Bs2JpsiPhi_mistagObservable_withTimeRes::getAperpAperpInt( double AperpAperp
                                                            , double expCosh_integral
                                                            , double expSinh_integral
                                                            , double expSin_integral
                                                            , double phis
                                                            , int Btype)
{       
        double v = expCosh_integral + cos(phis) * expSinh_integral
                 + sin(phis)*expSin_integral * Btype;

        return AperpAperp * v;
}

inline double Bs2JpsiPhi_mistagObservable_withTimeRes::getAparaAperpInt( double AparaAperp
                                                            , double expSinh_integral
                                                            , double expCos_integral
                                                            , double expSin_integral
                                                            , double phis
                                                            , int Btype
							    , double delta_perp
							    , double delta_para)
{
        double v = -cos( delta_perp - delta_para ) * sin(phis) * expSinh_integral
        	 +  sin( delta_perp - delta_para ) * expCos_integral * Btype;
        	 -  cos( delta_perp - delta_para ) * cos(phis) * expSin_integral * Btype;
	
	return AparaAperp * v;
}

inline double Bs2JpsiPhi_mistagObservable_withTimeRes::getAzeroAparaInt( double AzeroApara
                                                            , double expCosh_integral
                                                            , double expSinh_integral
                                                            , double expSin_integral
                                                            , double phis
                                                            , int Btype
							    , double delta_para)
{
        double v = ( expCosh_integral - cos(phis) * expSinh_integral
                   + sin(phis)*expSin_integral * Btype ) * cos(delta_para);

        return AzeroApara * v;
}

inline double Bs2JpsiPhi_mistagObservable_withTimeRes::getAzeroAperpInt( double AzeroAperp
                                                            , double expSinh_integral
                                                            , double expCos_integral
                                                            , double expSin_integral
                                                            , double phis
                                                            , int Btype
							    , double delta_perp)
{
        double v = -cos(delta_perp) * sin(phis) * expSinh_integral
        	 +  sin(delta_perp) * expCos_integral * Btype
        	 -  cos(delta_perp) * cos(phis) * expSin_integral * Btype;

        return AzeroAperp * v;
}


inline void Bs2JpsiPhi_mistagObservable_withTimeRes::getExpTrigIntegrals(
	  double & expCosh_integral
	, double & expSinh_integral
	, double & expCos_integral
	, double & expSin_integral
	, double t_low
	, double t_high
	, double tauL
	, double tauH
	, double deltaM
	, double mean_time_res1
	, double mean_time_res2
	, double sigma_time_res1
	, double sigma_time_res2
	, double frac_time_res1
)	
{
	double deltaCosh1 = 0.;
	double deltaSinh1 = 0.;
	double deltaCos1 = 0.;
	double deltaSin1 = 0.;
	getErfPart( t_low, t_high, mean_time_res1, sigma_time_res1, tauL, tauH
		  , deltaM, deltaCosh1, deltaSinh1, deltaCos1, deltaSin1); 
	
	double deltaCosh2 = 0.;
	double deltaSinh2 = 0.;
	double deltaCos2 = 0.;
	double deltaSin2 = 0.;
	getErfPart( t_low, t_high, mean_time_res2, sigma_time_res2, tauL, tauH
		  , deltaM, deltaCosh2, deltaSinh2, deltaCos2, deltaSin2); 
	
	expCosh_integral  = deltaCosh1 * frac_time_res1;
	expCosh_integral += deltaCosh2 * (1. - frac_time_res1);
	expSinh_integral  = deltaSinh1 * frac_time_res1;
	expSinh_integral += deltaSinh2 * (1. - frac_time_res1);
	expCos_integral   = deltaCos1  * frac_time_res1;
	expCos_integral  += deltaCos2  * (1. - frac_time_res1);
	expSin_integral   = deltaSin1  * frac_time_res1;
	expSin_integral  += deltaSin2  * (1. - frac_time_res1);

	return;
}

inline void Bs2JpsiPhi_mistagObservable_withTimeRes::getErfPart(   double t_low
						      , double t_high
						      , double mean
						      , double sigma
						      , double tau_L
						      , double tau_H
						      , double deltaM
						      , double & deltaCosh
						      , double & deltaSinh
						      , double & deltaCos
						      , double & deltaSin )
							
{
	double tau = 2. / ( 1./tau_L + 1./tau_H);
	double deltaExpH = 0.;
        double deltaExpL = 0.;

        static double root2(sqrt(2.));

        double cH = sigma/(root2 * tau_H);
        double xpminH = (t_low - mean)/tau_H;
        double xpmaxH = (t_high - mean)/tau_H;
        double cL = sigma/(root2 * tau_L);
        double xpminL = (t_low - mean)/tau_L;
        double xpmaxL = (t_high - mean)/tau_L;

        double umin = xpminH/(2. * cH);
        double umax = xpmaxH/(2. * cH);

        deltaExpH += -1. * tau_H * ( RooMath::erf(-umax) - RooMath::erf(-umin) +
                                     exp(cH*cH) * ( exp(-xpmaxH) * RooMath::erfc(-umax + cH)
                                                  - exp(-xpminH) * RooMath::erfc(-umin + cH) ));
        deltaExpH /= 2.;

        deltaExpL += -1. * tau_L * ( RooMath::erf(-umax) - RooMath::erf(-umin) +
                                     exp(cL*cL) * ( exp(-xpmaxL) * RooMath::erfc(-umax + cL)
                                                  - exp(-xpminL) * RooMath::erfc(-umin + cL) ));
        deltaExpL /= 2.;

        double c = sigma/(root2*tau);
        double wt = deltaM * tau;
        RooComplex evalDif(evalCerf(-wt, -umax, c) - evalCerf(-wt, -umin, c));

        deltaSin = -tau/(1. + wt*wt) * ( -evalDif.im() + wt*evalDif.re() + wt*(RooMath::erf(-umax) - RooMath::erf(-umin)) ) ;
        deltaCos = -tau/(1. + wt*wt) * (  evalDif.re() + wt*evalDif.im() + RooMath::erf(-umax) - RooMath::erf(-umin) );
        deltaSin /= 2.;
        deltaCos /= 2.;

        deltaCosh = (deltaExpH + deltaExpL)/2.;
        deltaSinh = (deltaExpH - deltaExpL)/2.;

	return;
}

void Bs2JpsiPhi_mistagObservable_withTimeRes::getPhysicsParameters( double & gamma
					, double & deltaGamma
					, double & deltaM
					, double & Phi_s
					, double & Azero_sq
					, double & Apara_sq
					, double & Aperp_sq
					, double & delta_zero 
					, double & delta_para
					, double & delta_perp
					, double & mean_time_res1
					, double & sigma_time_res1
                            		, double & mean_time_res2
					, double & sigma_time_res2
                            		, double & frac_time_res1)
{
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
        deltaGamma = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	deltaM     = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
	Phi_s      = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
	Azero_sq   = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
	Apara_sq   = 1 - Azero_sq - Aperp_sq;
	Aperp_sq   = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
	delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	mean_time_res1  = allParameters.GetPhysicsParameter( mean_time_res1Name )->GetValue();
	sigma_time_res1 = allParameters.GetPhysicsParameter( sigma_time_res1Name )->GetValue();
	mean_time_res2  = allParameters.GetPhysicsParameter( mean_time_res2Name )->GetValue();
	sigma_time_res2 = allParameters.GetPhysicsParameter( sigma_time_res2Name )->GetValue();
	frac_time_res1  = allParameters.GetPhysicsParameter( frac_time_res1Name )->GetValue();
	if (isnan(deltaM)){
		cout << "oh no, nan" << endl;
	}
	return;
}

double Bs2JpsiPhi_mistagObservable_withTimeRes::ExpCosh( double GammaL, double GammaH
              , double mean_time_res1, double sigma_time_res1
              , double mean_time_res2, double sigma_time_res2
              , double frac_time_res1, double time)
{
        double v = 0.;

        // Calculate the contribution from each resolution term
        v += erfc(GammaL, GammaH, mean_time_res1, sigma_time_res1, time, 1)*frac_time_res1;
        v += erfc(GammaL, GammaH, mean_time_res2, sigma_time_res2, time, 1)*(1-frac_time_res1);

        return v;
}


double Bs2JpsiPhi_mistagObservable_withTimeRes::ExpSinh( double GammaL, double GammaH
	      , double mean_time_res1, double sigma_time_res1
              , double mean_time_res2, double sigma_time_res2
              , double frac_time_res1, double time) 
{
  	double v = 0.;
   
	// Calculate the contribution from each resolution term
      	v += erfc(GammaL, GammaH, mean_time_res1, sigma_time_res1, time, -1)*frac_time_res1;
      	v += erfc(GammaL, GammaH, mean_time_res2, sigma_time_res2, time, -1)*(1-frac_time_res1);

        return v;
}

double Bs2JpsiPhi_mistagObservable_withTimeRes::erfc( double gammaL, double gammaH
					 , double mean, double sigma
					 , double time, int coshOrSinh )
{
        double xH = gammaH * ( time - mean );
        double cH = gammaH * sigma/sqrt(2.);
        double xL = gammaL * ( time - mean );
        double cL = gammaL * sigma/sqrt(2.);
        double u  = ( time - mean )/( sqrt(2.) * sigma );
        
	double delta =              exp( cH * cH - xH) * TMath::Erfc( cH - u )/2.
        	     + coshOrSinh * exp( cL * cL - xL) * TMath::Erfc( cL - u )/2.;
        delta /= 2.;

	return delta;
}

double Bs2JpsiPhi_mistagObservable_withTimeRes::ExpCos( double gamma, double deltaM
              , double mean_time_res1, double sigma_time_res1
              , double mean_time_res2, double sigma_time_res2
              , double frac_time_res1, double time)
{
	double v = 0.;
        v += cerfRe(gamma, deltaM, mean_time_res1, sigma_time_res1, time)*frac_time_res1;
        v += cerfRe(gamma, deltaM, mean_time_res2, sigma_time_res2, time)*(1-frac_time_res1);

  	return v;
}

double Bs2JpsiPhi_mistagObservable_withTimeRes::ExpSin( double gamma, double deltaM
              , double mean_time_res1, double sigma_time_res1
              , double mean_time_res2, double sigma_time_res2
              , double frac_time_res1, double time)
{
        double v = 0.;
        v += cerfIm(gamma, deltaM, mean_time_res1, sigma_time_res1, time)*frac_time_res1;
        v += cerfIm(gamma, deltaM, mean_time_res2, sigma_time_res2, time)*(1-frac_time_res1);
        
        return v;
}

double Bs2JpsiPhi_mistagObservable_withTimeRes::cerfRe( double gamma, double deltaM
                                           , double mean, double sigma
                                           , double time )
{
        double c  = gamma*sigma/sqrt(2.);
        double u  = (time - mean)/( sqrt(2.) * sigma );
        double wt = deltaM/gamma;
        double delta = evalCerfRe( -wt, -u, c)/2;

	return delta;
}

double Bs2JpsiPhi_mistagObservable_withTimeRes::cerfIm( double gamma, double deltaM
                                           , double mean, double sigma
                                           , double time)
{
        double c  = gamma*sigma/sqrt(2.);
        double u  = (time - mean)/( sqrt(2.) * sigma );
        double wt = deltaM/gamma;
        double delta = -1.*evalCerfRe( -wt, -u, c)/2;

        return delta;
}

// These functions are taken from Yuehong's code. You have been warned...
RooComplex Bs2JpsiPhi_mistagObservable_withTimeRes::evalCerfApprox( double swt, double u, double c )
{
  static double rootpi = sqrt( atan2(0.,-1.) );
  RooComplex z( swt * c, u + c);
  RooComplex zc( u + c , -swt * c);
  RooComplex zsq = z * z;
  RooComplex v = -zsq - u * u;

  return v.exp() * ( -zsq.exp() / (zc * rootpi) + 1) * 2;
}

double Bs2JpsiPhi_mistagObservable_withTimeRes::evalCerfRe( double swt, double u, double c)
{
    RooComplex z(swt * c, u + c);
    return (z.im() > -4.0) ? RooMath::FastComplexErrFuncRe(z) * exp( -u * u ) : evalCerfApprox( swt, u, c).re();
}

double Bs2JpsiPhi_mistagObservable_withTimeRes::evalCerfIm( double swt, double u, double c)  {
    RooComplex z( swt * c, u + c);
    return (z.im() > -4.0) ? RooMath::FastComplexErrFuncIm(z) * exp( -u * u ) : evalCerfApprox( swt, u, c).im();
}

RooComplex Bs2JpsiPhi_mistagObservable_withTimeRes::evalCerf( double swt, double u, double c ) {
    RooComplex z(swt * c, u + c);
    return (z.im() > -4.0) ? RooMath::FastComplexErrFunc(z) * exp( -u * u ) : evalCerfApprox( swt, u, c);
}


