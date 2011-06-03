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
Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc() : cachedAzeroAzeroIntB(), cachedAparaAparaIntB(),
cachedAperpAperpIntB(), cachedAparaAperpIntB(), cachedAzeroAparaIntB(),
cachedAzeroAperpIntB(), cachedAzeroAzeroIntBbar(), cachedAparaAparaIntBbar(),
cachedAperpAperpIntBbar(), cachedAparaAperpIntBbar(), cachedAzeroAparaIntBbar(),
cachedAzeroAperpIntBbar(), cachedAzero(), cachedApara(), cachedAperp(), cachedExpCosh(),
cachedExpSinh(), cachedExpCos(), cachedExpSin(), cachedSinDeltaPerpPara(),
cachedSinDeltaPerp(), cachedSinDeltaPara(), cachedCosDeltaPerpPara(),
cachedCosDeltaPerp(), cachedCosDeltaPara(), cachedSinPhis(), cachedCosPhis(),

	normalisationCacheValid(false)
	, evaluationCacheValid(false)

	// Physics parameters
	, gammaName     ( "gamma" )
	, deltaGammaName( "deltaGamma" )
	, deltaMName    ( "deltaM")
	, Phi_sName     ( "Phi_s")
	, Azero_sqName  ( "Azero_sq" )
	, Apara_sqName  ( "Apara_sq" )
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
	, timeRes1Name	( "timeResolution1" )
	, timeRes2Name	( "timeResolution2" )
	, timeRes1FractionName	( "timeResolution1Fraction" )

	// Observables
	, timeName	( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	( "phi" )
	, cosPsiName	( "cosPsi" )
	//, timeres	( "resolution" )
	, tagName	( "tag" )

	, gamma(), deltaGamma(), deltaMs(), Phi_s(), Azero_sq(), Apara_sq(), Aperp_sq(),
	AzeroApara(), AzeroAperp(), AparaAperp(), delta_zero(), delta_para(),
	delta_perp(), omega(), timeRes(), timeRes1(), timeRes2(), timeRes1Frac(),
	angAccI1(), angAccI2(), angAccI3(), angAccI4(), angAccI5(), angAccI6(), time(),
	cosTheta(), phi(), cosPsi(), q(), tlow(), thigh()
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
	parameterNames.push_back( timeRes1Name );
	parameterNames.push_back( timeRes2Name );
	parameterNames.push_back( timeRes1FractionName );
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
        omega    = allParameters.GetPhysicsParameter( mistagName )->GetValue();
        timeRes1 = allParameters.GetPhysicsParameter( timeRes1Name )->GetValue();
        timeRes2 = allParameters.GetPhysicsParameter( timeRes2Name )->GetValue();
        timeRes1Frac = allParameters.GetPhysicsParameter( timeRes1FractionName )->GetValue();
        angAccI1 = allParameters.GetPhysicsParameter( angAccI1Name )->GetValue();
        angAccI2 = allParameters.GetPhysicsParameter( angAccI2Name )->GetValue();
        angAccI3 = allParameters.GetPhysicsParameter( angAccI3Name )->GetValue();
        angAccI4 = allParameters.GetPhysicsParameter( angAccI4Name )->GetValue();
        angAccI5 = allParameters.GetPhysicsParameter( angAccI5Name )->GetValue();
        angAccI6 = allParameters.GetPhysicsParameter( angAccI6Name )->GetValue();

        Apara_sq = 1 - Azero_sq - Aperp_sq;
	if ( Apara_sq < 0.) return false;
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
	time = measurement->GetObservable( timeName )->GetValue();
        cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
        phi      = measurement->GetObservable( phiName )->GetValue();
        cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();
	q = (int)measurement->GetObservable( tagName )->GetValue();

        if(timeRes1Frac >= 0.9999)
	{
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                return buildPDFnumerator();
        }
        else
	{
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                double val1 = buildPDFnumerator();
                // Set the member variable for time resolution to the second value and calculate
                timeRes = timeRes2;
                double val2 = buildPDFnumerator();
                return timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
        }
}

double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::buildPDFnumerator()
{
	// The angular functions f1->f6 as defined in roadmap Table 1.
	double f1, f2, f3, f4, f5, f6;
	Mathematics::getBs2JpsiPhiAngularFunctions( f1, f2, f3, f4, f5, f6, cosTheta, phi, cosPsi );

	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59
	// First for the B
	double AzeroAzeroB, AparaAparaB, AperpAperpB;
	double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB;
	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB
			, ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
			, 1);

	// Now for the Bbar
	double AzeroAzeroBbar, AparaAparaBbar, AperpAperpBbar;
	double ImAparaAperpBbar, ReAzeroAparaBbar, ImAzeroAperpBbar;
	getTimeDependentAmplitudes( AzeroAzeroBbar, AparaAparaBbar, AperpAperpBbar
			, ImAparaAperpBbar, ReAzeroAparaBbar, ImAzeroAperpBbar
			, -1);

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
	q = (int)measurement->GetObservable( tagName )->GetValue();

	IConstraint * timeBound = boundary->GetConstraint("time");
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


        if(timeRes1Frac >= 0.9999)
        {
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                return buildPDFdenominator();
        }
        else
        {
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                double val1 = buildPDFdenominator();
		normalisationCacheValid = false;
                // Set the member variable for time resolution to the second value and calculate
                timeRes = timeRes2;
                double val2 = buildPDFdenominator();
                return timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
        }
}

double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::buildPDFdenominator()
{
	// Work out the mistag
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
				, 1);

		getTimeAmplitudeIntegrals(  cachedAzeroAzeroIntBbar
				, cachedAparaAparaIntBbar
				, cachedAperpAperpIntBbar
				, cachedAparaAperpIntBbar
				, cachedAzeroAparaIntBbar
				, cachedAzeroAperpIntBbar
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

	//double norm = (angAccI1 + angAccI2 + angAccI3 + angAccI4 + angAccI5 + angAccI6)/(3*32*TMath::Pi()/9);
	double norm = 1.0;
	return (w1*v1 + w2*v2)/norm;
}

void Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::getTimeDependentAmplitudes(  double & AzeroAzero
		, double & AparaApara
		, double & AperpAperp
		, double & ImAparaAperp
		, double & ReAzeroApara
		, double & ImAzeroAperp
		, int Btype )
{
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

	//cout << gamma << " " << deltaGamma << " " << Azero_sq << " " << Aperp_sq << endl;
	//cout << cachedExpCosh << " " << cachedExpSinh << " " << cachedExpCos << " " << cachedExpSin << endl;
	// Now calculate the amplitudes
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

void Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc::getTimeAmplitudeIntegrals( double & AzeroAzeroInt
		, double & AparaAparaInt
		, double & AperpAperpInt
		, double & AparaAperpInt
		, double & AzeroAparaInt
		, double & AzeroAperpInt
		, int Btype)
{

        double cosPhis = cos(Phi_s);
        double sinPhis = sin(Phi_s);
	double cosDeltaPerpMinusPara = cos(delta_perp - delta_para);
        double sinDeltaPerpMinusPara = sin(delta_perp - delta_para);

        double expCoshInt = Mathematics::ExpCoshInt( tlow, thigh, gamma, deltaGamma, timeRes );
        double expSinhInt = Mathematics::ExpSinhInt( tlow, thigh, gamma, deltaGamma, timeRes );
        double expSinInt  = Mathematics::ExpSinInt(  tlow, thigh, gamma, deltaMs, timeRes );
        double expCosInt  = Mathematics::ExpCosInt(  tlow, thigh, gamma, deltaMs, timeRes );

        AzeroAzeroInt = Azero_sq * ( expCoshInt - cosPhis * expSinhInt + Btype * sinPhis * expSinInt );
        AparaAparaInt = Apara_sq * ( expCoshInt - cosPhis * expSinhInt + Btype * sinPhis * expSinInt );
        AperpAperpInt = Aperp_sq * ( expCoshInt + cosPhis * expSinhInt - Btype * sinPhis * expSinInt );
	AparaAperpInt = AparaAperp * (         -cosDeltaPerpMinusPara * sinPhis * expSinhInt
				      + Btype * sinDeltaPerpMinusPara * expCosInt
				      - Btype * cosDeltaPerpMinusPara * cosPhis * expSinInt);
	AzeroAparaInt = AzeroApara * cos(delta_para) * ( expCoshInt - cosPhis * expSinhInt
                                      			+ Btype * sinPhis * expSinInt);
	AzeroAperpInt = AzeroAperp * (         -cos(delta_perp) * sinPhis * expSinhInt
                                      + Btype * sin(delta_perp) * expCosInt
                                      - Btype * cos(delta_perp) * cosPhis * expSinInt);
	return;
}

/*
AzeroAzero[t_] := Exp[-gamma*t] * (   Cosh[deltaGamma*t/2]
                                    - cosPhis * Sinh[deltaGamma*t/2]
                                    + Btype * sinPhis * Sin[deltaMs*t] );
AzeroAzeroInt = CForm[FullSimplify[ Integrate[ AzeroAzero[t], {t, tmin, tmax}]]]
*/
/*
AparaApara[t_] := Exp[-gamma*t] * (   Cosh[deltaGamma*t/2]
                                    - cosPhis * Sinh[deltaGamma*t/2]
                                    + Btype * sinPhis * Sin[deltaMs*t] );
AparaAparaInt = CForm[FullSimplify[ Integrate[ AparaApara[t], {t, tmin, tmax}]]]
*/
/*
AperpAperp[t_] := Exp[-gamma*t] * (   Cosh[deltaGamma*t/2]
                                    + cosPhis * Sinh[deltaGamma*t/2]
                                    - Btype * sinPhis * Sin[deltaMs*t] );
AperpAperpInt = CForm[FullSimplify[ Integrate[ AperpAperp[t], {t, tmin, tmax}]]]
*/
/*
AparaAperp[t_] := Exp[-gamma*t] * (-cosDeltaPerpMinusPara * sinPhis * Sinh[deltaGamma*t/2]
                                    + Btype * sinDeltaPerpMinusPara * Cos[deltaMs*t]
                                    - Btype * cosDeltaPerpMinusPara * cosPhis * Sinh[deltaGamma*t/2] );
AparaAperpInt = CForm[FullSimplify[ Integrate[ AparaAperp[t], {t, tmin, tmax}]]]
*/
/*
AzeroApara[t_] := Exp[-gamma*t] * cosDeltaPara * (Cosh[deltaGamma*t/2]
                                    - cosPhis * sinh[deltaGamma*t/2]
                                    + Btype * sinPhis * Sin[deltaMs*t] );
AzeroAparaInt = CForm[FullSimplify[ Integrate[ AzeroApara[t], {t, tmin, tmax}]]]
*/
/*
AzeroAperp[t_] := Exp[-gamma*t] * ( - cosDeltaPerp * sinPhis * Sinh[deltaGamma*t/2]
                                    + Btype * sinDeltaPerp * Cos[deltaMs*t]
                                    - Btype * cosDeltaPerp * cosPhis * Sin[deltaMs*t] );
AzeroAperpInt = CForm[FullSimplify[ Integrate[ AzeroAperp[t], {t, tmin, tmax}]]]
*/
