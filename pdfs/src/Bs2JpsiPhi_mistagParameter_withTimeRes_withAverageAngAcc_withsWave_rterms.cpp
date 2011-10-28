// $Id: Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms.cpp,v 1.0 cofitzpa
/** @class Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi with average angular acceptance input as a paramter including s-wave component
 *
 *  @author Conor Fitzpatrick
 *  @date 2011-01-28
 */

#include "Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms );

//Constructor
Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms( PDFConfigurator* configurator ) : cachedAzeroAzeroIntB(),
cachedAparaAparaIntB(), cachedAperpAperpIntB(), cachedAparaAperpIntB(),
cachedAzeroAparaIntB(), cachedAzeroAperpIntB(), cachedAsAsIntB(), cachedAsAparaIntB(),
cachedAsAperpIntB(), cachedAsAzeroIntB(), cachedAzeroAzeroIntBbar(),
cachedAparaAparaIntBbar(), cachedAperpAperpIntBbar(), cachedAparaAperpIntBbar(),
cachedAzeroAparaIntBbar(), cachedAzeroAperpIntBbar(), cachedAsAsIntBbar(),
cachedAsAparaIntBbar(), cachedAsAperpIntBbar(), cachedAsAzeroIntBbar(),
cachedSinDeltaPerpPara(), cachedCosDeltaPerpPara(), cachedSinDeltaPerp(),
cachedCosDeltaPerp(), cachedCosDeltaPara(), cachedSinDeltaPerpS(), cachedSinDeltaParaS(),
cachedCosDeltaParaS(), cachedSinDeltaZeroS(), cachedCosDeltaZeroS(), cachedSinPhis(),
cachedCosPhis(), cachedExpCosh(), cachedExpSinh(), cachedExpCos(), cachedExpSin(),



	//Cache flags
	  normalisationCacheValid(false)
	, evaluationCacheValid(false)

	// Physics parameters
	, gammaName     ( configurator->getName("gamma") )
	, deltaGammaName( configurator->getName("deltaGamma") )
	, deltaMName    ( configurator->getName("deltaM") )
	, Phi_sName     ( configurator->getName("Phi_s") )
	, R_alphaName  ( configurator->getName("R_alpha") )
	, R_betaName  ( configurator->getName("R_beta") )
	, R_gammaName	( configurator->getName("R_gamma") )
	, delta_zeroName( configurator->getName("delta_zero") )
	, delta_paraName( configurator->getName("delta_para") )
	, delta_perpName( configurator->getName("delta_perp") )
	, delta_sName	( configurator->getName("delta_s") )
	// Angular acceptance factors
	, angAccI1Name	( configurator->getName("angAccI1") )
	, angAccI2Name	( configurator->getName("angAccI2") )
	, angAccI3Name	( configurator->getName("angAccI3") )
	, angAccI4Name	( configurator->getName("angAccI4") )
	, angAccI5Name	( configurator->getName("angAccI5") )
	, angAccI6Name	( configurator->getName("angAccI6") )
	, angAccI7Name	( configurator->getName("angAccI7") )
	, angAccI8Name	( configurator->getName("angAccI8") )
	, angAccI9Name	( configurator->getName("angAccI9") )
	, angAccI10Name	( configurator->getName("angAccI10") )
	// Detector parameters
	, mistagName	( configurator->getName("mistag") )
	, timeRes1Name	( configurator->getName("timeResolution1") )
	, timeRes2Name	( configurator->getName("timeResolution2") )
	, timeRes1FractionName	( configurator->getName("timeResolution1Fraction") )
	// Observables
	, timeName	( configurator->getName("time") )
	, cosThetaName	( configurator->getName("cosTheta") )
	, phiName	( configurator->getName("phi") )
	, cosPsiName	( configurator->getName("cosPsi") )
	, tagName	( configurator->getName("tag") )

	, timeconstraintName (configurator->getName("time") )
	// Member variables that will contain the parameter values
	, gamma(), deltaGamma(), deltaMs(), Phi_s(), R_alpha(),
	R_beta(), R_gamma(), Azero_sq(), Apara_sq(), Aperp_sq(),
	As_sq(), AzeroApara(), AzeroAperp(), AparaAperp(),
	AsApara(), AsAperp(), AsAzero(), delta_zero(),
	delta_para(), delta_perp(), delta_s(), omega(),
	timeRes(), timeRes1(), timeRes2(), timeRes1Frac(),
	angAccI1(), angAccI2(), angAccI3(), angAccI4(),
	angAccI5(), angAccI6(), angAccI7(), angAccI8(),
	angAccI9(), angAccI10(), time(), cosTheta(), phi(),
	cosPsi(), q(), tlow(), thigh()
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::MakePrototypes()
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
	parameterNames.push_back( R_alphaName );
	parameterNames.push_back( R_betaName );
	parameterNames.push_back( R_gammaName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_zeroName );
	parameterNames.push_back( delta_sName );
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
	parameterNames.push_back( angAccI7Name );
	parameterNames.push_back( angAccI8Name );
	parameterNames.push_back( angAccI9Name );
	parameterNames.push_back( angAccI10Name );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::~Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	deltaGamma = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	deltaMs    = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
	Phi_s      = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
	R_alpha   = allParameters.GetPhysicsParameter( R_alphaName )->GetValue();
	R_beta   = allParameters.GetPhysicsParameter( R_betaName )->GetValue();
	R_gamma   = allParameters.GetPhysicsParameter( R_gammaName )->GetValue();
	delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	delta_s = allParameters.GetPhysicsParameter( delta_sName )->GetValue();
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
	angAccI7 = allParameters.GetPhysicsParameter( angAccI7Name )->GetValue();
	angAccI8 = allParameters.GetPhysicsParameter( angAccI8Name )->GetValue();
	angAccI9 = allParameters.GetPhysicsParameter( angAccI9Name )->GetValue();
	angAccI10 = allParameters.GetPhysicsParameter( angAccI10Name )->GetValue();
	
	Azero_sq = R_alpha;						//R_alpha = 0 means Azero_sq = 0
	Aperp_sq = (1.0 - R_alpha)*R_beta;				//R_beta = 0 means Aperp_sq = 0
	As_sq = (1.0 - R_alpha)*(1.0 - R_beta)*R_gamma;			//R_gamma = 0 means As_sq = 0
	Apara_sq = (1.0 - R_alpha)*(1.0 - R_beta)*(1.0 - R_gamma);
	double sum = Apara_sq + Azero_sq + Aperp_sq + As_sq;
	//Apara_sq = 1 - Azero_sq - Aperp_sq - As_sq;

	cout << "Azero_sq: " << Azero_sq << " Aperp_sq: " << Aperp_sq << " As_sq: " << As_sq << " Apara_sq: " << Apara_sq << " sum: " << sum << endl; 
	if (R_alpha > 1.0 || R_alpha < 0.0 || R_beta > 1.0 || R_beta < 0.0 || R_gamma > 1.0 || R_gamma < 0.0){
	cerr << "Warning! Rterms are not on range [0,1]!: " << R_alpha << "," << R_beta << "," <<R_gamma << endl;
	return false;
	}

	AparaAperp = sqrt(Apara_sq)*sqrt(Aperp_sq);
	AzeroApara = sqrt(Azero_sq)*sqrt(Apara_sq);
	AzeroAperp = sqrt(Azero_sq)*sqrt(Aperp_sq);
	AsApara = sqrt(As_sq)*sqrt(Apara_sq);
	AsAperp = sqrt(As_sq)*sqrt(Aperp_sq);
	AsAzero = sqrt(As_sq)*sqrt(Azero_sq);

	return isOK;
}

//Return a list of parameters not to be integrated
vector<string> Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::GetDoNotIntegrateList()
{
	vector<string> doNotIntList;
	doNotIntList.push_back(mistagName);
	return doNotIntList;
}

//Calculate the function value
double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::Evaluate(DataPoint * measurement)
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

double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::buildPDFnumerator()
{
	// The angular functions f1->f6 as defined in roadmap Table 1.
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
	Mathematics::getBs2JpsiPhiAngularFunctionsWithSwave( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, cosTheta, phi, cosPsi );

	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59
	// First for the B
	double AzeroAzeroB, AparaAparaB, AperpAperpB,  AsAsB;
	double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB, ReAsAparaB, ImAsAperpB, ReAsAzeroB;
	getTimeDependentAmplitudes( 
			AzeroAzeroB, 
			AparaAparaB, 
			AperpAperpB, 
			AsAsB, 
			ImAparaAperpB, 
			ReAzeroAparaB, 
			ImAzeroAperpB, 
			ReAsAparaB, 
			ImAsAperpB, 
			ReAsAzeroB, 
			1
			);

	// Now for the Bbar
	double AzeroAzeroBbar, AparaAparaBbar, AperpAperpBbar, AsAsBbar;
	double ImAparaAperpBbar, ReAzeroAparaBbar, ImAzeroAperpBbar, ReAsAparaBbar, ImAsAperpBbar, ReAsAzeroBbar;
	getTimeDependentAmplitudes( 
			AzeroAzeroBbar, 
			AparaAparaBbar, 
			AperpAperpBbar, 
			AsAsBbar, 
			ImAparaAperpBbar, 
			ReAzeroAparaBbar, 
			ImAzeroAperpBbar, 
			ReAsAparaBbar, 
			ImAsAperpBbar, 
			ReAsAzeroBbar,
			-1
			);

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

	return ( w1*v1 + w2*v2 );
}


double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	// Now need to know the tag and the mistag
	q = (int)measurement->GetObservable( tagName )->GetValue();

	IConstraint * timeBound = boundary->GetConstraint(timeconstraintName);
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

double Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::buildPDFdenominator()
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
				, cachedAsAsIntB
				, cachedAparaAperpIntB
				, cachedAzeroAparaIntB
				, cachedAzeroAperpIntB
				, cachedAsAparaIntB
				, cachedAsAperpIntB
				, cachedAsAzeroIntB
				, 1);

		getTimeAmplitudeIntegrals(  cachedAzeroAzeroIntBbar
				, cachedAparaAparaIntBbar
				, cachedAperpAperpIntBbar
				, cachedAsAsIntBbar
				, cachedAparaAperpIntBbar
				, cachedAzeroAparaIntBbar
				, cachedAzeroAperpIntBbar
				, cachedAsAparaIntBbar
				, cachedAsAperpIntBbar
				, cachedAsAzeroIntBbar
				, -1);

		normalisationCacheValid = true;
	}

	double v1 = cachedAzeroAzeroIntB * angAccI1
		+ cachedAparaAparaIntB * angAccI2
		+ cachedAperpAperpIntB * angAccI3
		+ cachedAparaAperpIntB * angAccI4
		+ cachedAzeroAparaIntB * angAccI5
		+ cachedAzeroAperpIntB * angAccI6
		+ cachedAsAsIntB * angAccI7
		+ cachedAsAparaIntB * angAccI8
		+ cachedAsAperpIntB * angAccI9
		+ cachedAsAzeroIntB * angAccI10
		;

	double v2 = cachedAzeroAzeroIntBbar * angAccI1
		+ cachedAparaAparaIntBbar * angAccI2
		+ cachedAperpAperpIntBbar * angAccI3
		+ cachedAparaAperpIntBbar * angAccI4
		+ cachedAzeroAparaIntBbar * angAccI5
		+ cachedAzeroAperpIntBbar * angAccI6
		+ cachedAsAsIntBbar * angAccI7
		+ cachedAsAparaIntBbar * angAccI8
		+ cachedAsAperpIntBbar * angAccI9
		+ cachedAsAzeroIntBbar * angAccI10
		;

	//double norm = (angAccI1 + angAccI2 + angAccI3 + angAccI4 + angAccI5 + angAccI6)/(3*32*TMath::Pi()/9);
	double norm = 1.0;
	return (w1*v1 + w2*v2)/norm;
}

void Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::getTimeDependentAmplitudes(  
		double & AzeroAzero
		, double & AparaApara
		, double & AperpAperp
		, double & AsAs
		, double & ImAparaAperp
		, double & ReAzeroApara
		, double & ImAzeroAperp
		, double & ReAsApara
		, double & ImAsAperp
		, double & ReAsAzero
		, int Btype )
{
	// Quantities depending only on physics parameters can be cached
	if ( !evaluationCacheValid ){evaluateCache();}

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
	double cachedSinPhisExpSinh = cachedSinPhis * cachedExpSinh;
	double cachedCosPhisExpSinh = cachedCosPhis * cachedExpSinh;
	double cachedSinPhisExpSin = cachedSinPhis * cachedExpSin;
	double cachedCosPhisExpSin = cachedCosPhis * cachedExpSin;
	double Amp1 = cachedExpCosh - cachedCosPhisExpSinh + Btype * cachedSinPhisExpSin;
	double Amp2 = cachedExpCosh + cachedCosPhisExpSinh - Btype * cachedSinPhisExpSin;

	AzeroAzero   = Azero_sq   * Amp1;
	AparaApara   = Apara_sq   * Amp1;
	AperpAperp   = Aperp_sq   * Amp2;
	ImAparaAperp = AparaAperp * ( - cachedCosDeltaPerpPara * cachedSinPhisExpSinh + Btype * cachedSinDeltaPerpPara * cachedExpCos - Btype * cachedCosDeltaPerpPara * cachedCosPhisExpSin);
	ReAzeroApara = AzeroApara * cachedCosDeltaPara  * Amp1;
	ImAzeroAperp = AzeroAperp * ( - cachedCosDeltaPerp     * cachedSinPhisExpSinh + Btype * cachedSinDeltaPerp     * cachedExpCos - Btype * cachedCosDeltaPerp     * cachedCosPhisExpSin);
	AsAs         = As_sq      * Amp2;
	ReAsApara    = AsApara    * ( - cachedSinDeltaParaS    * cachedSinPhisExpSinh + Btype * cachedCosDeltaParaS    * cachedExpCos - Btype * cachedSinDeltaParaS    * cachedCosPhisExpSin);
	ImAsAperp    = AsAperp    * cachedSinDeltaPerpS * Amp2;
	ReAsAzero    = AsAzero    * ( - cachedSinDeltaZeroS    * cachedSinPhisExpSinh + Btype * cachedCosDeltaZeroS    * cachedExpCos - Btype * cachedSinDeltaZeroS    * cachedCosPhisExpSin);
	return;
}

void Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::getTimeAmplitudeIntegrals( 
		double & AzeroAzeroInt
		, double & AparaAparaInt
		, double & AperpAperpInt
		, double & AsAsInt
		, double & AparaAperpInt
		, double & AzeroAparaInt
		, double & AzeroAperpInt
		, double & AsAparaInt
		, double & AsAperpInt
		, double & AsAzeroInt
		, int Btype)
{
	if(!evaluationCacheValid){evaluateCache();}
	double expCoshInt = Mathematics::ExpCoshInt( tlow, thigh, gamma, deltaGamma, timeRes );
	double expSinhInt = Mathematics::ExpSinhInt( tlow, thigh, gamma, deltaGamma, timeRes );
	double expSinInt  = Mathematics::ExpSinInt(  tlow, thigh, gamma, deltaMs, timeRes );
	double expCosInt  = Mathematics::ExpCosInt(  tlow, thigh, gamma, deltaMs, timeRes );

	double cachedSinPhisExpSinhInt = cachedSinPhis * expSinhInt;
        double cachedCosPhisExpSinhInt = cachedCosPhis * expSinhInt;
        double cachedSinPhisExpSinInt = cachedSinPhis * expSinInt;
        double cachedCosPhisExpSinInt = cachedCosPhis * expSinInt;
	double Int1 = expCoshInt - cachedCosPhisExpSinhInt + Btype * cachedSinPhisExpSinInt;
	double Int2 = expCoshInt + cachedCosPhisExpSinhInt - Btype * cachedSinPhisExpSinInt;

	AzeroAzeroInt = Azero_sq * Int1;
	AparaAparaInt = Apara_sq * Int1;
	AperpAperpInt = Aperp_sq * Int2;
	AparaAperpInt = AparaAperp * ( -cachedCosDeltaPerpPara * cachedSinPhisExpSinhInt + Btype * cachedSinDeltaPerpPara * expCosInt - Btype * cachedCosDeltaPerpPara * cachedCosPhisExpSinInt);
	AzeroAparaInt = AzeroApara * cachedCosDeltaPara  * Int1;
	AzeroAperpInt = AzeroAperp * ( -cachedCosDeltaPerp     * cachedSinPhisExpSinhInt + Btype * cachedSinDeltaPerp     * expCosInt - Btype * cachedCosDeltaPerp     * cachedCosPhisExpSinInt);
	AsAsInt       = As_sq    * Int2;	
	AsAparaInt = AsApara       * ( -cachedSinDeltaParaS    * cachedSinPhisExpSinhInt + Btype * cachedCosDeltaParaS    * expCosInt - Btype * cachedSinDeltaParaS    * cachedCosPhisExpSinInt);
	AsAperpInt = AsAperp       * cachedSinDeltaPerpS * Int2;
	AsAzeroInt = AsAzero       * ( -cachedSinDeltaZeroS    * cachedSinPhisExpSinhInt + Btype * cachedCosDeltaZeroS    * expCosInt - Btype * cachedSinDeltaZeroS    * cachedCosPhisExpSinInt);


	return;
}

void Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms::evaluateCache(){
	cachedSinDeltaPerpPara  = sin( delta_perp - delta_para );
	cachedCosDeltaPerpPara  = cos( delta_perp - delta_para );
	cachedSinDeltaPerp      = sin( delta_perp );
	cachedCosDeltaPerp      = cos( delta_perp );
	cachedCosDeltaPara      = cos( delta_para );
	cachedSinDeltaPerpS     = sin(delta_perp - delta_s);
	cachedSinDeltaParaS     = sin(delta_para - delta_s);
	cachedCosDeltaParaS     = cos(delta_para - delta_s);
	cachedSinDeltaZeroS     = sin(-delta_s);
	cachedCosDeltaZeroS     = cos(-delta_s);
	cachedSinPhis = sin( Phi_s );
	cachedCosPhis = cos( Phi_s );
	evaluationCacheValid = true;
}

