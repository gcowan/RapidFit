/** @class Bd2JpsiKstar_sWave_rTerms Bd2JpsiKstar_sWave.cpp
 *
 *  RapidFit PDF for Bd2JpsiPhi with average angular acceptance input as a paramter
 *
 *  @author Ailsa Sparkes asparkes@cern.ch
 *  @date 2011-01-26
 */

#include "TMath.h"
#include <cmath>

#include "Bd2JpsiKstar_sWave_rTerms.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"

PDF_CREATOR( Bd2JpsiKstar_sWave_rTerms );

//Constructor
Bd2JpsiKstar_sWave_rTerms::Bd2JpsiKstar_sWave_rTerms( PDFConfigurator* configurator ) : cachedAzeroAzeroIntB(), cachedAparaAparaIntB(), cachedAperpAperpIntB(), cachedAparaAperpIntB(),
	cachedAzeroAparaIntB(), cachedAzeroAperpIntB(), cachedAsAsIntB(), cachedAparaAsIntB(), cachedAperpAsIntB(), cachedAzeroAsIntB(), cachedSinDeltaPerpPara(),
	cachedCosDeltaPara(), cachedSinDeltaPerp(), cachedCosDeltaParaS(), cachedSinDeltaPerpS(), cachedCosDeltaS(), cachedAzero(), cachedApara(), cachedAperp(),
	cachedAs(),
         normalisationCacheValid(false)
	, evaluationCacheValid(false)
	// Physics parameters
	, gammaName     ( configurator->getName("gamma") )
	, deltaMName    ( configurator->getName("deltaM") )
	, R_alphaName   ( configurator->getName("R_alpha")  )
        , R_betaName    ( configurator->getName("R_beta")   )
        , R_gammaName   ( configurator->getName("R_gamma")  )
	, delta_zeroName( configurator->getName("delta_zero") )
 	, delta_paraName( configurator->getName("delta_para") )
        , delta_perpName( configurator->getName("delta_perp") )
	//, Aperp_sqName  ( configurator->getName("Aperp_sq") )
	//, Apara_sqName  ( configurator->getName("Apara_sq") )
	//, As_sqName	( configurator->getName("As_sq") )
        , delta_sName( configurator->getName("delta_s") )
	, angAccI1Name	( configurator->getName("angAccI1") )
	, angAccI2Name	( configurator->getName("angAccI2") )
	, angAccI3Name	( configurator->getName("angAccI3") )
	, angAccI4Name	( configurator->getName("angAccI4") )
	, angAccI5Name	( configurator->getName("angAccI5") )
	, angAccI6Name	( configurator->getName("angAccI6") )
	, angAccI7Name  ( configurator->getName("angAccI7") )
	, angAccI8Name  ( configurator->getName("angAccI8") )
	, angAccI9Name  ( configurator->getName("angAccI9") )
	, angAccI10Name  ( configurator->getName("angAccI10") )
	, timeRes1Name	( configurator->getName("timeResolution1") )
	, timeRes2Name	( configurator->getName("timeResolution2") )
	, timeRes1FractionName	( configurator->getName("timeResolution1Fraction") )

	, gamma(), deltaMs(), Azero_sq(), Apara_sq(), Aperp_sq(), As_sq(), R_alpha(), R_beta(), R_gamma(), AzeroApara(), AzeroAperp(), AparaAperp(), AparaAs(),
	AperpAs(), AzeroAs(), delta_zero(), delta_para(), delta_perp(), delta_s(), omega(), timeRes(), timeRes1(), timeRes2(), timeRes1Frac(), angAccI1(),
	angAccI2(), angAccI3(), angAccI4(), angAccI5(), angAccI6(), angAccI7(), angAccI8(), angAccI9(), angAccI10()

	// Observables (What we want to gain from the pdf after inserting physics parameter values)
	, timeName	( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	( "phi" )
	, cosPsiName	( "cosPsi" )
	, KstarFlavourName  ( "KstarFlavour" )

	, timeconstraintName( "time" )
	,time(), cosTheta(), phi(), cosPsi(), KstarFlavour(), tlow(), thigh()
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bd2JpsiKstar_sWave_rTerms::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );
	allObservables.push_back( KstarFlavourName );


	// Need to think about additional parameters like
	// event-by-event propertime resolution and acceptance.
	// This will require event-by-event PDF normalisation,
	// but we are already doing this for tagging.

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	//parameterNames.push_back( Apara_sqName );
	//parameterNames.push_back( Aperp_sqName );
	//parameterNames.push_back( As_sqName );
	parameterNames.push_back( R_alphaName );
        parameterNames.push_back( R_betaName );
        parameterNames.push_back( R_gammaName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_zeroName );
	parameterNames.push_back( deltaMName );
	parameterNames.push_back( delta_sName );
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
allParameters = ParameterSet(parameterNames);
}

//Destructor
Bd2JpsiKstar_sWave_rTerms::~Bd2JpsiKstar_sWave_rTerms()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bd2JpsiKstar_sWave_rTerms::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
        // Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
        gamma = allParameters.GetPhysicsParameter( gammaName )->GetValue();
//	gamma = 1/temp;
	deltaMs    = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
        //Aperp_sq   = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
        //Apara_sq   = allParameters.GetPhysicsParameter( Apara_sqName )->GetValue();
	//As_sq   = allParameters.GetPhysicsParameter( As_sqName )->GetValue();
	R_alpha   = allParameters.GetPhysicsParameter( R_alphaName )->GetValue();
        R_beta   = allParameters.GetPhysicsParameter( R_betaName )->GetValue();
        R_gamma   = allParameters.GetPhysicsParameter( R_gammaName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	delta_s = allParameters.GetPhysicsParameter( delta_sName )->GetValue();
	delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
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

//Azero_sq = 1 - Aperp_sq - Apara_sq;
  	Aperp_sq = R_alpha;                                             //R_alpha = 0 means Azero_sq = 0
        Apara_sq = (1.0 - R_alpha)*R_beta;                              //R_beta = 0 means Aperp_sq = 0
        Azero_sq = (1.0 - R_alpha)*(1.0 - R_beta)*R_gamma;                 //R_gamma = 0 means As_sq = 0
        As_sq = (1.0 - R_alpha)*(1.0 - R_beta)*(1.0 - R_gamma);
//        double sum = Apara_sq + Azero_sq + Aperp_sq + As_sq;
        //Apara_sq = 1 - Azero_sq - Aperp_sq - As_sq;
        //cout << "Azero_sq: " << Azero_sq << " Aperp_sq: " << Aperp_sq << " As_sq: " << As_sq << " Apara_sq: " << Apara_sq << " sum: " << sum << endl;
        if (R_alpha > 1.0 || R_alpha < 0.0 || R_beta > 1.0 || R_beta < 0.0 || R_gamma > 1.0 || R_gamma < 0.0){
        cerr << "Warning! Rterms are not on range [0,1]!: " << R_alpha << "," << R_beta << "," <<R_gamma << endl;
        return false;
        }


		AparaAperp = sqrt(Apara_sq)*sqrt(Aperp_sq);
        	AzeroApara = sqrt(Azero_sq)*sqrt(Apara_sq);
        	AzeroAperp = sqrt(Azero_sq)*sqrt(Aperp_sq);
		AparaAs    = sqrt(Apara_sq)*sqrt(As_sq);
		AperpAs	   = sqrt(Aperp_sq)*sqrt(As_sq);
		AzeroAs    = sqrt(Azero_sq)*sqrt(As_sq);

	return isOK;
}

//Return a list of parameters not to be integrated
vector<string> Bd2JpsiKstar_sWave_rTerms::GetDoNotIntegrateList()
{
	vector<string> doNotIntList;
	return doNotIntList;
}

double Bd2JpsiKstar_sWave_rTerms::q() const { return KstarFlavour;}


//Calculate the function value
double Bd2JpsiKstar_sWave_rTerms::Evaluate(DataPoint * measurement)
{
	time = measurement->GetObservable( timeName )->GetValue();
        cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
        phi      = measurement->GetObservable( phiName )->GetValue();
        cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();
	KstarFlavour = measurement->GetObservable( KstarFlavourName )->GetValue();

 double returnValue;

        if(timeRes1Frac >= 0.9999)
        {
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                returnValue =  buildPDFnumerator();
        //      return returnValue;
}
        else
        {
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                double val1 = buildPDFnumerator();
                // Set the member variable for time resolution to the second value and calculate
                timeRes = timeRes2;
                double val2 = buildPDFnumerator();

                returnValue = timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
        //      return returnValue; //added by Ailsa
        }


                if(  ( (returnValue <= 0.) && (time>0.) ) || std::isnan(returnValue) ) {
                cout << endl ;
                cout << " Bd2JpsiKstar_sWave_rTerms::evaluate() returns <=0 or nan :" << returnValue << endl ;
                cout << "   Aperp    " << Aperp_sq ;
                cout << "   APara    " << Apara_sq ;
                cout << "   A0    " << Azero_sq;
		cout << "   As   " << As_sq << endl;
                cout << "   Dperp    " << delta_perp;
                cout << "   Dpara    " << delta_para;
                cout << "   Dzero    " << delta_zero;
		cout << "   Ds     " << delta_s << endl;
                cout << "   gamma   " << gamma << endl;
        cout << " For event with: " << endl ;
                cout << "   time " << time << endl ;
                if( std::isnan(returnValue) ) exit(1) ;
        }

        return returnValue;

}

//Amplitudes Used in three angle PDF
double Bd2JpsiKstar_sWave_rTerms::AT() const {
        if( Aperp_sq <= 0. ) return 0. ;
        else return sqrt(Aperp_sq) ;
}
double Bd2JpsiKstar_sWave_rTerms::AP() const {
        if( Apara_sq <= 0. ) return 0. ;
        else return sqrt(Apara_sq) ;
}
double Bd2JpsiKstar_sWave_rTerms::A0() const {
        if( Azero_sq <= 0. ) return 0. ;
        else return sqrt(Azero_sq) ;
}



double Bd2JpsiKstar_sWave_rTerms::Gamma() const   {
        return gamma ;
}


double Bd2JpsiKstar_sWave_rTerms::buildPDFnumerator()
{
	// The angular functions f1->f6 as defined in roadmap Table 1.(same for Kstar)
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
	Mathematics::getBs2JpsiPhiAngularFunctionsWithSwave( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, cosTheta, phi, cosPsi );

	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59  //No tagging so only need 2 (h± pg 72)
	// First for the B
	double AzeroAzeroB, AparaAparaB, AperpAperpB, AsAsB;
	double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB;
	double ReAparaAsB, ImAperpAsB, ReAzeroAs;

	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB
			, ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
			, AsAsB, ReAparaAsB, ImAperpAsB, ReAzeroAs
							   );


	double v1 = f1 * AzeroAzeroB
		+ f2 * AparaAparaB
		+ f3 * AperpAperpB
		+ q() * f4 * ImAparaAperpB
		+ f5 * ReAzeroAparaB
		+ q() * f6 * ImAzeroAperpB
		+ f7 * AsAsB
		+ f8 * ReAparaAsB
		+ f9 * ImAperpAsB
		+ f10 * ReAzeroAs;

	/*
	if (std::isnan(v1))
	{
		cout << f1 << " " << f2 << " " << f3 << " " <<f4 << " " << f5 << " " << f6 << endl;
			cout << AzeroAzeroB << " " << AperpAperpB << " " << AparaAparaB << " " << ImAparaAperpB << " " << ReAzeroAparaB << " " << ImAzeroAperpB << endl;
	}
	 */
		return v1;
}


double Bd2JpsiKstar_sWave_rTerms::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{

        double returnValue;
        KstarFlavour = measurement->GetObservable( KstarFlavourName )->GetValue();
        IConstraint * timeBound = boundary->GetConstraint( timeconstraintName );
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
                returnValue =  buildPDFdenominator();
        }
        else
        {
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                double val1 = buildPDFdenominator();
                // Set the member variable for time resolution to the second value and calculate
                timeRes = timeRes2;
                double val2 = buildPDFdenominator();
                returnValue =  timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
        }

	      if( (returnValue <= 0.) || std::isnan(returnValue) ) {
                cout << " Bd2JpsiKstar_sWave_rTerms::Normalisation() returns <=0 or nan " << endl ;
                cout << " AT    " << Aperp_sq ;
                cout << " AP    " << Apara_sq ;
                cout << " A0    " << Azero_sq;
		cout << " As   " << As_sq;
                cout << "   Dperp    " << delta_perp;
                cout << "   Dpara    " << delta_para;
                cout << "   Dzero   " << delta_zero;
                cout << "   Ds     " << delta_s << endl;
		cout << "   gamma   " << gamma << endl;

                exit(1) ;
        }

        return returnValue ;

}

double Bd2JpsiKstar_sWave_rTerms::buildPDFdenominator()
{

	if (!normalisationCacheValid)
	{
		// The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
		getTimeAmplitudeIntegrals(  cachedAzeroAzeroIntB
				, cachedAparaAparaIntB
				, cachedAperpAperpIntB
				, cachedAparaAperpIntB
				, cachedAzeroAparaIntB
				, cachedAzeroAperpIntB
				, cachedAsAsIntB
				, cachedAparaAsIntB
				, cachedAperpAsIntB
				, cachedAzeroAsIntB
						);


		normalisationCacheValid = true;
	}

	//cout << gamma << endl;

	double v1 = cachedAzeroAzeroIntB * angAccI1
		+ cachedAparaAparaIntB * angAccI2
		+ cachedAperpAperpIntB * angAccI3
		+ q() * cachedAparaAperpIntB * angAccI4
		+ cachedAzeroAparaIntB * angAccI5
		+ q() * cachedAzeroAperpIntB * angAccI6
		+ cachedAsAsIntB * angAccI7
		+ cachedAparaAsIntB * angAccI8
		+ cachedAperpAsIntB * angAccI9
		+ cachedAzeroAsIntB * angAccI10;

	return v1;
}

void Bd2JpsiKstar_sWave_rTerms::getTimeDependentAmplitudes(
		  double & AzeroAzero
		, double & AparaApara
		, double & AperpAperp
		, double & ImAparaAperp
		, double & ReAzeroApara
		, double & ImAzeroAperp
		, double & AsAs
		, double & ReAparaAs
		, double & ImAperpAs
		, double & ReAzeroAs
		)
{
	(void)AsAs;

	// Quantities depending only on physics parameters can be cached
	if ( !evaluationCacheValid )
	{
		cachedAzero = sqrt( Azero_sq );
		cachedApara = sqrt( Apara_sq );
		cachedAperp = sqrt( Aperp_sq );
		cachedAs = sqrt (As_sq);

		cachedSinDeltaPerpPara	= sin( delta_perp - delta_para );
		cachedCosDeltaPara	= cos( delta_para );
		cachedSinDeltaPerp	= sin( delta_perp );
		cachedCosDeltaParaS	= cos( delta_para - delta_s );
		cachedSinDeltaPerpS 	= sin( delta_perp - delta_s );
		cachedCosDeltaS		= cos( delta_s);

		evaluationCacheValid = true;
	}


	//cout << gamma << " " << deltaGamma << " " << Azero_sq << " " << Aperp_sq << endl;

	//cout << cachedExpCosh << " " << cachedExpSinh << " " << cachedExpCos << " " << cachedExpSin << endl;


	// Now calculate the amplitudes

	double Exp = Mathematics::Exp(time, gamma, timeRes);

	AzeroAzero = Azero_sq * Exp;  // changed- see note 2009-015 eq 11-13
        AparaApara = Apara_sq * Exp;  //
        AperpAperp = Aperp_sq * Exp;  //
	AzeroAzero = Azero_sq * Exp;

        ImAparaAperp = cachedApara*cachedAperp * cachedSinDeltaPerpPara * Exp;    //See http://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=33933 page14

	ReAzeroApara = cachedAzero*cachedApara * cachedCosDeltaPara * Exp;

        ImAzeroAperp = cachedAzero*cachedAperp * cachedSinDeltaPerp * Exp;

	ReAparaAs = cachedApara*cachedAs * cachedCosDeltaParaS * Exp; //AILSA_ NOT SURE

	ImAperpAs = cachedAperp*cachedAs * cachedSinDeltaPerpS * Exp; //AILSA

	ReAzeroAs = cachedAzero*cachedAs * cachedCosDeltaS * Exp; //AILSA

	//if ( std::isnan(ImAparaAperp)) cout << Azero_sq << " " << Apara_sq << " " << Aperp_sq << " " << Exp << endl;

	return;
}

void Bd2JpsiKstar_sWave_rTerms::getTimeAmplitudeIntegrals(
		  double & AzeroAzeroInt
		, double & AparaAparaInt
		, double & AperpAperpInt
		, double & AparaAperpInt
		, double & AzeroAparaInt
		, double & AzeroAperpInt
		, double & AsAsInt
		, double & AparaAsInt
		, double & AperpAsInt
		, double & AzeroAsInt
	)
{

	double ExpInt = Mathematics::ExpInt(tlow, thigh, gamma, timeRes);

		AzeroAzeroInt = Azero_sq * ExpInt;
	        AparaAparaInt = Apara_sq  * ExpInt;
       		AperpAperpInt = Aperp_sq * ExpInt;
		AsAsInt = As_sq * ExpInt;

		AparaAperpInt = AparaAperp * cachedSinDeltaPerpPara * ExpInt;

		AzeroAparaInt = AzeroApara * cachedCosDeltaPara * ExpInt;

		AzeroAperpInt = AzeroAperp * cachedSinDeltaPerp * ExpInt;

		AparaAsInt = AparaAs * cachedCosDeltaParaS * ExpInt;

		AperpAsInt = AperpAs * cachedSinDeltaPerpS * ExpInt;

		AzeroAsInt = AzeroAs * cachedCosDeltaS * ExpInt;

	return;
}
