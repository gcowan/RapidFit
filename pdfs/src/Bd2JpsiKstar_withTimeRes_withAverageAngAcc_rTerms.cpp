// $Id: Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.cpp
 *
 *  RatagFit PDF for Bd2JpsiPhi with average angular acceptance input as a paramter
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2010-01-12
 */

#include "TMath.h"
#include <cmath>

#include "Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"

PDF_CREATOR( Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms );

//Constructor
Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms( PDFConfigurator* configurator ) :
	cachedAzeroAzeroIntB(), cachedAparaAparaIntB(),
	cachedAperpAperpIntB(), cachedAparaAperpIntB(),
	cachedAzeroAparaIntB(), cachedAzeroAperpIntB(),
	cachedSinDeltaPerpPara(), cachedCosDeltaPara(),
	cachedSinDeltaPerp(), cachedCosSqDeltaM(),
	cachedAzero(), cachedApara(), cachedAperp(),


	  normalisationCacheValid(false)
	, evaluationCacheValid(false)

	// Physics parameters
	,  gammaName     ( configurator->getName("gamma") )
	, deltaMName    ( configurator->getName("deltaM") )
	, R_alphaName  ( configurator->getName("R_alpha") )
        , R_betaName  ( configurator->getName("R_beta") )
	, delta_zeroName( configurator->getName("delta_zero") )
	, delta_paraName( configurator->getName("delta_para") )
	, delta_perpName( configurator->getName("delta_perp") )
//	, Aperp_sqName  ( configurator->getName("Aperp_sq") )
//	, Apara_sqName  ( configurator->getName("Apara_sq") )
	, angAccI1Name	( configurator->getName("angAccI1") )
	, angAccI2Name	( configurator->getName("angAccI2") )
	, angAccI3Name	( configurator->getName("angAccI3") )
	, angAccI4Name	( configurator->getName("angAccI4") )
	, angAccI5Name	( configurator->getName("angAccI5") )
	, angAccI6Name	( configurator->getName("angAccI6") )
	, timeRes1Name	( configurator->getName("timeResolution1") )
	, timeRes2Name	( configurator->getName("timeResolution2") )
	, timeRes1FractionName	( configurator->getName("timeResolution1Fraction") )

	// Observables (What we want to gain from the pdf after inserting physics parameter values)
	, timeName	( configurator->getName("time") )
	, cosThetaName	( configurator->getName("cosTheta") )
	, phiName	( configurator->getName("phi") )
	, cosPsiName	( configurator->getName("cosPsi") )
	, KstarFlavourName  ( configurator->getName("KstarFlavour") )

	, timeconstraintName( configurator->getName("time") )
	, gamma(), deltaMs(), R_alpha(), R_beta(), Azero_sq(),
	Apara_sq(), Aperp_sq(), AzeroApara(), AzeroAperp(),
	AparaAperp(), delta_zero(), delta_para(), delta_perp(),
	omega(), timeRes(), timeRes1(), timeRes2(),
	timeRes1Frac(), angAccI1(), angAccI2(), angAccI3(),
	angAccI4(), angAccI5(), angAccI6(), time(), cosTheta(),
	phi(), cosPsi(), KstarFlavour(), tlow(), thigh()
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::MakePrototypes()
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
	 parameterNames.push_back( R_alphaName );
        parameterNames.push_back( R_betaName );

	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_zeroName );
	parameterNames.push_back( deltaMName );
	parameterNames.push_back( timeRes1Name );
	parameterNames.push_back( timeRes2Name );
	parameterNames.push_back( timeRes1FractionName );
	parameterNames.push_back( angAccI1Name );
	parameterNames.push_back( angAccI2Name );
	parameterNames.push_back( angAccI3Name );
	parameterNames.push_back( angAccI4Name );
	parameterNames.push_back( angAccI5Name );
	parameterNames.push_back( angAccI6Name );
	allParameters = ParameterSet(parameterNames);
}

//Destructor
Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::~Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
        // Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
        gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
        deltaMs    = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
        //Aperp_sq   = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
        //Apara_sq   = allParameters.GetPhysicsParameter( Apara_sqName )->GetValue();
	R_alpha   = allParameters.GetPhysicsParameter( R_alphaName )->GetValue();
        R_beta   = allParameters.GetPhysicsParameter( R_betaName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
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

        /*Azero_sq = 1 - Aperp_sq - Apara_sq;
		 if ( Azero_sq < 0.) return false ;
		AparaAperp = sqrt(Apara_sq)*sqrt(Aperp_sq);
        	AzeroApara = sqrt(Azero_sq)*sqrt(Apara_sq);
        	AzeroAperp = sqrt(Azero_sq)*sqrt(Aperp_sq);

	Apara_sq = allParameters.GetPhysicsParameter( Apara_sqName )->GetValue();
        if( (Apara_sq < 0.) || (Apara_sq > 1.)  ) { cout << "Warning in Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::SetPhysicsParameters: Apara_sq <0 or >1 but left as is" <<  endl ;      }
        Aperp_sq = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
        if( (Aperp_sq < 0.) || (Aperp_sq > 1.)  ) { cout << "Warning in Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::SetPhysicsParameters: Aperp_sq <0 or >1 but left as is" <<  endl ;      }

        Azero_sq = (1. - Apara_sq - Aperp_sq ) ;
        if( Azero_sq < 0. ) {
                cout << "Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::SetPhysicsParameters: derived parameter Apara_sq <0  and so set to zero" <<  endl ;
                Azero_sq = 0. ;
        }
*/

	Aperp_sq = R_alpha;                                             //R_alpha = 0 means Azero_sq = 0
        Apara_sq = (1.0 - R_alpha)*R_beta;                              //R_beta = 0 means Aperp_sq = 0
        Azero_sq = (1.0 - R_alpha)*(1.0 - R_beta);
//        double sum = Apara_sq + Azero_sq + Aperp_sq;
        //Apara_sq = 1 - Azero_sq - Aperp_sq - As_sq;

//        cout << "Azero_sq: " << Azero_sq << " Aperp_sq: " << Aperp_sq << " Apara_sq: " << Apara_sq << " sum: " << sum << endl;

	if (R_alpha > 1.0 || R_alpha < 0.0 || R_beta > 1.0 || R_beta < 0.0 ){
        cerr << "Warning! Rterms are not on range [0,1]!: " << R_alpha << "," << R_beta << endl;
        return false;
        }


	return isOK;
}

double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::q() const { return KstarFlavour;}


//Return a list of parameters not to be integrated
vector<string> Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::GetDoNotIntegrateList()
{
	vector<string> doNotIntList;
	return doNotIntList;
}

//Calculate the function value
double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::Evaluate(DataPoint * measurement)
{
	time = measurement->GetObservable( timeName )->GetValue();
        cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
        phi      = measurement->GetObservable( phiName )->GetValue();
        cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();
	KstarFlavour = measurement->GetObservable( KstarFlavourName )->GetValue();
	//cout << gamma << " " << Aperp_sq << " " << Azero_sq << endl;

	double returnValue;

        if(timeRes1Frac >= 0.9999)
	{
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                returnValue =  buildPDFnumerator();
	//	return returnValue;
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
	//	return returnValue; //added by Ailsa
        }


	        if(  ( (returnValue <= 0.) && (time>0.) ) || std::isnan(returnValue) ) {
                cout << endl ;
                cout << " Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::evaluate() returns <=0 or nan :" << returnValue << endl ;
                cout << "   Aperp    " << Aperp_sq ;
                cout << "   APara    " << Apara_sq ;
                cout << "   A0    " << Azero_sq << endl ;
		cout << "   Dperp    " << delta_perp;
		cout << "   Dpara    " << delta_para;
		cout << "   Dzero    " << delta_zero << endl;
		cout << "   gamma   " << gamma << endl;
        cout << " For event with: " << endl ;
                cout << "   time " << time << endl ;
                if( std::isnan(returnValue) ) exit(1) ;
        }

	return returnValue;

}

//Amplitudes Used in three angle PDF
double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::AT() const {
        if( Aperp_sq <= 0. ) return 0. ;
        else return sqrt(Aperp_sq) ;
}
double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::AP() const {
        if( Apara_sq <= 0. ) return 0. ;
        else return sqrt(Apara_sq) ;
}
double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::A0() const {
        if( Azero_sq <= 0. ) return 0. ;
        else return sqrt(Azero_sq) ;
}

double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::DT() const {
        if( delta_perp <= 0. ) return 0. ;
        else return delta_perp ;
}
double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::DP() const {
        if( delta_para <= 0. ) return 0. ;
        else return delta_para ;
}
double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::D0() const {
        if( delta_zero <= 0. ) return 0. ;
        else return delta_zero ;
}


double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::Gamma() const   {
	return gamma ;
}


double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::buildPDFnumerator()
{
	// The angular functions f1->f6 as defined in roadmap Table 1.(same for Kstar)
	double f1, f2, f3, f4, f5, f6;
	Mathematics::getBs2JpsiPhiAngularFunctions( f1, f2, f3, f4, f5, f6, cosTheta, phi, cosPsi );

	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59  //No tagging so only need 2 (h± pg 72)
	// First for the B
	double AzeroAzeroB, AparaAparaB, AperpAperpB;
	double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB;

	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB
			, ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
							   );

//	double q = KstarFlavour;
	//W+
	double v1 = f1 * AzeroAzeroB
		+ f2 * AparaAparaB
		+ f3 * AperpAperpB
		+ q() * f4 * ImAparaAperpB
		+ f5 * ReAzeroAparaB
		+ q() * f6 * ImAzeroAperpB;

	/*
	if (std::isnan(v1))
	{
		cout << f1 << " " << f2 << " " << f3 << " " <<f4 << " " << f5 << " " << f6 << endl;
			cout << AzeroAzeroB << " " << AperpAperpB << " " << AparaAparaB << " " << ImAparaAperpB << " " << ReAzeroAparaB << " " << ImAzeroAperpB << endl;
	}
	 */
		return v1;
}


double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	double returnValue;
	KstarFlavour = measurement->GetObservable( KstarFlavourName )->GetValue();
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
                cout << " Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::Normalisation() returns <=0 or nan " << endl ;
                cout << " AT    " << Aperp_sq ;
                cout << " AP    " << Apara_sq ;
                cout << " A0    " << Azero_sq ;
		cout << "   Dperp    " << delta_perp;
                cout << "   Dpara    " << delta_para;
                cout << "   Dzero   " << delta_zero << endl;
		cout << "   gamma   " << gamma << endl;

		exit(1) ;
        }

        return returnValue ;

}

double Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::buildPDFdenominator()
{

	if (!normalisationCacheValid)
	{
		// The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
		getTimeAmplitudeIntegrals(  cachedAzeroAzeroIntB
				, cachedAparaAparaIntB
				, cachedAperpAperpIntB
				, cachedAparaAperpIntB
				, cachedAzeroAparaIntB
								  , cachedAzeroAperpIntB);


		normalisationCacheValid = true;
	}

	 //double q = KstarFlavour;

	double v1 = cachedAzeroAzeroIntB * angAccI1
		+ cachedAparaAparaIntB * angAccI2
		+ cachedAperpAperpIntB * angAccI3
		+ q() * cachedAparaAperpIntB * angAccI4
		+ cachedAzeroAparaIntB * angAccI5
		+ q() * cachedAzeroAperpIntB * angAccI6;

	return v1;
}

void Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::getTimeDependentAmplitudes(  double & AzeroAzero
		, double & AparaApara
		, double & AperpAperp
		, double & ImAparaAperp
		, double & ReAzeroApara
		, double & ImAzeroAperp		)
{
	// Quantities depending only on physics parameters can be cached
	if ( !evaluationCacheValid )
	{
		cachedAzero = sqrt( Azero_sq );
		cachedApara = sqrt( Apara_sq );
		cachedAperp = sqrt( Aperp_sq );

		cachedSinDeltaPerpPara	= sin( delta_perp - delta_para );
		cachedCosDeltaPara	= cos( delta_para );
		cachedSinDeltaPerp	= sin( delta_perp );

	//	cachedCosSqDeltaM	= cos((deltaMs*timeRes)/2)*cos((deltaMs*timeRes)/2);


		evaluationCacheValid = true;
	}


	//cout << gamma << " " << deltaGamma << " " << Azero_sq << " " << Aperp_sq << endl;
	//cout << cachedExpCosh << " " << cachedExpSinh << " " << cachedExpCos << " " << cachedExpSin << endl;
	// Now calculate the amplitudes

	double Exp = Mathematics::Exp(time, gamma, timeRes);


	AzeroAzero = Azero_sq * Exp;  // changed- see note 2009-015 eq 11-13
        AparaApara = Apara_sq * Exp;  //
        AperpAperp = Aperp_sq * Exp;  //

        ImAparaAperp = cachedApara*cachedAperp * cachedSinDeltaPerpPara * Exp;    //See http://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=33933 page14

	ReAzeroApara = cachedAzero*cachedApara * cachedCosDeltaPara * Exp;

        ImAzeroAperp = cachedAzero*cachedAperp * cachedSinDeltaPerp * Exp;

	//if ( std::isnan(ImAparaAperp)) cout << Azero_sq << " " << Apara_sq << " " << Aperp_sq << " " << Exp << endl;

	return;
}

void Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms::getTimeAmplitudeIntegrals( double & AzeroAzeroInt
		, double & AparaAparaInt
		, double & AperpAperpInt
		, double & AparaAperpInt
		, double & AzeroAparaInt
		, double & AzeroAperpInt	)
{

	double ExpInt = Mathematics::ExpInt(tlow, thigh, gamma, timeRes);

		AzeroAzeroInt = Azero_sq * ExpInt;
	        AparaAparaInt = Apara_sq  * ExpInt;
       		AperpAperpInt = Aperp_sq * ExpInt;

		AparaAperpInt = AparaAperp * cachedSinDeltaPerpPara * ExpInt;

		AzeroAparaInt = AzeroApara * cachedCosDeltaPara * ExpInt;

		AzeroAperpInt = AzeroAperp * cachedSinDeltaPerp * ExpInt;


	return;
}
