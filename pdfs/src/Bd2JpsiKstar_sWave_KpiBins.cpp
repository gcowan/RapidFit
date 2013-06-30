/** @class Bd2JpsiKstar_sWave_KpiBins Bd2JpsiKstar_sWave_KpiBins.cpp
 *
 *  RapidFit PDF for Bd2JpsiPhi with average angular acceptance input as a paramter
 *
 *  @author Ailsa Sparkes asparkes@cern.ch
 *  @date 2011-01-26
 */

#include "TMath.h"
#include <cmath>

#include "Bd2JpsiKstar_sWave_KpiBins.h"
#include "Mathematics.h"
#include "SlicedAcceptance.h"
#include <iostream>
//#include "math.h"
//#include "TMath.h"

PDF_CREATOR( Bd2JpsiKstar_sWave_KpiBins );

//Constructor
Bd2JpsiKstar_sWave_KpiBins::Bd2JpsiKstar_sWave_KpiBins(PDFConfigurator* configurator) :
	cachedAzeroAzeroIntB(), cachedAparaAparaIntB(), cachedAperpAperpIntB(), cachedAparaAperpIntB(), cachedAzeroAparaIntB(), cachedAzeroAperpIntB(),
	cachedAsAsIntB(), cachedAparaAsIntB(), cachedAperpAsIntB(), cachedAzeroAsIntB(), AzeroAzeroB(), AparaAparaB(), AperpAperpB(), AsAsB(), ImAparaAperpB(),
	ReAzeroAparaB(), ImAzeroAperpB(), ReAparaAsB(), ImAperpAsB(), ReAzeroAsB(), cachedSinDeltaPerpPara(), cachedCosDeltaPara(), cachedSinDeltaPerp(),
	cachedCosDeltaParaS(), cachedSinDeltaPerpS(), cachedCosDeltaS(), cachedAzero(), cachedApara(), cachedAperp(), cachedAs(), timeAcc(NULL),t(),
	// Physics parameters
	gammaName             ( configurator->getName("gamma" ))
	, Rzero_sqName ( configurator->getName("Rzero_sq"))
	, Rpara_sqName  ( configurator->getName("Rpara_sq" ))
	, Rperp_sqName  ( configurator->getName("Rperp_sq" ))
	, As_sqName     ( configurator->getName("As_sq" ))
	, delta_zeroName( configurator->getName("delta_zero" ))
	, delta_paraName( configurator->getName("delta_para" ))
	, delta_perpName( configurator->getName("delta_perp" ))
	, delta_sName   ( configurator->getName("delta_s" ))
	, angAccI1Name  ( configurator->getName("angAccI1" ))
	, angAccI2Name  ( configurator->getName("angAccI2" ))
	, angAccI3Name  ( configurator->getName("angAccI3" ))
	, angAccI4Name  ( configurator->getName("angAccI4" ))
	, angAccI5Name  ( configurator->getName("angAccI5" ))
	, angAccI6Name  ( configurator->getName("angAccI6" ))
	, angAccI7Name  ( configurator->getName("angAccI7" ))
	, angAccI8Name  ( configurator->getName("angAccI8" ))
	, angAccI9Name  ( configurator->getName("angAccI9" ))
	, angAccI10Name ( configurator->getName("angAccI10" ))
	, timeRes1Name  ( configurator->getName("timeResolution1") )
	, timeRes2Name  ( configurator->getName("timeResolution2" ))
	, timeRes1FractionName  ( configurator->getName("timeResolution1Fraction" ))
	, _useTimeAcceptance(false)

	// Observables (What we want to gain from the pdf after inserting physics parameter values)
	, normalisationCacheValid(false)
, evaluationCacheValid(false)
	, timeName      ( configurator->getName("time" ))
	, cosThetaName  ( configurator->getName("cosTheta" ))
	, phiName       ( configurator->getName("phi" ))
	, cosPsiName    ( configurator->getName("cosPsi" ))
	, KstarFlavourName  ( configurator->getName("KstarFlavour" ))
	//, timeres     ( "resolution" )

	, timeconstraintName( configurator->getName("time" ))
, gamma(), Rzero_sq(), Rpara_sq(), Rperp_sq(), As_sq(), AzeroApara(), AzeroAperp(), AparaAperp(), AparaAs(), AperpAs(), AzeroAs(),
	delta_zero(), delta_para(), delta_perp(), delta_s(), omega(), timeRes(), timeRes1(), timeRes2(), timeRes1Frac(), angAccI1(), angAccI2(),
	angAccI3(), angAccI4(), angAccI5(), angAccI6(), angAccI7(), angAccI8(), angAccI9(), angAccI10(), Ap_sq(), Ap(), time(), cosTheta(), phi(),
	cosPsi(), KstarFlavour(), tlo(), thi()

{
	MakePrototypes();
	_useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;

	if( useTimeAcceptance() ) {
		timeAcc = new SlicedAcceptance( 0., 14.0, 0.0171 ) ;
		cout << "Bd2JpsiKstar_sWave_KpiBins:: Constructing timeAcc: Upper time acceptance beta=0.0171 [0 < t < 14] " << endl ;
	}
	else {
		timeAcc = new SlicedAcceptance( 0., 14. ) ;
		cout << "Bd2JpsiKstar_sWave_KpiBins:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl ;
	}

}

//Make the data point and parameter set
void Bd2JpsiKstar_sWave_KpiBins::MakePrototypes()
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
	parameterNames.push_back( Rpara_sqName );
	parameterNames.push_back( Rperp_sqName );
	parameterNames.push_back( As_sqName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
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
Bd2JpsiKstar_sWave_KpiBins::~Bd2JpsiKstar_sWave_KpiBins()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bd2JpsiKstar_sWave_KpiBins::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	// Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
	gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	Rpara_sq   = allParameters.GetPhysicsParameter( Rpara_sqName )->GetValue();
	Rperp_sq   = allParameters.GetPhysicsParameter( Rperp_sqName )->GetValue();
	As_sq   = allParameters.GetPhysicsParameter( As_sqName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	delta_s = allParameters.GetPhysicsParameter( delta_sName )->GetValue();
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

	Rzero_sq = 1 - Rperp_sq - Rpara_sq;
	Aperp_sq = Rperp_sq*(1-As_sq);
	Apara_sq = Rpara_sq*(1-As_sq);
	Azero_sq = Rzero_sq*(1-As_sq);

	AparaAperp = sqrt(Apara_sq)*sqrt(Aperp_sq);
	AzeroApara = sqrt(Azero_sq)*sqrt(Apara_sq);
	AzeroAperp = sqrt(Azero_sq)*sqrt(Aperp_sq);
	AparaAs    = sqrt(Apara_sq)*sqrt(As_sq);
	AperpAs	   = sqrt(Aperp_sq)*sqrt(As_sq);
	AzeroAs    = sqrt(Azero_sq)*sqrt(As_sq);

	return isOK;
}

//Return a list of parameters not to be integrated
vector<string> Bd2JpsiKstar_sWave_KpiBins::GetDoNotIntegrateList()
{
	vector<string> doNotIntList;
	return doNotIntList;
}

double Bd2JpsiKstar_sWave_KpiBins::q() const { return KstarFlavour;}

//Calculate the function value
double Bd2JpsiKstar_sWave_KpiBins::Evaluate(DataPoint * measurement)
{

	double returnValue;
	time = measurement->GetObservable( timeName )->GetValue();
	cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
	phi      = measurement->GetObservable( phiName )->GetValue();
	cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();
	KstarFlavour = measurement->GetObservable( KstarFlavourName )->GetValue();

	//cout << gamma << " " << Aperp_sq << " " << Azero_sq << endl;

	if(timeRes1Frac >= 0.9999)
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		returnValue =  buildPDFnumerator();
	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		double val1 = buildPDFnumerator();
		// Set the member variable for time resolution to the second value and calculate
		timeRes = timeRes2;
		double val2 = buildPDFnumerator();
		//return timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
		returnValue = timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
	}

	if( (returnValue <= 0.) || std::isnan(returnValue) ) {
		cout << " Bd2JpsiKstar_sWave_KpiBins::Evaluate() returns <=0 or nan " << endl ;
		cout << " AT    " << Aperp_sq ;
		cout << " AP    " << Apara_sq ;
		cout << " A0    " << Azero_sq;
		cout << " As   " << As_sq;
		cout << "   Dperp    " << delta_perp;
		cout << "   Dpara    " << delta_para;
		cout << "   Ds     " << delta_s << endl;
		cout << "   gamma   " << gamma << endl;

		throw 10 ;


	}

	return returnValue;
}


double Bd2JpsiKstar_sWave_KpiBins::buildPDFnumerator()
{
	// The angular functions f1->f6 as defined in roadmap Table 1.(same for Kstar)
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
	Mathematics::getBs2JpsiPhiAngularFunctionsWithSwave( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, cosTheta, phi, cosPsi );

	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59  //No tagging so only need 2 (h± pg 72)
	// First for the B
	double AzeroAzeroB, AparaAparaB, AperpAperpB, AsAsB;
	double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB;
	double ReAparaAsB, ImAperpAsB, ReAzeroAsB;

	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB
			, ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
			, AsAsB, ReAparaAsB, ImAperpAsB, ReAzeroAsB
			);

	//q() tags the K* flavour - it changes the sign of f4, f6 and f9
	double v1 = f1 * AzeroAzeroB
		+ f2 * AparaAparaB
		+ f3 * AperpAperpB
		+ f4 * ImAparaAperpB * q()
		+ f5 * ReAzeroAparaB
		+ f6 * ImAzeroAperpB * q()
		+ f7 * AsAsB
		+ f8 * ReAparaAsB
		+ f9 * ImAperpAsB * q()
		+ f10 * ReAzeroAsB
		;
	if( useTimeAcceptance() ) v1  = v1 * timeAcc->getValue(time);
	return v1;
}


double Bd2JpsiKstar_sWave_KpiBins::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	double returnValue;
	IConstraint * timeBound = boundary->GetConstraint(timeconstraintName);

	time = measurement->GetObservable( timeName )->GetValue();
	KstarFlavour = measurement->GetObservable( KstarFlavourName )->GetValue();


	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		return -1.;
	}
	else
	{
		tlo = timeBound->GetMinimum();
		thi = timeBound->GetMaximum();
	}


	if(timeRes1Frac >= 0.9999)
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		returnValue =  buildCompositePDFdenominator();

	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		double val1 = buildCompositePDFdenominator();
		// Set the member variable for time resolution to the second value and calculate
		timeRes = timeRes2;
		double val2 = buildCompositePDFdenominator();
		//return timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
		returnValue = timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;

	}


	if( (returnValue <= 0.) || std::isnan(returnValue) ) {
		cout << " Bd2JpsiKstar_sWave_KpiBins::Normalisation() returns <=0 or nan " << endl ;
		cout << " AT    " << Aperp_sq ;
		cout << " AP    " << Apara_sq ;
		cout << " A0    " << Azero_sq;
		cout << " As   " << As_sq;
		cout << "   Dperp    " << delta_perp;
		cout << "   Dpara    " << delta_para;
		cout << "   Ds     " << delta_s << endl;
		cout << "   gamma   " << gamma << endl;

		exit(1) ;
	}

	return returnValue;
}



//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,

double Bd2JpsiKstar_sWave_KpiBins::buildCompositePDFdenominator( )
{
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;
	double returnValue = 0;

	if( true /*useTimeAcceptance()*/ ) {                // Set to true because seleting false makes a single slice for 0 --> 14.
		//This loops over each time slice, does the normalisation between the limits, and accumulates
		for( unsigned int islice = 0; islice < timeAcc->numberOfSlices(); ++islice )
		{
			//Set the time integrals
			tlo = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow() ;
			thi = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh() ;
			if( thi> tlo ) returnValue+= this->buildPDFdenominator(  ) * timeAcc->getSlice(islice)->height() ;
		}
	}
	else {
		returnValue = this->buildPDFdenominator() ;
	}

	tlo = tlo_boundary;
	thi = thi_boundary ;
	return returnValue ;
}




double Bd2JpsiKstar_sWave_KpiBins::NormAnglesOnlyForAcceptanceWeights(DataPoint* measurement, PhaseSpaceBoundary * boundary)
{
	(void) boundary;
	double returnValue;
	time = measurement->GetObservable( timeName )->GetValue();
	KstarFlavour = measurement->GetObservable( KstarFlavourName )->GetValue();

	//First job for any new set of parameters is to Cache the time integrals


	if(timeRes1Frac >= 0.9999)
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		returnValue =  buildPDFdenominatorAngles();
	}
	else
	{
		// Set the member variable for time resolution to the first value and calculate
		timeRes = timeRes1;
		double val1 = buildPDFdenominatorAngles();
		// Set the member variable for time resolution to the second value and calculate
		timeRes = timeRes2;
		double val2 = buildPDFdenominatorAngles();
		//return timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
		returnValue = timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;

	}


	if( (returnValue <= 0.) || std::isnan(returnValue) ) {
		cout << " Bd2JpsiKstar_sWave_KpiBins::Normalisation() returns <=0 or nan " << endl ;
		cout << " AT    " << Aperp_sq ;
		cout << " AP    " << Apara_sq ;
		cout << " A0    " << Azero_sq;
		cout << " As   " << As_sq;
		cout << "   Dperp    " << delta_perp;
		cout << "   Dpara    " << delta_para;
		cout << "   Ds     " << delta_s << endl;
		cout << "   gamma   " << gamma << endl;

		exit(1) ;
	}

	return returnValue;
}



double Bd2JpsiKstar_sWave_KpiBins::buildPDFdenominator()
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




	double v1 = cachedAzeroAzeroIntB * angAccI1
		+ cachedAparaAparaIntB * angAccI2
		+ cachedAperpAperpIntB * angAccI3
		+ cachedAparaAperpIntB * angAccI4* q()
		+ cachedAzeroAparaIntB * angAccI5
		+ cachedAzeroAperpIntB * angAccI6 * q()
		+ cachedAsAsIntB * angAccI7
		+ cachedAparaAsIntB * angAccI8
		+ cachedAperpAsIntB * angAccI9 * q()
		+ cachedAzeroAsIntB * angAccI10
		;
	return v1;




}

double Bd2JpsiKstar_sWave_KpiBins::buildPDFdenominatorAngles()  //test method
{

	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
	Mathematics::getBs2JpsiPhiAngularFunctionsWithSwave( f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, cosTheta, phi, cosPsi );


	// The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB,
			ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
			, AsAsB, ReAparaAsB, ImAperpAsB, ReAzeroAsB
			);


	double v1 =  AzeroAzeroB
		+ AparaAparaB
		+ AperpAperpB
		+ AsAsB
		;
	return v1;

}


void Bd2JpsiKstar_sWave_KpiBins::getTimeDependentAmplitudes(
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



	// Now calculate the amplitudes

	double Exp = Mathematics::Exp(time, gamma, timeRes);

	AzeroAzero = Azero_sq * Exp;  // changed- see note 2009-015 eq 11-13
	AparaApara = Apara_sq * Exp;  //
	AperpAperp = Aperp_sq * Exp;  //
	AsAs = As_sq * Exp;

	ImAparaAperp = cachedApara*cachedAperp * cachedSinDeltaPerpPara * Exp;    //See http://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=33933 page14
	ReAzeroApara = cachedAzero*cachedApara * cachedCosDeltaPara * Exp;
	ImAzeroAperp = cachedAzero*cachedAperp * cachedSinDeltaPerp * Exp;
	ReAparaAs = cachedApara*cachedAs * cachedCosDeltaParaS * Exp; //AILSA_ NOT SURE
	ImAperpAs = cachedAperp*cachedAs * cachedSinDeltaPerpS * Exp; //AILSA
	ReAzeroAs = cachedAzero*cachedAs * cachedCosDeltaS * Exp; //AILSA


	//if ( std::isnan(ImAparaAperp)) cout << Azero_sq << " " << Apara_sq << " " << Aperp_sq << " " << Exp << endl;

	return;
}

void Bd2JpsiKstar_sWave_KpiBins::getTimeAmplitudeIntegrals(
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

	if ( !evaluationCacheValid )
	{


		cachedAzero = sqrt( Azero_sq );
		cachedApara = sqrt( Apara_sq );
		cachedAperp = sqrt( Aperp_sq );
		cachedAs = sqrt (As_sq);

		cachedSinDeltaPerpPara  = sin( delta_perp - delta_para );
		cachedCosDeltaPara      = cos( delta_para );
		cachedSinDeltaPerp      = sin( delta_perp );
		cachedCosDeltaParaS     = cos( delta_para - delta_s );
		cachedSinDeltaPerpS     = sin( delta_perp - delta_s );
		cachedCosDeltaS         = cos( delta_s);

		evaluationCacheValid = true;
	}




	double ExpInt = Mathematics::ExpInt(tlo, thi, gamma, timeRes);


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
