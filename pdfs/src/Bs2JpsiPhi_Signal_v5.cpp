// $Id: Bs2JpsiPhi_Signal_v5.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_Signal_v5 Bs2JpsiPhi_Signal_v5.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-13
 */

#include "Mathematics.h"
#include "Bs2JpsiPhi_Angluar_Terms.h"
#include "Bs2JpsiPhi_Signal_v5.h"

#include <iostream>
#include <cmath>
// #include "TF1.h"

#define DEBUGFLAG true

using namespace::std;

PDF_CREATOR( Bs2JpsiPhi_Signal_v5 );

//......................................
//Constructor(s)
//New one with configurator
Bs2JpsiPhi_Signal_v5::Bs2JpsiPhi_Signal_v5(PDFConfigurator* configurator) :
	// Physics parameters
	gammaName				( configurator->getName("gamma") )
	, deltaGammaName		( configurator->getName("deltaGamma") )
	, deltaMName			( configurator->getName("deltaM") )
	, Phi_sName				( configurator->getName("Phi_s") )
	, Azero_sqName			( configurator->getName("Azero_sq") )
	, Apara_sqName			( configurator->getName("Apara_sq") )
	, Aperp_sqName			( configurator->getName("Aperp_sq") )
	, delta_zeroName		( configurator->getName("delta_zero") )
	, delta_paraName		( configurator->getName("delta_para") )
	, delta_perpName		( configurator->getName("delta_perp") )
	, As_sqName				( configurator->getName("F_s") )
	, delta_sName			( configurator->getName("delta_s") )
	, CspName				( configurator->getName("Csp") )
	, cosdparName			( configurator->getName("cosdpar") ) //PELC-COSDPAR Special for fitting cosdpar separately
	, cosphisName			( configurator->getName("cosphis") )
	, sinphisName			( configurator->getName("sinphis") )
	, lambdaName			( configurator->getName("lambda") )
	// tagging parameters
	, mistagName			( configurator->getName("mistag") )
	, mistagP1Name			( configurator->getName("mistagP1") )
	, mistagP0Name			( configurator->getName("mistagP0") )
	, mistagSetPointName	( configurator->getName("mistagSetPoint") )
	, mistagDeltaP1Name		( configurator->getName("mistagDeltaP1") )
	, mistagDeltaP0Name		( configurator->getName("mistagDeltaP0") )
	, mistagDeltaSetPointName ( configurator->getName("mistagDeltaSetPoint") )
	// Detector parameters
	, resScaleName			( configurator->getName("timeResolutionScale") )
	, eventResolutionName	( configurator->getName("eventResolution") )
	, res1Name				( configurator->getName("timeResolution1") )
	, res2Name				( configurator->getName("timeResolution2") )
	, res3Name				( configurator->getName("timeResolution3") )
	, res2FractionName		( configurator->getName("timeResolution2Fraction") )
	, res3FractionName		( configurator->getName("timeResolution3Fraction") )
	, timeOffsetName		( configurator->getName("timeOffset") )
	// Angular acceptance factors
	, angAccI1Name			( configurator->getName("angAccI1") )
	, angAccI2Name			( configurator->getName("angAccI2") )
	, angAccI3Name			( configurator->getName("angAccI3") )
	, angAccI4Name			( configurator->getName("angAccI4") )
	, angAccI5Name			( configurator->getName("angAccI5") )
	, angAccI6Name			( configurator->getName("angAccI6") )
	, angAccI7Name			( configurator->getName("angAccI7") )
	, angAccI8Name			( configurator->getName("angAccI8") )
	, angAccI9Name			( configurator->getName("angAccI9") )
	, angAccI10Name			( configurator->getName("angAccI10") )
	// Observables
	, timeName				( configurator->getName("time") )
	, cosThetaName			( configurator->getName("cosTheta") )
	, cosPsiName			( configurator->getName("cosPsi") )
	, phiName				( configurator->getName("phi") )
	, cthetakName 			( configurator->getName("helcosthetaK") )
	, cthetalName			( configurator->getName("helcosthetaL") )
	, phihName				( configurator->getName("helphi") )
	, tagName				( configurator->getName("tag") )
	// Other things
	, _useEventResolution(false)
	, _useTimeAcceptance(false)
	, _useHelicityBasis(false)
	, _numericIntegralForce(false)
	, _numericIntegralTimeOnly(false)
	, _useCosAndSin(false)
	, _useCosDpar(false)
	, _usePunziMistag(false)
	, _usePunziSigmat(false)
	, allowNegativeAsSq(false)
	//objects
	,t(), ctheta_tr(), phi_tr(), ctheta_1(), ctheta_k(), phi_h(), ctheta_l(), tag(),
	_gamma(), dgam(), Aperp_sq(), Apara_sq(), Azero_sq(), As_sq(), delta_para(),
	delta_perp(), delta_zero(), delta_s(), delta1(), delta2(), delta_ms(), phi_s(), _cosphis(), _sinphis(), _mistag(), _mistagP1(), _mistagP0(), _mistagSetPoint(),
	resolutionScale(), resolution1(), resolution2(), resolution3(), resolution2Fraction(), resolution3Fraction(), timeOffset(),
	angAccI1(), angAccI2(), angAccI3(), angAccI4(), angAccI5(), angAccI6(), angAccI7(), angAccI8(), angAccI9(), angAccI10(),
	tlo(), thi(), expL_stored(), expH_stored(), expSin_stored(), expCos_stored(),
	intExpL_stored(), intExpH_stored(), intExpSin_stored(), intExpCos_stored(), timeAcc(NULL), normalisationCacheValid(false),
	CachedA1(), CachedA2(), CachedA3(), CachedA4(), CachedA5(), CachedA6(), CachedA7(), CachedA8(), CachedA9(), CachedA10(),
	resolution(), eventResolution(),timeIntegralCacheValid(), storeExpL(), storeExpH(), storeExpSin(), storeExpCos(), normalisationCacheUntagged()
{
	componentIndex = 0;

	std::cout << "Constructing PDF: Bs2JpsiPhi_Signal_v5 " << std::endl ;

	//...........................................
	// Configure  options
	_numericIntegralForce    = configurator->isTrue( "NumericIntegralForce") ;
	_numericIntegralTimeOnly = configurator->isTrue( "NumericIntegralTimeOnly" ) ;
	_useEventResolution = configurator->isTrue( "UseEventResolution" ) ;
	_useCosAndSin = configurator->isTrue( "UseCosAndSin" ) ;
	_useCosDpar = configurator->isTrue( "UseCosDpar" ) ;
	_useHelicityBasis = configurator->isTrue( "UseHelicityBasis" ) ;
	_usePunziSigmat = configurator->isTrue( "UsePunziSigmat" ) ;
	_usePunziMistag = configurator->isTrue( "UsePunziMistag" ) ;
	allowNegativeAsSq = configurator->isTrue( "AllowNegativeAsSq" ) ;


	//...............................................
	// Configure to use angular acceptance machinery
	string angAccFile = configurator->getConfigurationValue( "AngularAcceptanceFile" ) ;
	_angAccIgnoreNumerator = configurator->isTrue( "AngularAcceptanceIgnoreNumerator" ) ;
	if( angAccFile == "" ) cout << "Bs2JpsiPhi_Signal_v5:: Using flat angular acceptance " << endl ;
	else cout << "Bs2JpsiPhi_Signal_v5:: Constructing angAcc using file: " << angAccFile << endl ;
	angAcc = new AngularAcceptance( angAccFile, _useHelicityBasis ) ;
	angAccI1 = angAcc->af1() ;  cout << "  af1 = " << angAccI1 << endl ;
	angAccI2 = angAcc->af2() ;	cout << "  af2 = " << angAccI2 << endl ;
	angAccI3 = angAcc->af3() ;	cout << "  af3 = " << angAccI3 << endl ;
	angAccI4 = angAcc->af4() ;	cout << "  af4 = " << angAccI4 << endl ;
	angAccI5 = angAcc->af5() ;	cout << "  af5 = " << angAccI5 << endl ;
	angAccI6 = angAcc->af6() ;	cout << "  af6 = " << angAccI6 << endl ;
	angAccI7 = angAcc->af7() ;	cout << "  af7 = " << angAccI7 << endl ;
	angAccI8 = angAcc->af8() ;	cout << "  af8 = " << angAccI8 << endl ;
	angAccI9 = angAcc->af9() ;	cout << "  af9 = " << angAccI9 << endl ;
	angAccI10 = angAcc->af10();	cout << "  af10 = " << angAccI10 << endl ;
	if( _angAccIgnoreNumerator ) cout << "Bs2JpsiPhi_Signal_v5:: Ignoring angular acceptance numerator " << endl ;

	//...........................................
	// Configure to use time acceptance machinery
	_useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;
	if( useTimeAcceptance() ) {
		if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) ) {
			timeAcc = new SlicedAcceptance( 0., 14.0, /*0.0157*/ 0.0112) ;
			cout << "Bs2JpsiPhi_Signal_v5:: Constructing timeAcc: Upper time acceptance beta=0.0112 [0 < t < 14] " << endl ;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" ) {
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ) ) ;
			cout << "Bs2JpsiPhi_Signal_v5:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl ;
		}
	}
	else {
		timeAcc = new SlicedAcceptance( 0., 14. ) ;
		cout << "Bs2JpsiPhi_Signal_v5:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl ;
	}

	//Make empty Cache for the time integrals. This has to be done now after the SlicedAcceptance is created
	vector<double> empty ;
	for( unsigned int islice = 0; islice < timeAcc->numberOfSlices(); ++islice ) empty.push_back(0.0) ;
	for( int ires=0; ires < 4 ; ++ires ) {
		storeExpL.push_back( empty ) ;
		storeExpH.push_back( empty ) ;
		storeExpSin.push_back( empty ) ;
		storeExpCos.push_back( empty ) ;
	}


	this->TurnCachingOff();

	this->SetNumericalNormalisation( false );

	//........................
	// Now do some actual work
	this->MakePrototypes();
	this->SetupAngularTerms();

	//PELC  - debug to plot the distribution of PDF values for each event
	//histOfPdfValues = new TH1D( "HistOfPdfValue" ,  "HistOfPdfValue" , 110, -0.00001, 0.00001 ) ;
	//c0  = new TCanvas;
	//histCounter = 0;
	//~PELC


}

void Bs2JpsiPhi_Signal_v5::SetupAngularTerms()
{
	A0A0_Obs = PseudoObservable( "A0A0" );
	APAP_Obs = PseudoObservable( "APAP" );
	ATAT_Obs = PseudoObservable( "ATAT" );
	ASAS_Obs = PseudoObservable( "ASAS" );
	ImAPAT_Obs = PseudoObservable( "ImAPAT" );
	ReA0AP_Obs = PseudoObservable( "ReA0AP" );
	ImA0AT_Obs = PseudoObservable( "ImA0AT" );
	ReASAP_Obs = PseudoObservable( "ReASAP" );
	ImASAT_Obs = PseudoObservable( "ImASAT" );
	ReASA0_Obs = PseudoObservable( "ReASA0" );

	if( _useHelicityBasis )
	{
		//      Required to let the PseudoObservables for the Angular terms Know what they require
		angularTermDependencies.push_back( cthetakName );
		angularTermDependencies.push_back( cthetalName );
		angularTermDependencies.push_back( phihName );
	}
	else
	{
		//      Required to let the PseudoObservables for the Angular terms Know what they require
		angularTermDependencies.push_back( cosThetaName );
		angularTermDependencies.push_back( cosPsiName );
		angularTermDependencies.push_back( phiName );
	}

	A0A0_Obs.AddDependencies( angularTermDependencies );
	APAP_Obs.AddDependencies( angularTermDependencies );
	ATAT_Obs.AddDependencies( angularTermDependencies );
	ASAS_Obs.AddDependencies( angularTermDependencies );
	ImAPAT_Obs.AddDependencies( angularTermDependencies );
	ReA0AP_Obs.AddDependencies( angularTermDependencies );
	ImA0AT_Obs.AddDependencies( angularTermDependencies );
	ReASAP_Obs.AddDependencies( angularTermDependencies );
	ImASAT_Obs.AddDependencies( angularTermDependencies );
	ReASA0_Obs.AddDependencies( angularTermDependencies );

	if( _useHelicityBasis )
	{
		A0A0_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::HangleFactorA0A0 );
		APAP_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::HangleFactorAPAP );
		ATAT_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::HangleFactorATAT );
		ASAS_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::HangleFactorASAS );
		ImAPAT_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::HangleFactorImAPAT );
		ReA0AP_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::HangleFactorReA0AP );
		ImA0AT_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::HangleFactorImA0AT );
		ReASAP_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::HangleFactorReASAP );
		ImASAT_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::HangleFactorImASAT );
		ReASA0_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::HangleFactorReASA0 );
	}
	else
	{
		A0A0_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::TangleFactorA0A0 );
		APAP_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::TangleFactorAPAP );
		ATAT_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::TangleFactorATAT );
		ASAS_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::TangleFactorASAS );
		ImAPAT_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::TangleFactorImAPAT );
		ReA0AP_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::TangleFactorReA0AP );
		ImA0AT_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::TangleFactorImA0AT );
		ReASAP_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::TangleFactorReASAP );
		ImASAT_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::TangleFactorImASAT );
		ReASA0_Obs.AddFunction( Bs2JpsiPhi_Angular_Terms::TangleFactorReASA0 );
	}

	//	The cached Values per data point will ALWAYS BE VALID
	A0A0_Obs.SetValid( true );
	APAP_Obs.SetValid( true );
	ATAT_Obs.SetValid( true );
	ASAS_Obs.SetValid( true );
	ImAPAT_Obs.SetValid( true );
	ReA0AP_Obs.SetValid( true );
	ImA0AT_Obs.SetValid( true );
	ReASAP_Obs.SetValid( true );
	ImASAT_Obs.SetValid( true );
	ReASA0_Obs.SetValid( true );
}


//........................................................
//Destructor
Bs2JpsiPhi_Signal_v5::~Bs2JpsiPhi_Signal_v5() {}


//......................................
//Make the data point and parameter set
void Bs2JpsiPhi_Signal_v5::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );

	if( _useHelicityBasis ) {
		allObservables.push_back( cthetakName );
		allObservables.push_back( cthetalName );
		allObservables.push_back( phihName );
	}
	else {
		allObservables.push_back( cosThetaName );
		allObservables.push_back( phiName );
		allObservables.push_back( cosPsiName );
	}

	allObservables.push_back( tagName );
	allObservables.push_back( mistagName );

	if(useEventResolution()) allObservables.push_back( eventResolutionName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( deltaGammaName );
	parameterNames.push_back( Aperp_sqName );
	parameterNames.push_back( Azero_sqName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_zeroName );
	parameterNames.push_back( As_sqName );
	parameterNames.push_back( delta_sName );
	parameterNames.push_back( CspName );
	if( _useCosDpar ) parameterNames.push_back( cosdparName ); //PELC-COSDPAR Special for fitting cosdpar separately
	parameterNames.push_back( deltaMName );

	if( _useCosAndSin ) {
		parameterNames.push_back( cosphisName );
		parameterNames.push_back( sinphisName );
	}
	else{
		parameterNames.push_back( Phi_sName );
	}
	parameterNames.push_back( lambdaName );

	parameterNames.push_back( mistagP1Name );
	parameterNames.push_back( mistagP0Name );
	parameterNames.push_back( mistagSetPointName );
	parameterNames.push_back( mistagDeltaP1Name );
	parameterNames.push_back( mistagDeltaP0Name );
	parameterNames.push_back( mistagDeltaSetPointName );

	parameterNames.push_back( resScaleName );
	if( ! useEventResolution() ) {
		parameterNames.push_back( res1Name );
		parameterNames.push_back( res2Name );
		parameterNames.push_back( res3Name );
		parameterNames.push_back( res2FractionName );
		parameterNames.push_back( res3FractionName );
	}
	parameterNames.push_back( timeOffsetName );

	allParameters = ParameterSet(parameterNames);
}


//.........................................................
//Return a list of observables not to be integrated
vector<string> Bs2JpsiPhi_Signal_v5::GetDoNotIntegrateList()
{
	vector<string> list;

	if( ! _usePunziMistag) list.push_back(mistagName) ;
	if( useEventResolution() && ! _usePunziSigmat) list.push_back(eventResolutionName) ;

	if( _numericIntegralTimeOnly ) {
		if( _useHelicityBasis ) {
			list.push_back( cthetakName );
			list.push_back( cthetalName ) ;
			list.push_back( phihName ) ;
		}
		else {
			list.push_back( cosThetaName );
			list.push_back( cosPsiName ) ;
			list.push_back( phiName ) ;
		}
	}
	return list;
}


//........................................................
//Set the physics parameters into member variables

bool Bs2JpsiPhi_Signal_v5::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;  //This is only used for the untagged events and only if not useing event resolution
	timeIntegralCacheValid = false;   //This cannot be used if event resolution is used

	bool result = allParameters.SetPhysicsParameters(NewParameterSet);

	// Physics parameters.
	_gamma  = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	dgam      = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();

	Azero_sq = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
	if( (Azero_sq < 0.) || (Azero_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_Signal_v5::SetPhysicsParameters: Azero_sq <0 or >1 but left as is" <<  endl ;	}
	Aperp_sq = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
	if( (Aperp_sq < 0.) || (Aperp_sq > 1.)  ) { cout << "Warning in Bs2JpsiPhi_Signal_v5::SetPhysicsParameters: Aperp_sq <0 or >1 but left as is" <<  endl ;	}

	Apara_sq = (1. - Azero_sq - Aperp_sq ) ;
	if( Apara_sq < 0. ) {
		cout << "Warning in Bs2JpsiPhi_Signal_v5::SetPhysicsParameters: derived parameter Apara_sq <0  and so set to zero" <<  endl ;
		Apara_sq = 0. ;
	}

	double fs = allParameters.GetPhysicsParameter( As_sqName )->GetValue();
	if( (fs < 0.) || (fs >= 1.) ) { cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_dev::SetPhysicsParameters: As_sq <0 or >=1 but left as is" <<  endl ;	}
	As_sq = fs / (1. - fs ) ;

	Csp = allParameters.GetPhysicsParameter( CspName )->GetValue();

	delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	//	delta_s	   = allParameters.GetPhysicsParameter( delta_sName )->GetValue();  /// This is the original
	delta_s = allParameters.GetPhysicsParameter( delta_sName )->GetValue() +delta_perp ;   // This is the ambiguity inspired one with less corrn to delta_perp
	delta1 = delta_perp -  delta_para ;
	delta2 = delta_perp -  delta_zero ;

	if( _useCosDpar ) cosdpar = allParameters.GetPhysicsParameter( cosdparName )->GetValue(); //PELC-COSDPAR Special for fitting cosdpar separately

	delta_ms = allParameters.GetPhysicsParameter( deltaMName )->GetValue();

	if(_useCosAndSin){
		_cosphis = allParameters.GetPhysicsParameter( cosphisName )->GetValue();
		_sinphis = allParameters.GetPhysicsParameter( sinphisName )->GetValue();
	}
	else{
		phi_s     = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
		_cosphis = cos(phi_s) ;
		_sinphis = sin(phi_s) ;
	}
	lambda = allParameters.GetPhysicsParameter( lambdaName )->GetValue();

	// Mistag parameters
	_mistagP1		= allParameters.GetPhysicsParameter( mistagP1Name )->GetValue();
	_mistagP0		= allParameters.GetPhysicsParameter( mistagP0Name )->GetValue();
	_mistagSetPoint = allParameters.GetPhysicsParameter( mistagSetPointName )->GetValue();
	_mistagDeltaP1		= allParameters.GetPhysicsParameter( mistagDeltaP1Name )->GetValue();
	_mistagDeltaP0		= allParameters.GetPhysicsParameter( mistagDeltaP0Name )->GetValue();
	_mistagDeltaSetPoint = allParameters.GetPhysicsParameter( mistagDeltaSetPointName )->GetValue();

	// Detector parameters
	resolutionScale		= allParameters.GetPhysicsParameter( resScaleName )->GetValue();
	if( ! useEventResolution() ) {
		resolution1         = allParameters.GetPhysicsParameter( res1Name )->GetValue();
		resolution2         = allParameters.GetPhysicsParameter( res2Name )->GetValue();
		resolution3         = allParameters.GetPhysicsParameter( res3Name )->GetValue();
		resolution2Fraction = allParameters.GetPhysicsParameter( res2FractionName )->GetValue();
		resolution3Fraction = allParameters.GetPhysicsParameter( res3FractionName )->GetValue();
	}
	timeOffset          = allParameters.GetPhysicsParameter( timeOffsetName )->GetValue();

	// New: Prepare the coefficients of all of the time dependent terms (C,D,S etc)
	this->prepareCDS() ;
	this->prepareTimeFac();
	return result;
}

void Bs2JpsiPhi_Signal_v5::prepareTimeFac()
{
	//	Evaluate time functions
	_expLObs = PseudoObservable( "expL_time" );
	_expLObs.AddFunction( Mathematics::Exp_Wrapper );
	_expHObs = PseudoObservable( "expH_time" );
	_expHObs.AddFunction( Mathematics::Exp_Wrapper );
	_expSinObs = PseudoObservable( "expSin_time" );
	_expSinObs.AddFunction( Mathematics::ExpSin_Wrapper );
	_expCosObs = PseudoObservable( "expCos_time" );
	_expCosObs.AddFunction( Mathematics::ExpCos_Wrapper );

	//	Normalise time functions
	if( !useTimeAcceptance() )
	{
		_intexpLObs = PseudoObservable( "intexpL_time" );
		_intexpLObs.AddFunction( Mathematics::ExpInt_Wrapper );
		_intexpHObs = PseudoObservable( "intexpH_time" );
		_intexpHObs.AddFunction( Mathematics::ExpInt_Wrapper );
		_intexpSinObs = PseudoObservable( "intexpSin_time" );
		_intexpSinObs.AddFunction( Mathematics::ExpSinInt_Wrapper );
		_intexpCosObs = PseudoObservable( "intexpCos_time" );
		_intexpCosObs.AddFunction( Mathematics::ExpCosInt_Wrapper );
	}
	/*
	   else
	   {

	// This Code DOES WORK, it just ends up consuming a HUGE amount of memory and as such ends up being a worse bottleneck than the one it was intended to remove

	for( unsigned int i=0; i< timeAcc->numberOfSlices(); ++i )
	{
	TString _intexpLObs_vec_name="intexpL_time_"; _intexpLObs_vec_name+=i;
	_intexpLObs_vec.push_back( PseudoObservable( _intexpLObs_vec_name.Data() ) );
	_intexpLObs_vec.back().AddFunction( Mathematics::ExpInt_Wrapper );
	TString _intexpHObs_vec_name="intexpH_time"; _intexpHObs_vec_name+=i;
	_intexpHObs_vec.push_back( PseudoObservable( _intexpHObs_vec_name.Data() ) );
	_intexpHObs_vec.back().AddFunction( Mathematics::ExpInt_Wrapper );
	TString _intexpSinObs_vec_name="intexpSin_time"; _intexpSinObs_vec_name+=i;
	_intexpSinObs_vec.push_back( PseudoObservable( _intexpSinObs_vec_name.Data() ) );
	_intexpSinObs_vec.back().AddFunction( Mathematics::ExpSinInt_Wrapper );
	TString _intexpCosObs_vec_name="intexpCos_time"; _intexpCosObs_vec_name+=i;
	_intexpCosObs_vec.push_back( PseudoObservable( _intexpCosObs_vec_name.Data() ) );
	_intexpCosObs_vec.back().AddFunction( Mathematics::ExpCosInt_Wrapper );
	}

	}
	*/
}


//.............................................................
//Calculate the PDF value for a given set of observables for use by numeric integral

double Bs2JpsiPhi_Signal_v5::EvaluateForNumericIntegral(DataPoint * measurement)
{
	if( _numericIntegralTimeOnly ) return this->EvaluateTimeOnly(measurement) ;
	else return this->Evaluate(measurement) ;
}


//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_Signal_v5::Evaluate(DataPoint * measurement)
{
	_datapoint = measurement;

	A0A0_value = measurement->GetPseudoObservable( A0A0_Obs );
	APAP_value = measurement->GetPseudoObservable( APAP_Obs );
	ATAT_value = measurement->GetPseudoObservable( ATAT_Obs );
	ASAS_value = measurement->GetPseudoObservable( ASAS_Obs );
	ImAPAT_value = measurement->GetPseudoObservable( ImAPAT_Obs );
	ReA0AP_value = measurement->GetPseudoObservable( ReA0AP_Obs );
	ImA0AT_value = measurement->GetPseudoObservable( ImA0AT_Obs );
	ReASAP_value = measurement->GetPseudoObservable( ReASAP_Obs );
	ImASAT_value = measurement->GetPseudoObservable( ImASAT_Obs );
	ReASA0_value = measurement->GetPseudoObservable( ReASA0_Obs );

	double angAcceptanceFactor = 0 ;

	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;

	if( _useHelicityBasis )
	{
		Observable* thetaK_obs = measurement->GetObservable( cthetakName );
		Observable* thetaL_obs = measurement->GetObservable( cthetalName );
		Observable* hphi_obs = measurement->GetObservable( phihName );
		ctheta_k   = thetaK_obs->GetValue();
		phi_h      = hphi_obs->GetValue();
		ctheta_l   = thetaL_obs->GetValue();
		angAcceptanceFactor = angAcc->getValue( thetaK_obs, thetaL_obs, hphi_obs );  // Histogram is generated in PDF basis!
	}
	else
	{
		Observable* theta_obs = measurement->GetObservable( cosThetaName );
		Observable* psi_obs = measurement->GetObservable( cosPsiName );
		Observable* phi_obs = measurement->GetObservable( phiName );
		ctheta_tr = theta_obs->GetValue();
		phi_tr    = phi_obs->GetValue();
		ctheta_1  = psi_obs->GetValue();
		angAcceptanceFactor = angAcc->getValue( psi_obs, theta_obs, phi_obs );
	}

	tag = (int)measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue();

	if( useEventResolution() ) eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();

	//Cache amplitues and angles terms used in cross section
	this->CacheAmplitudesAndAngles() ;

	//***** This will be returned*****
	//It gets constructed according to what you want to do with resolution
	double returnValue ;


	if( resolutionScale <= 0. ) {
		//This is the "code" to run with resolution=0
		resolution = 0. ;
		returnValue = this->diffXsec( );
	}
	else if( useEventResolution() ) {
		// Event-by-event resolution has been selected
		resolution = eventResolution * resolutionScale ;
		returnValue = this->diffXsec( );
	}
	else {
		//Good old fashioned triple Gaussian resolution is selected
		double val1=0. , val2=0., val3=0. ;
		double resolution1Fraction = 1. - resolution2Fraction - resolution3Fraction ;
		if(resolution1Fraction > 0 ) {
			resolution = resolution1 * resolutionScale ;
			val1 = this->diffXsec( );
		}
		if(resolution2Fraction > 0 ) {
			resolution = resolution2 * resolutionScale ;
			val2 = this->diffXsec( );
		}
		if(resolution3Fraction > 0 ) {
			resolution = resolution3 * resolutionScale ;
			val3 = this->diffXsec( );
		}
		returnValue = resolution1Fraction*val1 + resolution2Fraction*val2 + resolution3Fraction*val3 ;
	}


	//conditions to throw exception
	bool c1 = isnan(returnValue) ;
	bool c2 = (resolutionScale> 0.) && (returnValue <= 0.) ;
	bool c3 = (resolutionScale<=0.) && (t>0.) && (returnValue <= 0.)  ;
	if( DEBUGFLAG && (c1 || c2 || c3)  ) {
		this->DebugPrint( " Bs2JpsiPhi_Signal_v5::Evaluate() returns <=0 or nan :" , returnValue ) ;
		if( isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}


	if( _angAccIgnoreNumerator ) return returnValue ;
	else return returnValue  * angAcceptanceFactor ;
}


//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_Signal_v5::EvaluateTimeOnly(DataPoint * measurement)
{
	_datapoint = measurement;

	A0A0_value = measurement->GetPseudoObservable( A0A0_Obs );
	APAP_value = measurement->GetPseudoObservable( APAP_Obs );
	ATAT_value = measurement->GetPseudoObservable( ATAT_Obs );
	ASAS_value = measurement->GetPseudoObservable( ASAS_Obs );
	ImAPAT_value = measurement->GetPseudoObservable( ImAPAT_Obs );
	ReA0AP_value = measurement->GetPseudoObservable( ReA0AP_Obs );
	ImA0AT_value = measurement->GetPseudoObservable( ImA0AT_Obs );
	ReASAP_value = measurement->GetPseudoObservable( ReASAP_Obs );
	ImASAT_value = measurement->GetPseudoObservable( ImASAT_Obs );
	ReASA0_value = measurement->GetPseudoObservable( ReASA0_Obs );

	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue();
	if( useEventResolution() ) eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();

	double returnValue;


	if( resolutionScale <= 0. ) {
		resolution = 0. ;
		returnValue = this->diffXsecTimeOnly( );
	}
	else if( useEventResolution() ) {
		resolution = eventResolution * resolutionScale ;
		returnValue = this->diffXsecTimeOnly( );
	}
	else {
		double val1=0. , val2=0., val3=0. ;
		double resolution1Fraction = 1. - resolution2Fraction - resolution3Fraction ;
		if(resolution1Fraction > 0 ) {
			resolution = resolution1 * resolutionScale ;
			val1 = this->diffXsecTimeOnly( );
		}
		if(resolution2Fraction > 0 ) {
			resolution = resolution2 * resolutionScale ;
			val2 = this->diffXsecTimeOnly( );
		}
		if(resolution3Fraction > 0 ) {
			resolution = resolution3 * resolutionScale ;
			val3 = this->diffXsecTimeOnly( );
		}
		returnValue = resolution1Fraction*val1 + resolution2Fraction*val2 + resolution3Fraction*val3 ;
	}


	//conditions to throw exception
	bool c1 = isnan(returnValue) ;
	bool c2 = (resolutionScale> 0.) && (returnValue <= 0.) ;
	bool c3 = (resolutionScale<=0.) && (t>0.) && (returnValue <= 0.)  ;
	if( DEBUGFLAG && (c1 || c2 || c3)  ) {
		this->DebugPrint( " Bs2JpsiPhi_Signal_v5::EvaluateTimeOnly() returns <=0 or nan :" , returnValue ) ;
		if( isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}


	return returnValue ;

}


//...............................................................
//Calculate the normalisation for a given set of physics parameters and boundary

double Bs2JpsiPhi_Signal_v5::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	_datapoint = measurement;

	if( _numericIntegralForce ) return -1. ;

	// Get observables into member variables
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue() ;
	if( useEventResolution() ) eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();

	// Get time boundaries into member variables
	IConstraint * timeBound = boundary->GetConstraint( timeName );
	if ( timeBound->GetUnit() == "NameNotFoundError" ) {
		cerr << "Bound on time not provided" << endl;
		return 0;
	}
	else {
		tlo = timeBound->GetMinimum();
		thi = timeBound->GetMaximum();
	}


	//*** This is what will be returned.***
	//How it is calcualted depends upon how resolution is treated
	double returnValue=0 ;


	// If we are going to use event-by-event resolution, then all the caching is irrelevant and will be bypassed.
	// You cant cache either untagged or tagged since the resolution will change from event to event and affect the normalisation.
	if( useEventResolution() )
	{
		if( resolutionScale <=0. ) resolution = 0. ;
		else resolution = eventResolution * resolutionScale ;
		returnValue = this->diffXsecCompositeNorm1( 0 );
	}
	else //We are not going to use event by event resolution so we can use all the caching machinery
	{

		//First job for any new set of parameters is to Cache the time integrals
		if( ! timeIntegralCacheValid )
		{
			CacheTimeIntegrals() ;
			timeIntegralCacheValid = true ;
			// if( ! useEventResolution() ) timeIntegralCacheValid = true ;
		}


		//If this is an untagged event and the result has been cached, then it can be used
		// Otherwise must calculate the normalisation
		if( (tag==0) && normalisationCacheValid )
		{
			returnValue = normalisationCacheUntagged ;
		}
		else //So we need to calculate the normalisation
		{
			if( resolutionScale <= 0. )
			{
				resolution = 0. ;
				returnValue = this->diffXsecCompositeNorm1( 0 );
			}
			else
			{
				double val1=0. , val2=0., val3=0. ;
				double resolution1Fraction = 1. - resolution2Fraction - resolution3Fraction ;
				if(resolution1Fraction > 0 )
				{
					resolution = resolution1 * resolutionScale ;
					val1 = this->diffXsecCompositeNorm1( 1 );
				}
				if(resolution2Fraction > 0 )
				{
					resolution = resolution2 * resolutionScale ;
					val2 = this->diffXsecCompositeNorm1( 2 );
				}
				if(resolution3Fraction > 0 )
				{
					resolution = resolution3 * resolutionScale ;
					val3 = this->diffXsecCompositeNorm1( 3 );
				}
				returnValue = resolution1Fraction*val1 + resolution2Fraction*val2 + resolution3Fraction*val3 ;
			}
		}

		// If this is an untagged event then the normaisation is invariant and so can be cached
		if( (tag==0) && !normalisationCacheValid )
		{
			normalisationCacheUntagged = returnValue ;
			normalisationCacheValid = true ;
		}

	}

	// Conditions to throw exception
	bool c1 = isnan(returnValue)  ;
	bool c2 = (returnValue <= 0.) ;
	if( DEBUGFLAG && (c1 || c2 ) ) {
		this->DebugPrint( " Bs2JpsiPhi_Signal_v5::Normalisation() returns <=0 or nan :" , returnValue ) ;
		if( isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}


	return returnValue ;
}



//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_Signal_v5::preCalculateTimeFactors()
{
	if( useEventResolution() )
	{
		vector<double> Exp_Input_1(3,0.); Exp_Input_1[0]=t; Exp_Input_1[1]=gamma_l(); Exp_Input_1[2]=resolution;
		vector<double> Exp_Input_2(3,0.); Exp_Input_2[0]=t; Exp_Input_2[1]=gamma_h(); Exp_Input_2[2]=resolution;
		vector<double> Input_2(4,0.); Input_2[0]=t; Input_2[1]=gamma(); Input_2[2]=delta_ms; Input_2[3]=resolution;

		expL_stored = _datapoint->GetPseudoObservable( _expLObs, Exp_Input_1 );
		expH_stored = _datapoint->GetPseudoObservable( _expHObs, Exp_Input_2 );
		expSin_stored = _datapoint->GetPseudoObservable( _expSinObs, Input_2 );
		expCos_stored = _datapoint->GetPseudoObservable( _expCosObs, Input_2 );
	}
	else
	{
		expL_stored = Mathematics::Exp( t, gamma_l(), resolution );
		expH_stored = Mathematics::Exp( t, gamma_h(), resolution );
		expSin_stored = Mathematics::ExpSin( t, gamma(), delta_ms, resolution );
		expCos_stored = Mathematics::ExpCos( t, gamma(), delta_ms, resolution );
	}
	return;
}


//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_Signal_v5::preCalculateTimeIntegrals()
{
	if( !useTimeAcceptance() && useEventResolution() )
	{
		vector<double> Exp_Input_1(4,0.); Exp_Input_1[0]=tlo; Exp_Input_1[1]=thi; Exp_Input_1[2]=gamma_l(); Exp_Input_1[3]=resolution;
		vector<double> Exp_Input_2(4,0.); Exp_Input_2[0]=tlo; Exp_Input_2[1]=thi; Exp_Input_2[2]=gamma_h(); Exp_Input_2[3]=resolution;
		vector<double> Input_2(5,0.); Input_2[0]=tlo; Input_2[1]=thi; Input_2[2]=gamma(); Input_2[3]=delta_ms; Input_2[4]=resolution;

		intExpL_stored = _datapoint->GetPseudoObservable( _intexpLObs, Exp_Input_1 );
		intExpH_stored = _datapoint->GetPseudoObservable( _intexpHObs, Exp_Input_2 );
		intExpSin_stored = _datapoint->GetPseudoObservable( _intexpSinObs, Input_2 );
		intExpCos_stored = _datapoint->GetPseudoObservable( _intexpCosObs, Input_2 );
	}
	else
	{
		intExpL_stored = Mathematics::ExpInt( tlo, thi, gamma_l(), resolution );
		intExpH_stored = Mathematics::ExpInt( tlo, thi, gamma_h(), resolution );
		intExpSin_stored = Mathematics::ExpSinInt( tlo, thi, gamma(), delta_ms, resolution );
		intExpCos_stored = Mathematics::ExpCosInt( tlo, thi, gamma(), delta_ms, resolution );
	}


	/*else if( useTimeAcceptance() && useEventResolution() )
	  {
	//	This works, it just ends up storing the result and input for 160 calculations per datapoint. This requires a HUGE amount of i/o between processor and RAM
	//	even on machines where this is good this isn't as fast as simply lining up the calculation on the CPU.
	//	It also drinks in excess of 1Gb of RAM in it's current form (this isn't likely to reduce) so we can't be clever in the case of timeAcceptance && eventResolution

	vector<double> Exp_Input_1(4,0.); Exp_Input_1[0]=tlo; Exp_Input_1[1]=thi; Exp_Input_1[2]=gamma_l(); Exp_Input_1[3]=resolution;
	vector<double> Exp_Input_2(4,0.); Exp_Input_2[0]=tlo; Exp_Input_2[1]=thi; Exp_Input_2[2]=gamma_h(); Exp_Input_2[3]=resolution;
	vector<double> Input_2(5,0.); Input_2[0]=tlo; Input_2[1]=thi; Input_2[2]=gamma(); Input_2[3]=delta_ms; Input_2[4]=resolution;

	intExpL_stored = _datapoint->GetPseudoObservable( _intexpLObs_vec[ timeBinNum ], Exp_Input_1 );
	intExpH_stored = _datapoint->GetPseudoObservable( _intexpHObs_vec[ timeBinNum ], Exp_Input_2 );
	intExpSin_stored = _datapoint->GetPseudoObservable( _intexpSinObs_vec[ timeBinNum ], Input_2 );
	intExpCos_stored = _datapoint->GetPseudoObservable( _intexpCosObs_vec[ timeBinNum ], Input_2 );
	}*/
	return;
}

vector<string> Bs2JpsiPhi_Signal_v5::PDFComponents()
{
	vector<string> component_list;
	component_list.push_back( "CP-Even" );
	component_list.push_back( "CP-Odd" );
	component_list.push_back( "As" );
	component_list.push_back( "0" );
	return component_list;
}

double Bs2JpsiPhi_Signal_v5::EvaluateComponent( DataPoint* input, ComponentRef* Component )
{
	componentIndex = Component->getComponentNumber();
	if( componentIndex == -1 )
	{
		string ComponentName = Component->getComponentName();
		if( ComponentName.compare( "CP-Even" ) == 0 )
		{
			Component->setComponentNumber( 1 );
			componentIndex = 1;
		}
		else if( ComponentName.compare( "CP-Odd" ) == 0 )
		{
			Component->setComponentNumber( 2 );
			componentIndex = 2;
		}
		else if( ComponentName.compare( "As" ) == 0 )
		{
			Component->setComponentNumber( 3 );
			componentIndex = 3;
		}
		else
		{
			Component->setComponentNumber( 0 );
			componentIndex = 0;
		}
	}

	double return_value = this->Evaluate( input );
	componentIndex = 0;

	return return_value;
}

//...................................
// Main Diff cross section

double Bs2JpsiPhi_Signal_v5::diffXsec()
{
	preCalculateTimeFactors();

	double xsec=-1.;
	switch( componentIndex )
	{
		case 1:		//	CP-Even		CP-Odd=0 && S-Wave=0
			xsec = CachedA1 * timeFactorA0A0(  );
			xsec += CachedA2 * timeFactorAPAP(  );
			xsec += CachedA5 * timeFactorReA0AP(  );
			break;
		case 2:		//	CP-Odd		CP-Even=0 && S-Wave=0
			xsec = CachedA3 * timeFactorATAT(  );
			break;
		case 3:		//	S-Wave		CP-Even=0 && CP-Odd=0
			xsec = CachedA7 * timeFactorASAS(  );
			break;
		default:	//	Everything
			xsec =

				/*
				   A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) +
				   AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) +
				   AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) +

				   AP()*AT() * timeFactorImAPAT(  ) * angleFactorImAPAT( ) +
				   A0()*AP() * timeFactorReA0AP(  ) * angleFactorReA0AP( ) +
				   A0()*AT() * timeFactorImA0AT(  ) * angleFactorImA0AT( ) +

				   AS()*AS() * timeFactorASAS(  ) * angleFactorASAS( ) +

				   AS()*AP() * timeFactorReASAP(  ) * angleFactorReASAP( ) +
				   AS()*AT() * timeFactorImASAT(  ) * angleFactorImASAT( ) +
				   AS()*A0() * timeFactorReASA0(  ) * angleFactorReASA0( ) ;
				   */

				CachedA1 * timeFactorA0A0(  ) +
				CachedA2 * timeFactorAPAP(  ) +
				CachedA3 * timeFactorATAT(  ) +

				CachedA4 * timeFactorImAPAT(  ) +
				CachedA5 * timeFactorReA0AP(  ) +
				CachedA6 * timeFactorImA0AT(  ) +

				CachedA7 * timeFactorASAS(  ) +

				CachedA8 * timeFactorReASAP(  ) +
				CachedA9 * timeFactorImASAT(  ) +
				CachedA10 * timeFactorReASA0(  );

			Observable* timeObs = _datapoint->GetObservable( timeName );
			if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue( timeObs, timeOffset );
			if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2JpsiPhi_Signal_v5_v1::diffXsec( ) : return value < 0 = ", xsec ) ;

			//PELC - This turned out to be an important debugging tool
			//switch it on to see the values of PDF being returend.  If ANY go negative, it means there is a sign wrong in one or more of the terms
			//You need to enable in the .h file as well
			//histOfPdfValues->Fill(xsec) ;
			//histCounter++ ;
			//if( histCounter > 10000 ) {
			//	histOfPdfValues->Draw() ;
			//	c0->Update() ;
			//	c0->SaveAs( "histOfPdfValues-from-Evaluate.eps" ) ;
			//	histCounter = 0 ;
			//}
			break;
	}

	return xsec;
}


//...................................
// Integral over angles only for a fixed time.

double Bs2JpsiPhi_Signal_v5::diffXsecTimeOnly()
{
	preCalculateTimeFactors() ;

	double xsec =

		A0()*A0() * timeFactorA0A0(  ) * angAccI1 +
		AP()*AP() * timeFactorAPAP(  ) * angAccI2 +
		AT()*AT() * timeFactorATAT(  ) * angAccI3 +

		AP()*AT() * timeFactorImAPAT(  ) * angAccI4 +
		A0()*AP() * timeFactorReA0AP(  ) * angAccI5 +
		A0()*AT() * timeFactorImA0AT(  ) * angAccI6 +

		AS()*AS() * timeFactorASAS(  ) * angAccI7 +

		//AS()*AP() * timeFactorReASAP(  ) * angAccI8 +
		//AS()*AT() * timeFactorImASAT(  ) * angAccI9 +
		//AS()*A0() * timeFactorReASA0(  ) * angAccI10 ;
		ASint()*AP() * timeFactorReASAP(  ) * angAccI8 +
		ASint()*AT() * timeFactorImASAT(  ) * angAccI9 +
		ASint()*A0() * timeFactorReASA0(  ) * angAccI10 ;

	Observable* timeObs = _datapoint->GetObservable( timeName );
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue( timeObs, timeOffset );

	if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2JpsiPhi_Signal_v5_v1::diffXsecTimeOnly( ) : return value < 0 = ", xsec ) ;

	return xsec ;
}




//...................................
// Integral over all variables: t + angles

double Bs2JpsiPhi_Signal_v5::diffXsecNorm1()
{
	//preCalculateTimeIntegrals() ;  Replaced by new Caching mechanism , but this cant be used when event resolution is selected
	if( useEventResolution() ) preCalculateTimeIntegrals() ;

	double norm =

		A0()*A0() * timeFactorA0A0Int(  ) * angAccI1   +
		AP()*AP() * timeFactorAPAPInt(  ) * angAccI2   +
		AT()*AT() * timeFactorATATInt(  ) * angAccI3   +

		AP()*AT() * timeFactorImAPATInt(  ) * angAccI4 +
		A0()*AP() * timeFactorReA0APInt(  ) * angAccI5 +
		A0()*AT() * timeFactorImA0ATInt(  ) * angAccI6 +

		AS()*AS() * timeFactorASASInt(  ) * angAccI7   +

		//AS()*AP() * timeFactorReASAPInt(  ) * angAccI8 +
		//AS()*AT() * timeFactorImASATInt(  ) * angAccI9 +
		//AS()*A0() * timeFactorReASA0Int(  ) * angAccI10 ;
		ASint()*AP() * timeFactorReASAPInt(  ) * angAccI8 +
		ASint()*AT() * timeFactorImASATInt(  ) * angAccI9 +
		ASint()*A0() * timeFactorReASA0Int(  ) * angAccI10 ;

	if( DEBUGFLAG && (norm < 0) ) this->DebugPrintNorm( " Bs2JpsiPhi_Signal_v5_v1::diffXsecNorm1( ) : return value < 0 = ", norm ) ;

	return norm ;
}



//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,

double Bs2JpsiPhi_Signal_v5::diffXsecCompositeNorm1( int resolutionIndex )
{
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;
	double returnValue = 0;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		timeBinNum = islice;
		//De cache the time integrals  (unles using event Resolution
		/*if( ! useEventResolution() )*/ this->deCacheTimeIntegrals( (unsigned)resolutionIndex, islice ) ;

		tlo = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow() ;
		thi = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh() ;
		if( thi > tlo ) returnValue+= this->diffXsecNorm1(  ) * timeAcc->getSlice(islice)->height() ;
	}

	tlo = tlo_boundary;
	thi = thi_boundary;
	return returnValue;
}


//.......................................................
// New speed up method to Cache time integrals
void Bs2JpsiPhi_Signal_v5::CacheAmplitudesAndAngles()
{

	CachedA1 = A0()*A0() * A0A0_value;
	CachedA2 = AP()*AP() * APAP_value;
	CachedA3 = AT()*AT() * ATAT_value;

	CachedA4 = AP()*AT() * ImAPAT_value;
	CachedA5 = A0()*AP() * ReA0AP_value;
	CachedA6 = A0()*AT() * ImA0AT_value;

	CachedA7 = AS()*AS() * ASAS_value;

	CachedA8 = ASint()*AP() * ReASAP_value;
	CachedA9 = ASint()*AT() * ImASAT_value;
	CachedA10= ASint()*A0() * ReASA0_value;

	/*
	   CachedA1 = A0()*A0() * angleFactorA0A0();
	   CachedA2 = AP()*AP() * angleFactorAPAP();
	   CachedA3 = AT()*AT() * angleFactorATAT();

	   CachedA4 = AP()*AT() * angleFactorImAPAT();
	   CachedA5 = A0()*AP() * angleFactorReA0AP();
	   CachedA6 = A0()*AT() * angleFactorImA0AT();

	   CachedA7 = AS()*AS() * angleFactorASAS();

	   CachedA8 = AS()*AP() * angleFactorReASAP();
	   CachedA9 = AS()*AT() * angleFactorImASAT();
	   CachedA10= AS()*A0() * angleFactorReASA0();
	   CachedA8 = ASint()*AP() * angleFactorReASAP();
	   CachedA9 = ASint()*AT() * angleFactorImASAT();
	   CachedA10= ASint()*A0() * angleFactorReASA0();
	   */
}

//.......................................................
// New speed up method to Cache time integrals
void Bs2JpsiPhi_Signal_v5::CacheTimeIntegrals() {

	// This need to know  (and be modified for)
	//  --> Number of resolutions
	//  --> Number fo slices
	//  --->  tlo and thi for each slice

	double tlo_boundary = tlo ;
	double thi_boundary = thi ;

	if( useEventResolution() ) {
		unsigned int ires = 0 ;
		resolution = eventResolution * resolutionScale ;
		for( unsigned int islice = 0; islice < (unsigned)timeAcc->numberOfSlices(); ++islice ) {
			tlo = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow() ;
			thi = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh() ;
			if( thi > tlo ) {
				this->preCalculateTimeIntegrals() ;
				//cout << " >>>>> caching time integrals / " << intExpL_stored << "  /  "<< intExpH_stored << "  /  "<< intExpSin_stored << "  /  "<< intExpCos_stored << "  /  " << endl ;

				storeExpL[ires][islice] = intExpL_stored ;
				storeExpH[ires][islice] = intExpH_stored ;
				storeExpSin[ires][islice] = intExpSin_stored ;
				storeExpCos[ires][islice] = intExpCos_stored ;
			}
			else {
				storeExpL[ires][islice] = 0 ;
				storeExpH[ires][islice] = 0 ;
				storeExpSin[ires][islice] = 0 ;
				storeExpCos[ires][islice] = 0 ;
			}

		}
	}

	else {
		for( unsigned int ires=0; ires < 4 ; ++ires ) {

			if( ires==0 ) resolution = 0.0 ;
			if( ires==1 ) resolution = resolution1 * resolutionScale ;
			if( ires==2 ) resolution = resolution2 * resolutionScale ;
			if( ires==3 ) resolution = resolution3 * resolutionScale ;

			for( unsigned int islice = 0; islice < (unsigned)timeAcc->numberOfSlices(); ++islice ) {
				tlo = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow() ;
				thi = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh() ;
				if( thi > tlo ) {
					this->preCalculateTimeIntegrals() ;
					//cout << " >>>>> caching time integrals / " << intExpL_stored << "  /  "<< intExpH_stored << "  /  "<< intExpSin_stored << "  /  "<< intExpCos_stored << "  /  " << endl ;

					storeExpL[ires][islice] = intExpL_stored ;
					storeExpH[ires][islice] = intExpH_stored ;
					storeExpSin[ires][islice] = intExpSin_stored ;
					storeExpCos[ires][islice] = intExpCos_stored ;
				}
				else {
					storeExpL[ires][islice] = 0 ;
					storeExpH[ires][islice] = 0 ;
					storeExpSin[ires][islice] = 0 ;
					storeExpCos[ires][islice] = 0 ;
				}

			}
		}
	}

	tlo = tlo_boundary;
	thi = thi_boundary ;


}


//.......................................................
// New speed up method to Cache time integrals
void Bs2JpsiPhi_Signal_v5::deCacheTimeIntegrals( unsigned int ires, unsigned int islice )
{
	//Time integrals are stored under
	// ires =  resolution integral
	// islice = acceptance slice

	intExpL_stored   = storeExpL[ires][islice];
	intExpH_stored   = storeExpH[ires][islice];
	intExpSin_stored = storeExpSin[ires][islice];
	intExpCos_stored = storeExpCos[ires][islice];
}



//....................................................
// New to prepare all of the coeefficients needed in the time dependen terms
void Bs2JpsiPhi_Signal_v5::prepareCDS()
{

	double F1 = 2.0*lambda / (1.0 + lambda*lambda);
	double F2 = (1.0 - lambda*lambda) / (1.0 + lambda*lambda);

	_SS = _sinphis * F1;
	_DD = _cosphis * F1;
	_CC = F2;

}



//===========================================================================================
// Debug printout
//===========================================================================================


void Bs2JpsiPhi_Signal_v5::DebugPrint( string message, double value )  const
{
	PDF_THREAD_LOCK

		(void) message; (void) value;
	cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_Signal_v5::DebugPrint ***************************" << endl ;
	cout << message << value << endl <<endl ;

	cout << endl ;
	cout << "   gamma " << gamma() << endl ;
	cout << "   gl    " << gamma_l() << endl ;
	cout << "   gh    " << gamma_h()  << endl;
	cout << "   AT^2    " << AT()*AT() << endl;
	cout << "   AP^2    " << AP()*AP() << endl;
	cout << "   A0^2    " << A0()*A0() << endl ;
	cout << "   AS^2    " << AS()*AS() << endl ;
	cout << "   ATOTAL  " << AS()*AS()+A0()*A0()+AP()*AP()+AT()*AT() << endl ;
	cout << "   delta_ms       " << delta_ms << endl ;
	cout << "   mistag         " << mistag() << endl ;
	cout << "   mistagP1       " << _mistagP1 << endl ;
	cout << "   mistagP0       " << _mistagP0 << endl ;
	cout << "   mistagSetPoint " << _mistagSetPoint << endl ;
	cout << "   resolution " << resolution << endl ;
	cout << " For event with:  " << endl ;
	cout << "   time      " << t << endl ;
	cout << "   ctheta_tr " << ctheta_tr << endl ;
	cout << "   ctheta_1 " << ctheta_1 << endl ;
	cout << "   phi_tr " << phi_tr << endl ;

	PDF_THREAD_UNLOCK
}


void Bs2JpsiPhi_Signal_v5::DebugPrintXsec( string message, double value )  const
{
	PDF_THREAD_LOCK

		(void) message; (void) value;
	cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_Signal_v5::DebugPrintXsec ***************************" << endl ;
	cout << message << value << endl <<endl ;
	cout << "   A0()*A0() term: " <<  A0()*A0() * timeFactorA0A0(  ) * A0A0_value << endl ;
	cout << "   AP()*AP() term: " <<AP()*AP() * timeFactorAPAP(  ) * APAP_value << endl ;
	cout << "   AT()*AT() term: " <<AT()*AT() * timeFactorATAT(  ) * ATAT_value << endl << endl ;

	cout << "   AP()*AT() term: " <<AP()*AT() * timeFactorImAPAT(  ) * ImAPAT_value << endl ;
	cout << "                 : " <<AP()*AT() <<" / "<<  timeFactorImAPAT( )  <<" / "<<  ImAPAT_value << endl ;
	cout << "   A0()*AP() term: " <<A0()*AP() * timeFactorReA0AP(  ) * ReA0AP_value << endl ;
	cout << "                 : " <<A0()*AP() <<" / "<<  timeFactorReA0AP(  ) <<" / "<<  ReA0AP_value << endl ;
	cout << "   A0()*AT() term: " <<A0()*AT() * timeFactorImA0AT(  ) * ImA0AT_value << endl << endl;
	cout << "                 : " <<A0()*AT() <<" / "<<  timeFactorImA0AT(  ) <<" / "<<  ImA0AT_value << endl << endl;

	cout << "   AS()*AS() term: " <<AS()*AS() * timeFactorASAS(  ) * ASAS_value << endl << endl ;

	cout << "   AS()*AP() term: " <<AS()*AP() * timeFactorReASAP(  ) * ReASAP_value << endl ;
	cout << "                 : " <<AS()*AP() <<" / "<<   timeFactorReASAP(  ) <<" / "<<  ReASAP_value << endl ;
	cout << "   AS()*AT() term: " <<AS()*AT() * timeFactorImASAT(  ) * ImASAT_value << endl ;
	cout << "                 : " <<AS()*AT() <<" / "<<   timeFactorImASAT(  ) <<" / "<<   ImASAT_value << endl ;
	cout << "   AS()*A0() term: " <<AS()*A0() * timeFactorReASA0(  ) * ReASA0_value<< endl ;
	cout << "                 : " <<AS()*A0() <<" / "<<   timeFactorReASA0(  ) <<" / "<<  ReASA0_value << endl << endl ;

	double PwaveTot =
		A0()*A0() * timeFactorA0A0(  ) * A0A0_value +
		AP()*AP() * timeFactorAPAP(  ) * APAP_value +
		AT()*AT() * timeFactorATAT(  ) * ATAT_value +
		AP()*AT() * timeFactorImAPAT(  ) * ImAPAT_value +
		A0()*AP() * timeFactorReA0AP(  ) * ReA0AP_value +
		A0()*AT() * timeFactorImA0AT(  ) * ImA0AT_value ;

	double SwaveAdditions =
		AS()*AS() * timeFactorASAS(  ) * ASAS_value +
		AS()*AP() * timeFactorReASAP(  ) * ReASAP_value +
		AS()*AT() * timeFactorImASAT(  ) * ImASAT_value +
		AS()*A0() * timeFactorReASA0(  ) * ReASA0_value ;

	cout << "   Pwave Only : " << PwaveTot << endl ;
	cout << "   Swave add : " <<  SwaveAdditions << endl ;

	PDF_THREAD_UNLOCK
}

void Bs2JpsiPhi_Signal_v5::DebugPrintNorm( string message, double value )  const
{
	PDF_THREAD_LOCK

		(void) message; (void) value;
	cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_Signal_v5::DebugPrintNorm ***************************" << endl ;
	cout << message << value << endl <<endl ;

	cout << endl ;
	cout <<  A0()*A0() * timeFactorA0A0Int(  )* angAccI1  << endl ;
	cout <<  AP()*AP() * timeFactorAPAPInt(  )* angAccI2 << endl ;
	cout <<  AT()*AT() * timeFactorATATInt(  )* angAccI3 << endl << endl ;
	cout <<  AP()*AT() * timeFactorImAPATInt(  ) * angAccI4<< endl ;
	cout <<  A0()*AP() * timeFactorReA0APInt(  ) * angAccI5<< endl ;
	cout <<  A0()*AT() * timeFactorImA0ATInt(  ) * angAccI6<< endl << endl;
	cout <<  AS()*AS() * timeFactorASASInt(  ) * angAccI7 << endl ;
	cout <<  AS()*AP() * timeFactorReASAPInt(  ) * angAccI8<< endl ;
	cout <<  AS()*AT() * timeFactorImASATInt(  ) * angAccI9<< endl ;
	cout <<  AS()*A0() * timeFactorReASA0Int(  ) * angAccI10<< endl ;

	PDF_THREAD_UNLOCK
}

