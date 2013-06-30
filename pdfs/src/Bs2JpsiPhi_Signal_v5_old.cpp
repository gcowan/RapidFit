// $Id: Bs2JpsiPhi_Signal_v5_old.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_Signal_v5_old Bs2JpsiPhi_Signal_v5.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-13
 */

#include "TMath.h"
#include <cmath>

#include "Bs2JpsiPhi_Signal_v5_old.h"
#include <iostream>
#include "math.h"
#include "TF1.h"
#include "RooMath.h"
#include "Mathematics.h"

#define DEBUGFLAG true

PDF_CREATOR( Bs2JpsiPhi_Signal_v5_old );

//......................................
//Constructor(s)
//New one with configurator
Bs2JpsiPhi_Signal_v5_old::Bs2JpsiPhi_Signal_v5_old(PDFConfigurator* configurator) : 
// Physics parameters
gammaName				( configurator->getName("gamma") )
, deltaGammaName		( configurator->getName("deltaGamma") )
, deltaMName			( configurator->getName("deltaM") )
, Phi_sName				( configurator->getName("Phi_s") )
, Azero_sqName			( configurator->getName("Azero_sq") )
, Apara_sqName			( configurator->getName("Apara_sq") )
, Aperp_sqName			( configurator->getName("Aperp_sq") )
, As_sqName				( configurator->getName("F_s") )
, CspName			( configurator->getName("Csp") )
, delta_zeroName		( configurator->getName("delta_zero") )
, delta_paraName		( configurator->getName("delta_para") )
, delta_perpName		( configurator->getName("delta_perp") )
, delta_sName			( configurator->getName("delta_s") )
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

	_useHelicityBasis = configurator->isTrue( "UseHelicityBasis" ) ;

	std::cout << "Constructing PDF: Bs2JpsiPhi_Signal_v5_old " << std::endl ;
	
	//...............................................
	// Configure to use angular acceptance machinery
	string angAccFile = configurator->getConfigurationValue( "AngularAcceptanceFile" ) ;
	_angAccIgnoreNumerator = configurator->isTrue( "AngularAcceptanceIgnoreNumerator" ) ;
	if( angAccFile == "" ) cout << "Bs2JpsiPhi_Signal_v5_old:: Using flat angular acceptance " << endl ;
	else cout << "Bs2JpsiPhi_Signal_v5_old:: Constructing angAcc using file: " << angAccFile << endl ;
	angAcc = new AngularAcceptance( angAccFile, _useHelicityBasis );
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
	if( _angAccIgnoreNumerator ) cout << "Bs2JpsiPhi_Signal_v5_old:: Ignoring angular acceptance numerator " << endl ;
	
	//...........................................
	// Configure to use time acceptance machinery 
	_useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;
	if( useTimeAcceptance() ) {
		if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) ) {
			timeAcc = new SlicedAcceptance( 0., 14.0, /*0.0157*/ 0.0112) ;
			cout << "Bs2JpsiPhi_Signal_v5_old:: Constructing timeAcc: Upper time acceptance beta=0.0112 [0 < t < 14] " << endl ;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" ) {
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ) ) ;
			cout << "Bs2JpsiPhi_Signal_v5_old:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl ;
		}
	}
	else {
		timeAcc = new SlicedAcceptance( 0., 14. ) ;
		cout << "Bs2JpsiPhi_Signal_v5_old:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl ;
	}
	
	//Make empty Cache for the time integrals. This has to be done now after the SlicedAcceptance is created
	vector<double> empty ;
	for( int islice = 0; islice < timeAcc->numberOfSlices(); ++islice ) empty.push_back(0.0) ;
	for( int ires=0; ires < 4 ; ++ires ) {			
		storeExpL.push_back( empty ) ;
		storeExpH.push_back( empty ) ;
		storeExpSin.push_back( empty ) ;
		storeExpCos.push_back( empty ) ;
	}
	
	//...........................................
	// Configure numerical integration options 
	_numericIntegralForce    = configurator->isTrue( "NumericIntegralForce") ;
	_numericIntegralTimeOnly = configurator->isTrue( "NumericIntegralTimeOnly" ) ;
	
	//...........................................
	// Configure other options 
	_useEventResolution = configurator->isTrue( "UseEventResolution" ) ;
	_useCosAndSin = configurator->isTrue( "UseCosAndSin" ) ;
	_useCosDpar = configurator->isTrue( "UseCosDpar" ) ;
	_usePunziSigmat = configurator->isTrue( "UsePunziSigmat" ) ;
	_usePunziMistag = configurator->isTrue( "UsePunziMistag" ) ;
	allowNegativeAsSq = configurator->isTrue( "AllowNegativeAsSq" ) ;

	this->TurnCachingOff();

	this->SetNumericalNormalisation( false );
	
	//........................
	// Now do some actual work
	this->MakePrototypes();	
	

	//PELC  - debug to plot the distribution of PDF values for each event 
	//histOfPdfValues = new TH1D( "HistOfPdfValue" ,  "HistOfPdfValue" , 110, -0.00001, 0.00001 ) ;
	//c0  = new TCanvas;
	//histCounter = 0;
	//~PELC
	
	
}


//........................................................
//Destructor
Bs2JpsiPhi_Signal_v5_old::~Bs2JpsiPhi_Signal_v5_old() {}


//......................................
//Make the data point and parameter set
void Bs2JpsiPhi_Signal_v5_old::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	
	if( _useHelicityBasis )
	{
		allObservables.push_back( cthetakName );
		allObservables.push_back( cthetalName );
		allObservables.push_back( phihName );
	}
	else
	{
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
	parameterNames.push_back( As_sqName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( delta_zeroName );
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
vector<string> Bs2JpsiPhi_Signal_v5_old::GetDoNotIntegrateList()
{
	vector<string> list;
	
	//list.push_back(mistagName) ;
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

bool Bs2JpsiPhi_Signal_v5_old::SetPhysicsParameters( ParameterSet * NewParameterSet )
{

	normalisationCacheValid = false;  //This is only used for the untagged events and only if not useing event resolution
	timeIntegralCacheValid = false;   //This cannot be used if event resolution is used
	
	bool result = allParameters.SetPhysicsParameters(NewParameterSet);
	
	// Physics parameters. 
	_gamma  = allParameters.GetPhysicsParameter( gammaName )->GetValue();
	dgam      = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	
	Azero_sq = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
	if( (Azero_sq < 0.) || (Azero_sq > 1.)  )
	{
		PDF_THREAD_LOCK
		cout << "Warning in Bs2JpsiPhi_Signal_v5_old::SetPhysicsParameters: Azero_sq <0 or >1 but left as is" <<  endl ;
		PDF_THREAD_UNLOCK
	}
	Aperp_sq = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
	if( (Aperp_sq < 0.) || (Aperp_sq > 1.)  )
	{
		PDF_THREAD_LOCK
		cout << "Warning in Bs2JpsiPhi_Signal_v5_old::SetPhysicsParameters: Aperp_sq <0 or >1 but left as is" <<  endl ;
		PDF_THREAD_UNLOCK
	}

	Apara_sq = (1. - Azero_sq - Aperp_sq ) ;
	if( Apara_sq < 0. )
	{
		PDF_THREAD_LOCK
		cout << "Warning in Bs2JpsiPhi_Signal_v5_old::SetPhysicsParameters: derived parameter Apara_sq <0  and so set to zero" <<  endl ;
		Apara_sq = 0.;
		PDF_THREAD_UNLOCK
	}	
	
	double fs = allParameters.GetPhysicsParameter( As_sqName )->GetValue();	
	if( (fs < 0.) || (fs >= 1.) )
	{
		PDF_THREAD_LOCK
		cout << "Warning in Bs2JpsiPhi_SignalAlt_MP_dev::SetPhysicsParameters: As_sq <0 or >=1 but left as is" <<  endl;
		PDF_THREAD_UNLOCK
	}
	As_sq = fs / (1. - fs ) ;
	Csp = allParameters.GetPhysicsParameter( CspName )->GetValue();
	
	delta_zero = allParameters.GetPhysicsParameter( delta_zeroName )->GetValue();
	delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
	delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
	delta_s	   = allParameters.GetPhysicsParameter( delta_sName )->GetValue() + delta_perp;
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
	if( ! useEventResolution() )
	{
		resolution1         = allParameters.GetPhysicsParameter( res1Name )->GetValue();
		resolution2         = allParameters.GetPhysicsParameter( res2Name )->GetValue();
		resolution3         = allParameters.GetPhysicsParameter( res3Name )->GetValue();
		resolution2Fraction = allParameters.GetPhysicsParameter( res2FractionName )->GetValue();
		resolution3Fraction = allParameters.GetPhysicsParameter( res3FractionName )->GetValue();
	}
	timeOffset          = allParameters.GetPhysicsParameter( timeOffsetName )->GetValue();
	
	
	// New: Prepare the coefficients of all of the time dependent terms (C,D,S etc)
	this->prepareCDS() ;
	
	return result;
}


//.............................................................
//Calculate the PDF value for a given set of observables for use by numeric integral

double Bs2JpsiPhi_Signal_v5_old::EvaluateForNumericIntegral(DataPoint * measurement) 
{
	if( _numericIntegralTimeOnly ) return this->EvaluateTimeOnly(measurement) ;
	else return this->Evaluate(measurement) ;
}


//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_Signal_v5_old::Evaluate(DataPoint * measurement)
{

	double angAcceptanceFactor = 0 ;
	
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;

	if( _useHelicityBasis ) {
		ctheta_k   = measurement->GetObservable( cthetakName )->GetValue();
		phi_h      = TMath::Pi() + measurement->GetObservable( phihName )->GetValue();  // Pi offset is difference between angle calculator and "Our Paper"
		ctheta_l   = measurement->GetObservable( cthetalName )->GetValue();
		angAcceptanceFactor = angAcc->getValue( ctheta_l, ctheta_k, phi_h );
	}
	else {
		ctheta_tr = measurement->GetObservable( cosThetaName )->GetValue();
		phi_tr      = measurement->GetObservable( phiName )->GetValue();
		ctheta_1   = measurement->GetObservable( cosPsiName )->GetValue();
		angAcceptanceFactor = angAcc->getValue( ctheta_1, ctheta_tr, phi_tr );
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
	bool c1 = std::isnan(returnValue) ;
	bool c2 = (resolutionScale> 0.) && (returnValue <= 0.) ;
	bool c3 = (resolutionScale<=0.) && (t>0.) && (returnValue <= 0.)  ;	
	if( DEBUGFLAG && (c1 || c2 || c3)  ) {
		this->DebugPrint( " Bs2JpsiPhi_Signal_v5_old::Evaluate() returns <=0 or nan :" , returnValue ) ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	
			
	if( _angAccIgnoreNumerator ) return returnValue ;
	else return returnValue  * angAcceptanceFactor ;	
}


//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_Signal_v5_old::EvaluateTimeOnly(DataPoint * measurement)
{
	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue();
	if( useEventResolution() ) eventResolution = measurement->GetObservable( eventResolutionName )->GetValue();

	double returnValue ;
	

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
	bool c1 = std::isnan(returnValue) ;
	bool c2 = (resolutionScale> 0.) && (returnValue <= 0.) ;
	bool c3 = (resolutionScale<=0.) && (t>0.) && (returnValue <= 0.)  ;	
	if( DEBUGFLAG && (c1 || c2 || c3)  ) {
		this->DebugPrint( " Bs2JpsiPhi_Signal_v5_old::EvaluateTimeOnly() returns <=0 or nan :" , returnValue ) ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
		
	
	return returnValue ;
	
}


//...............................................................
//Calculate the normalisation for a given set of physics parameters and boundary

double Bs2JpsiPhi_Signal_v5_old::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{

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
	if( useEventResolution() )  {
		if( resolutionScale <=0. ) resolution = 0. ;
		else resolution = eventResolution * resolutionScale ;
		returnValue = this->diffXsecCompositeNorm1( 0 );
	}
		
	
	//We are not going to use event by event resolution so we can use all the caching machinery
	else {
		
		//First job for any new set of parameters is to Cache the time integrals
		if( ! timeIntegralCacheValid ) {
			CacheTimeIntegrals() ;
			timeIntegralCacheValid = true ;
			// if( ! useEventResolution() ) timeIntegralCacheValid = true ;
		}
	
	
		//If this is an untagged event and the result has been cached, then it can be used
		// Otherwise must calculate the normalisation	
		if( (tag==0) && normalisationCacheValid ) {
			returnValue = normalisationCacheUntagged ;
		}
	
		
		//So we need to calculate the normalisation 
		else {
			if( resolutionScale <= 0. ) {
				resolution = 0. ;
				returnValue = this->diffXsecCompositeNorm1( 0 );
			}
			else {						
				double val1=0. , val2=0., val3=0. ;		
				double resolution1Fraction = 1. - resolution2Fraction - resolution3Fraction ;
				if(resolution1Fraction > 0 ) {
					resolution = resolution1 * resolutionScale ;
					val1 = this->diffXsecCompositeNorm1( 1 );
				}
				if(resolution2Fraction > 0 ) {
					resolution = resolution2 * resolutionScale ;
					val2 = this->diffXsecCompositeNorm1( 2 );
				}
				if(resolution3Fraction > 0 ) {
					resolution = resolution3 * resolutionScale ;
					val3 = this->diffXsecCompositeNorm1( 3 );
				}
				returnValue = resolution1Fraction*val1 + resolution2Fraction*val2 + resolution3Fraction*val3 ;	
			}
		}
	
		// If this is an untagged event then the normaisation is invariant and so can be cached
		//if( (tag==0) && !normalisationCacheValid && !useEventResolution() )  {
		if( (tag==0) && !normalisationCacheValid )  {
			normalisationCacheUntagged = returnValue ;
			normalisationCacheValid = true ;
		}
		
	}
	
	// Conditions to throw exception
	bool c1 = std::isnan(returnValue)  ;
	bool c2 = (returnValue <= 0.) ;	
	if( DEBUGFLAG && (c1 || c2 ) ) {
		this->DebugPrint( " Bs2JpsiPhi_Signal_v5_old::Normalisation() returns <=0 or nan :" , returnValue ) ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
		
	
	return returnValue ;
}



//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_Signal_v5_old::preCalculateTimeFactors( ) const
{
	expL_stored = Mathematics::Exp( t, gamma_l(), resolution ) ;
	expH_stored = Mathematics::Exp( t, gamma_h(), resolution ) ;
	expSin_stored = Mathematics::ExpSin( t, gamma(), delta_ms, resolution ) ;
	expCos_stored = Mathematics::ExpCos( t, gamma(), delta_ms, resolution ) ;
	return ;
}


//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_Signal_v5_old::preCalculateTimeIntegrals( ) const
{
	intExpL_stored = Mathematics::ExpInt( tlo, thi, gamma_l(), resolution )  ;
	intExpH_stored = Mathematics::ExpInt( tlo, thi, gamma_h(), resolution )  ;
	intExpSin_stored = Mathematics::ExpSinInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
	intExpCos_stored = Mathematics::ExpCosInt( tlo, thi, gamma(), delta_ms, resolution ) ; 
	return ;
}

vector<string> Bs2JpsiPhi_Signal_v5_old::PDFComponents()
{
	vector<string> component_list;
	component_list.push_back( "CP-Even" );
	component_list.push_back( "CP-Odd" );
	component_list.push_back( "As" );
	component_list.push_back( "0" );
	return component_list;
}

double Bs2JpsiPhi_Signal_v5_old::EvaluateComponent( DataPoint* input, ComponentRef* Component )
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

double Bs2JpsiPhi_Signal_v5_old::diffXsec(  )  const
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

			if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue(t);
			if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2JpsiPhi_Signal_v5_old::diffXsec( ) : return value < 0 = ", xsec ) ;

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

double Bs2JpsiPhi_Signal_v5_old::diffXsecTimeOnly(  ) const
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
	
	ASint()*AP() * timeFactorReASAP(  ) * angAccI8 +
	ASint()*AT() * timeFactorImASAT(  ) * angAccI9 +
	ASint()*A0() * timeFactorReASA0(  ) * angAccI10 ;
	
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue(t);
	
	if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2JpsiPhi_Signal_v5_old::diffXsecTimeOnly( ) : return value < 0 = ", xsec ) ;
	
	return xsec ;
}




//...................................
// Integral over all variables: t + angles

double Bs2JpsiPhi_Signal_v5_old::diffXsecNorm1(  ) const
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
	
	ASint()*AP() * timeFactorReASAPInt(  ) * angAccI8 +  
	ASint()*AT() * timeFactorImASATInt(  ) * angAccI9 +  
	ASint()*A0() * timeFactorReASA0Int(  ) * angAccI10 ;  
	
	if( DEBUGFLAG && (norm < 0) ) this->DebugPrintNorm( " Bs2JpsiPhi_Signal_v5_old::diffXsecNorm1( ) : return value < 0 = ", norm ) ;
	
	return norm ;
}



//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,

double Bs2JpsiPhi_Signal_v5_old::diffXsecCompositeNorm1( int resolutionIndex )  
{   
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;
	double returnValue = 0;
	
	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		//De cache the time integrals  (unles using event Resolution
		if( ! useEventResolution() ) this->deCacheTimeIntegrals( (unsigned)resolutionIndex, islice ) ;
			
		tlo = tlo_boundary > timeAcc->getSlice((int)islice)->tlow() ? tlo_boundary : timeAcc->getSlice((int)islice)->tlow() ;
		thi = thi_boundary < timeAcc->getSlice((int)islice)->thigh() ? thi_boundary : timeAcc->getSlice((int)islice)->thigh() ;			
		if( thi > tlo ) returnValue+= this->diffXsecNorm1(  ) * timeAcc->getSlice((int)islice)->height() ;
	}		
	
	tlo = tlo_boundary;
	thi = thi_boundary ;
	return returnValue ;
}



//.......................................................
// New speed up method to Cache time integrals
void Bs2JpsiPhi_Signal_v5_old::CacheAmplitudesAndAngles() {
	
	CachedA1 = A0()*A0() * angleFactorA0A0( ) ;
	CachedA2 = AP()*AP() * angleFactorAPAP( ) ;
	CachedA3 = AT()*AT() * angleFactorATAT( ) ;
	
	CachedA4 = AP()*AT() * angleFactorImAPAT( ) ;
	CachedA5 = A0()*AP() * angleFactorReA0AP( ) ;
	CachedA6 = A0()*AT() * angleFactorImA0AT( ) ;
	
	CachedA7 = AS()*AS() * angleFactorASAS( ) ;
	
	CachedA8 = ASint()*AP() * angleFactorReASAP( ) ;
	CachedA9 = ASint()*AT() * angleFactorImASAT( ) ;
	CachedA10= ASint()*A0() * angleFactorReASA0( ) ;	
	
}

//.......................................................
// New speed up method to Cache time integrals
void Bs2JpsiPhi_Signal_v5_old::CacheTimeIntegrals() {
	
	// This need to know  (and be modified for)
	//  --> Number of resolutions
	//  --> Number fo slices
	//  --->  tlo and thi for each slice
	
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;

	if( useEventResolution() ) {
	    int ires = 0 ;
		resolution = eventResolution * resolutionScale ;
		for( unsigned int islice = 0; islice < (unsigned)timeAcc->numberOfSlices(); ++islice ) {
			tlo = tlo_boundary > timeAcc->getSlice((int)islice)->tlow() ? tlo_boundary : timeAcc->getSlice((int)islice)->tlow() ;
			thi = thi_boundary < timeAcc->getSlice((int)islice)->thigh() ? thi_boundary : timeAcc->getSlice((int)islice)->thigh() ;
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
				tlo = tlo_boundary > timeAcc->getSlice((int)islice)->tlow() ? tlo_boundary : timeAcc->getSlice((int)islice)->tlow() ;
				thi = thi_boundary < timeAcc->getSlice((int)islice)->thigh() ? thi_boundary : timeAcc->getSlice((int)islice)->thigh() ;
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
void Bs2JpsiPhi_Signal_v5_old::deCacheTimeIntegrals( unsigned int ires, unsigned int islice ) {
	
	//Time integrals are stored under 
	// ires =  resolution integral
	// islice = acceptance slice
	
	intExpL_stored   = storeExpL[ires][islice]  ;
	intExpH_stored   = storeExpH[ires][islice]  ;
	intExpSin_stored = storeExpSin[ires][islice]  ;
	intExpCos_stored = storeExpCos[ires][islice]  ;
	
	//cout << " <<<<< de-caching time integrals / " << intExpL_stored << "  /  "<< intExpH_stored << "  /  "<< intExpSin_stored << "  /  "<< intExpCos_stored << "  /  " << endl ;
	
}		



//....................................................
// New to prepare all of the coeefficients needed in the time dependen terms
void Bs2JpsiPhi_Signal_v5_old::prepareCDS()
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


void Bs2JpsiPhi_Signal_v5_old::DebugPrint( string message, double value )  const
{
	PDF_THREAD_LOCK
	(void) message; (void) value;
	cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_Signal_v5_old::DebugPrint ***************************" << endl ;
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


void Bs2JpsiPhi_Signal_v5_old::DebugPrintXsec( string message, double value )  const
{
	PDF_THREAD_LOCK
	(void) message; (void) value;
	cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_Signal_v5_old::DebugPrintXsec ***************************" << endl ;
	cout << message << value << endl <<endl ;
	cout << "   A0()*A0() term: " <<  A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) << endl ;
	cout << "   AP()*AP() term: " <<AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) << endl ;
	cout << "   AT()*AT() term: " <<AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) << endl << endl ;
	
	cout << "   AP()*AT() term: " <<AP()*AT() * timeFactorImAPAT(  ) * angleFactorImAPAT( )<< endl ;
	cout << "                 : " <<AP()*AT() <<" / "<<  timeFactorImAPAT( )  <<" / "<<  angleFactorImAPAT( )<< endl ;
	cout << "   A0()*AP() term: " <<A0()*AP() * timeFactorReA0AP(  ) * angleFactorReA0AP( )<< endl ;
	cout << "                 : " <<A0()*AP() <<" / "<<  timeFactorReA0AP(  ) <<" / "<<  angleFactorReA0AP( )<< endl ;
	cout << "   A0()*AT() term: " <<A0()*AT() * timeFactorImA0AT(  ) * angleFactorImA0AT( )<< endl << endl;
	cout << "                 : " <<A0()*AT() <<" / "<<  timeFactorImA0AT(  ) <<" / "<<  angleFactorImA0AT( )<< endl << endl;
	
	cout << "   AS()*AS() term: " <<AS()*AS() * timeFactorASAS(  ) * angleFactorASAS( ) << endl << endl ;
	
	cout << "   AS()*AP() term: " <<ASint()*AP() * timeFactorReASAP(  ) * angleFactorReASAP( )<< endl ;
	cout << "                 : " <<ASint()*AP() <<" / "<<   timeFactorReASAP(  ) <<" / "<<  angleFactorReASAP( )<< endl ;
	cout << "   AS()*AT() term: " <<ASint()*AT() * timeFactorImASAT(  ) * angleFactorImASAT( )<< endl ;
	cout << "                 : " <<ASint()*AT() <<" / "<<   timeFactorImASAT(  ) <<" / "<<   angleFactorImASAT( )<< endl ;
	cout << "   AS()*A0() term: " <<ASint()*A0() * timeFactorReASA0(  ) * angleFactorReASA0( )<< endl ;
	cout << "                 : " <<ASint()*A0() <<" / "<<   timeFactorReASA0(  ) <<" / "<<  angleFactorReASA0( )<< endl << endl ;
	
	double PwaveTot =
	A0()*A0() * timeFactorA0A0(  ) * angleFactorA0A0( ) +
	AP()*AP() * timeFactorAPAP(  ) * angleFactorAPAP( ) +
	AT()*AT() * timeFactorATAT(  ) * angleFactorATAT( ) +		
	AP()*AT() * timeFactorImAPAT(  ) * angleFactorImAPAT( ) +
	A0()*AP() * timeFactorReA0AP(  ) * angleFactorReA0AP( ) +
	A0()*AT() * timeFactorImA0AT(  ) * angleFactorImA0AT( ) ;
	
	double SwaveAdditions =
	AS()*AS() * timeFactorASAS(  ) * angleFactorASAS( ) +
	ASint()*AP() * timeFactorReASAP(  ) * angleFactorReASAP( ) +
	ASint()*AT() * timeFactorImASAT(  ) * angleFactorImASAT( ) +
	ASint()*A0() * timeFactorReASA0(  ) * angleFactorReASA0( ) ;
	
	cout << "   Pwave Only : " << PwaveTot << endl ;
	cout << "   Swave add : " <<  SwaveAdditions << endl ;
	PDF_THREAD_UNLOCK
}

void Bs2JpsiPhi_Signal_v5_old::DebugPrintNorm( string message, double value )  const
{
	PDF_THREAD_LOCK
	(void) message; (void) value;
	cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_Signal_v5_old::DebugPrintNorm ***************************" << endl ;
	cout << message << value << endl <<endl ;
	
	cout << endl ;
	cout <<  A0()*A0() * timeFactorA0A0Int(  )* angAccI1  << endl ;
	cout <<  AP()*AP() * timeFactorAPAPInt(  )* angAccI2 << endl ;
	cout <<  AT()*AT() * timeFactorATATInt(  )* angAccI3 << endl << endl ;
	cout <<  AP()*AT() * timeFactorImAPATInt(  ) * angAccI4<< endl ;
	cout <<  A0()*AP() * timeFactorReA0APInt(  ) * angAccI5<< endl ;
	cout <<  A0()*AT() * timeFactorImA0ATInt(  ) * angAccI6<< endl << endl;
	cout <<  AS()*AS() * timeFactorASASInt(  ) * angAccI7 << endl ;
	cout <<  ASint()*AP() * timeFactorReASAPInt(  ) * angAccI8<< endl ;
	cout <<  ASint()*AT() * timeFactorImASATInt(  ) * angAccI9<< endl ;
	cout <<  ASint()*A0() * timeFactorReASA0Int(  ) * angAccI10<< endl ;
	PDF_THREAD_UNLOCK
}

