// $Id: Bs2Jpsifzero_Signal_v6.cpp,v 1.1 2013/07/17 Dianne Ferguson Exp $
/** @class Bs2Jpsifzero_Signal_v6 Bs2Jpsifzero_Signal_v6.cpp
 *
 *  RapidFit PDF for Bs2Jpsifzero
 *
 *  @author Dianne Ferguson dferguso@cern.ch
 *  @date 2013-07-17
 */

#include "Mathematics.h"
#include "Bs2Jpsifzero_Signal_v6.h"
#include <iostream>
#include <cmath>
#include <iomanip>

#define DEBUGFLAG true

using namespace::std;

PDF_CREATOR( Bs2Jpsifzero_Signal_v6 );


Bs2Jpsifzero_Signal_v6::Bs2Jpsifzero_Signal_v6( const Bs2Jpsifzero_Signal_v6& input ) : BasePDF( (BasePDF&) input )
        , gammaName(input.gammaName), deltaGammaName(input.deltaGammaName), deltaMName(input.deltaMName), Phi_sName(input.Phi_sName)

	, Aperp_sqName(input.Aperp_sqName)

	, lambdaName(input.lambdaName), mistagName(input.mistagName)

        , mistagP1Name(input.mistagP1Name), mistagP0Name(input.mistagP0Name), mistagSetPointName(input.mistagSetPointName), mistagDeltaP1Name(input.mistagDeltaP1Name)

	, mistagDeltaP0Name(input.mistagDeltaP0Name), mistagDeltaSetPointName(input.mistagDeltaSetPointName)

        , timeName(input.timeName)

	, tagName(input.tagName)

	, _useEventResolution(input._useEventResolution), _useTimeAcceptance(input._useTimeAcceptance)

	, _numericIntegralForce(input._numericIntegralForce), _numericIntegralTimeOnly(input._numericIntegralTimeOnly)

        , _useCosAndSin(input._useCosAndSin), _usePunziMistag(input._usePunziMistag), _usePunziSigmat(input._usePunziSigmat)

	, t(input.t)

	, tag(input.tag), _gamma(input._gamma), dgam(input.dgam)

	, Aperp_sq(input.Aperp_sq)

	, delta_ms(input.delta_ms)

	, phi_s(input.phi_s), _cosphis(input._cosphis), _sinphis(input._sinphis), _mistag(input._mistag), _mistagP1(input._mistagP1), _mistagP0(input._mistagP0)

	, _mistagSetPoint(input._mistagSetPoint)

	, tlo(input.tlo), thi(input.thi), expL_stored(input.expL_stored), stored_AT(input.stored_AT)

	, expH_stored(input.expH_stored), expSin_stored(input.expSin_stored), expCos_stored(input.expCos_stored), intExpL_stored(input.intExpL_stored)

	, intExpH_stored(input.intExpH_stored), intExpSin_stored(input.intExpSin_stored), intExpCos_stored(input.intExpCos_stored), timeAcc(NULL)

	, normalisationCacheValid(input.normalisationCacheValid)

	, timeIntegralCacheValid(input.timeIntegralCacheValid)

	, normalisationCacheUntagged(input.normalisationCacheUntagged)

	, _expLObs(input._expLObs), _expHObs(input._expHObs), _expSinObs(input._expSinObs), _expCosObs(input._expCosObs), _intexpLObs(input._intexpLObs), _intexpHObs(input._intexpHObs)

	, _intexpSinObs(input._intexpSinObs), _intexpCosObs(input._intexpCosObs), _intexpLObs_vec(input._intexpLObs_vec), _intexpHObs_vec(input._intexpHObs_vec)
	
	, _intexpSinObs_vec(input._intexpSinObs_vec), _intexpCosObs_vec(input._intexpCosObs_vec), timeBinNum(input.timeBinNum), _datapoint(NULL)

	, Csp(input.Csp), lambda(input.lambda), _CC(input._CC), _DD(input._DD), _SS(input._SS), _mistagDeltaP1(input._mistagDeltaP1)

	, _mistagDeltaP0(input._mistagDeltaP0), _mistagDeltaSetPoint(input._mistagDeltaSetPoint), stored_gammal(input.stored_gammal), stored_gammah(input.stored_gammah)

{
	if( input.timeAcc != NULL ) timeAcc = new SlicedAcceptance( *(input.timeAcc) );
    resolutionModel = new ResolutionModel( *(input.resolutionModel) ) ;
}

//......................................
//Constructor(s)
//New one with configurator
Bs2Jpsifzero_Signal_v6::Bs2Jpsifzero_Signal_v6(PDFConfigurator* configurator) : BasePDF(),
	// Physics parameters
	gammaName				( configurator->getName("gamma") )
	, deltaGammaName		( configurator->getName("deltaGamma") )
	, deltaMName			( configurator->getName("deltaM") )
	, Phi_sName				( configurator->getName("Phi_s") )
	, Aperp_sqName			( configurator->getName("Aperp_sq") )
	, CspName				( configurator->getName("Csp") )
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
	// Observables
	, timeName				( configurator->getName("time") )
	, tagName				( configurator->getName("tag") )
	// Other things
	, _useEventResolution(false)
	, _useTimeAcceptance(false)
	, _numericIntegralForce(false)
	, _numericIntegralTimeOnly(false)
	, _useCosAndSin(false)
	, _usePunziMistag(false)
	, _usePunziSigmat(false)
	//objects
	,t(), tag(),
	_gamma(), dgam(), Aperp_sq(), 
	 delta_ms(), phi_s(), _cosphis(), _sinphis(), _mistag(), _mistagP1(), _mistagP0(), _mistagSetPoint(),
	tlo(), thi(), expL_stored(), expH_stored(), expSin_stored(), expCos_stored(),
	intExpL_stored(), intExpH_stored(), intExpSin_stored(), intExpCos_stored(), timeAcc(NULL), normalisationCacheValid(false),
	timeIntegralCacheValid(), normalisationCacheUntagged()
{

	std::cout << "Constructing PDF: Bs2Jpsifzero_Signal_v6 " << endl;

	//...........................................
	// Configure  options
	_numericIntegralForce    = configurator->isTrue( "NumericIntegralForce") ;
	_numericIntegralTimeOnly = configurator->isTrue( "NumericIntegralTimeOnly" ) ;
	_useEventResolution = configurator->isTrue( "UseEventResolution" ) ;
	_useCosAndSin = configurator->isTrue( "UseCosAndSin" ) ;
	_usePunziSigmat = configurator->isTrue( "UsePunziSigmat" ) ;
	_usePunziMistag = configurator->isTrue( "UsePunziMistag" ) ;

	//...........................................
	// Configure to use time acceptance machinery
	_useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;
	if( useTimeAcceptance() ) {
		if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) ) {
			timeAcc = new SlicedAcceptance( 0., 14.0, /*0.0157*/ 0.0112) ;
			cout << "Bs2Jpsifzero_Signal_v6:: Constructing timeAcc: Upper time acceptance beta=0.0112 [0 < t < 14] " << endl ;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" ) {
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ) ) ;
			cout << "Bs2Jpsifzero_Signal_v6:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl ;
		}
	}
	else {
		timeAcc = new SlicedAcceptance( 0., 14. ) ;
		cout << "Bs2Jpsifzero_Signal_v6:: Constructing timeAcc: DEFAULT FLAT [0 < t < 14]  " << endl ;
	}
    
	resolutionModel = new ResolutionModel( configurator ) ;
	if( resolutionModel->isPerEvent()  ) this->TurnCachingOff();

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
Bs2Jpsifzero_Signal_v6::~Bs2Jpsifzero_Signal_v6()
{
	if( timeAcc != NULL ) delete timeAcc;
        if( resolutionModel != NULL ) delete resolutionModel;
}


//......................................
//Make the data point and parameter set
void Bs2Jpsifzero_Signal_v6::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( tagName );
	allObservables.push_back( mistagName );

	resolutionModel->addObservables( allObservables ) ;

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( deltaGammaName );
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
	resolutionModel->addObservables( allObservables ) ;
	resolutionModel->addParameters( parameterNames ) ;
	allParameters = ParameterSet(parameterNames);
}


//.........................................................
//Return a list of observables not to be integrated
vector<string> Bs2Jpsifzero_Signal_v6::GetDoNotIntegrateList()
{
	vector<string> list;

	if( ! _usePunziMistag) list.push_back(mistagName) ;
	list.push_back("eventResolution") ;
	return list;
}


//........................................................
//Set the physics parameters into member variables

bool Bs2Jpsifzero_Signal_v6::SetPhysicsParameters( ParameterSet* NewParameterSet )
{
	bool result = allParameters.SetPhysicsParameters(NewParameterSet);
	resolutionModel->setParameters( allParameters ) ;

	// Physics parameters.
	_gamma  = allParameters.GetPhysicsParameter( gammaName )->GetValue(); 
	dgam      = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();
	Aperp_sq = 1.0;
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

	// New: Prepare the coefficients of all of the time dependent terms (C,D,S etc)
	this->prepareCDS() ;

	stored_AT = Aperp_sq > 0. ? sqrt(Aperp_sq) : 0.;
	stored_gammal = (gamma() + ( dgam *0.5 )) > 0. ? (gamma() + ( dgam *0.5 )) : 0.;
	stored_gammah = (gamma() - ( dgam *0.5 )) > 0. ? (gamma() - ( dgam *0.5 )) : 0.;
	return result;
}

//.............................................................
//Calculate the PDF value for a given set of observables for use by numeric integral

double Bs2Jpsifzero_Signal_v6::EvaluateForNumericIntegral(DataPoint * measurement)
{
	if( _numericIntegralTimeOnly ) return this->EvaluateTimeOnly(measurement) ;
	else return this->Evaluate(measurement) ;
}


//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2Jpsifzero_Signal_v6::Evaluate(DataPoint * measurement)
{
	_datapoint = measurement;

	resolutionModel->setObservables( measurement ) ;
	// Get observables into member variables
        ATAT_value = measurement->GetPseudoObservable( ATAT_Obs );
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue();

	double returnValue  = this->diffXsec();

	//conditions to throw exception
	bool c1 = std::isnan(returnValue) ;
	bool c3 =  (t>0.) && (returnValue <= 0.)  ;
	if( DEBUGFLAG && (c1 ||  c3)  ) {
		this->DebugPrint( " Bs2Jpsifzero_Signal_v6::Evaluate() returns <=0 or nan :" , returnValue ) ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	return returnValue;
}


//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2Jpsifzero_Signal_v6::EvaluateTimeOnly(DataPoint * measurement)
{
	_datapoint = measurement;
	// Get observables into member variables
        ATAT_value = measurement->GetPseudoObservable( ATAT_Obs );
	t = measurement->GetObservable( timeName )->GetValue() - timeOffset ;
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue();

	double returnValue;
	returnValue = this->diffXsecTimeOnly();

	//conditions to throw exception
	bool c1 = std::isnan(returnValue) ;
	bool c3 =  (t>0.) && (returnValue <= 0.)  ;
	if( DEBUGFLAG && (c1 ||  c3)  ) {
		this->DebugPrint( " Bs2Jpsifzero_Signal_v6::EvaluateTimeOnly() returns <=0 or nan :" , returnValue ) ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	return returnValue ;
}


//...............................................................
//Calculate the normalisation for a given set of physics parameters and boundary

double Bs2Jpsifzero_Signal_v6::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	_datapoint = measurement;
	resolutionModel->setObservables( measurement ) ;
	if( _numericIntegralForce ) return -1. ;

	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue();
	tag = (int)measurement->GetObservable( tagName )->GetValue();
	_mistag = measurement->GetObservable( mistagName )->GetValue() ;

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
	returnValue = this->diffXsecCompositeNorm1();
	// Conditions to throw exception
	bool c1 = std::isnan(returnValue)  ;
	bool c2 = (returnValue <= 0.) ;
	if( DEBUGFLAG && (c1 || c2 ) ) {
		this->DebugPrint( " Bs2Jpsifzero_Signal_v6::Normalisation :" , returnValue) ;
		if( std::isnan(returnValue) ) throw 10 ;
		if( returnValue <= 0. ) throw 10 ;
	}
	return returnValue ;
}



//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2Jpsifzero_Signal_v6::preCalculateTimeFactors()
{
        expL_stored = resolutionModel->Exp( t, gamma_l() );
        expH_stored = resolutionModel->Exp( t, gamma_h() );
        expSin_stored = resolutionModel->ExpSin( t, gamma(), delta_ms );
        expCos_stored = resolutionModel->ExpCos( t, gamma(), delta_ms );
	return;
}


//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2Jpsifzero_Signal_v6::preCalculateTimeIntegrals()
{
        intExpL_stored = resolutionModel->ExpInt( tlo, thi, gamma_l() );
        intExpH_stored = resolutionModel->ExpInt( tlo, thi, gamma_h() );
        intExpSin_stored = resolutionModel->ExpSinInt( tlo, thi, gamma(), delta_ms );
        intExpCos_stored = resolutionModel->ExpCosInt( tlo, thi, gamma(), delta_ms );
	return;
}

//...................................
// Main Diff cross section

double Bs2Jpsifzero_Signal_v6::diffXsec()
{
	preCalculateTimeFactors();
	double xsec = AT() * AT() * timeFactorATAT();

	Observable* timeObs = _datapoint->GetObservable( timeName );
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue( timeObs, timeOffset );
	if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2Jpsifzero_Signal_v6_v1::diffXsec( ) : return value < 0 = ", xsec ) ;

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
//			break;
//	}
	return xsec;
}


//...................................
// Integral over angles only for a fixed time.

double Bs2Jpsifzero_Signal_v6::diffXsecTimeOnly()
{
	preCalculateTimeFactors() ;
	double xsec = AT()*AT() * timeFactorATAT(  );

	Observable* timeObs = _datapoint->GetObservable( timeName );
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue( timeObs, timeOffset );

	if( DEBUGFLAG && (xsec < 0) ) this->DebugPrintXsec( " Bs2Jpsifzero_Signal_v6_v1::diffXsecTimeOnly( ) : return value < 0 = ", xsec ) ;

	return xsec ;
}




//...................................
// Integral over all variables: t + angles

double Bs2Jpsifzero_Signal_v6::diffXsecNorm1()
{
	preCalculateTimeIntegrals() ;//  Replaced by new Caching mechanism , but this cant be used when event resolution is selected
	double norm = AT()*AT() * timeFactorATATInt(  );

	//if( DEBUGFLAG && ((norm < 0)||(std::isnan(timeFactorATATInt()))) ) {
	if( DEBUGFLAG && ((norm < 0)||(std::isnan(norm))) ) {
		this->DebugPrintNorm( " Bs2Jpsifzero_Signal_v6_v1::diffXsecNorm1( )  ", norm ) ;

	     cout << "XXXXXXX  AT()= " <<  AT()  << "      /    timeint=   " << timeFactorATATInt(  ) << endl ;
	
	}
	return norm ;
}



//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,

double Bs2Jpsifzero_Signal_v6::diffXsecCompositeNorm1( )
{
	double tlo_boundary = tlo ;
	double thi_boundary = thi ;
	double returnValue = 0;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		timeBinNum = islice;
		tlo = tlo_boundary > timeAcc->getSlice(islice)->tlow() ? tlo_boundary : timeAcc->getSlice(islice)->tlow() ;
		thi = thi_boundary < timeAcc->getSlice(islice)->thigh() ? thi_boundary : timeAcc->getSlice(islice)->thigh() ;
		if( thi > tlo ) returnValue+= this->diffXsecNorm1(  ) * timeAcc->getSlice(islice)->height() ;
	}

	tlo = tlo_boundary;
	thi = thi_boundary;
	if( DEBUGFLAG && (std::isnan(returnValue)) ) {
		this->DebugPrintNorm( " Bs2Jpsifzero_Signal_v6_v1::diffXsecCompositeNorm1( ) : ", returnValue ) ;
	}
	return returnValue;
}


//....................................................
// New to prepare all of the coeefficients needed in the time dependen terms
void Bs2Jpsifzero_Signal_v6::prepareCDS()
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


void Bs2Jpsifzero_Signal_v6::DebugPrint( string message, double value )  const
{
	PDF_THREAD_LOCK

		(void) message; (void) value;
	cout << "*************DEBUG OUTPUT FROM Bs2Jpsifzero_Signal_v6::DebugPrint ***************************" << endl ;
	cout << message << value << endl <<endl ;

	cout << endl ;
	cout << "   gamma " << gamma() << endl ;
	cout << "   gl    " << gamma_l() << endl ;
	cout << "   gh    " << gamma_h()  << endl;
	cout << "   AT^2    " << AT()*AT() << endl;
	cout << "   delta_ms       " << delta_ms << endl ;
	cout << "   mistag         " << mistag() << endl ;
	cout << "   mistagP1       " << _mistagP1 << endl ;
	cout << "   mistagP0       " << _mistagP0 << endl ;
	cout << "   mistagSetPoint " << _mistagSetPoint << endl ;
	cout << " For event with:  " << endl ;
	cout << "   time      " << t << endl ;
	PDF_THREAD_UNLOCK
}


void Bs2Jpsifzero_Signal_v6::DebugPrintXsec( string message, double value )  const
{
	PDF_THREAD_LOCK

		(void) message; (void) value;
	cout << "*************DEBUG OUTPUT FROM Bs2Jpsifzero_Signal_v6::DebugPrintXsec ***************************" << endl ;
	cout << message << value << endl <<endl ;
	cout << "   AT()*AT() term: " <<AT()*AT() * timeFactorATAT(  ) * ATAT_value << endl << endl ;
	PDF_THREAD_UNLOCK
}

void Bs2Jpsifzero_Signal_v6::DebugPrintNorm( string message, double value )  const
{
	PDF_THREAD_LOCK

		(void) message; (void) value;
	cout << "*************DEBUG OUTPUT FROM Bs2Jpsifzero_Signal_v6::DebugPrintNorm ***************************" << endl ;
	cout << message << value << endl <<endl ;

	cout << endl ;
	cout <<  AT()*AT() * timeFactorATATInt(  )<< endl; 

	PDF_THREAD_UNLOCK
}

