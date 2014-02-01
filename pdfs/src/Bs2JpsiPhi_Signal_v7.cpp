// $Id: Bs2JpsiPhi_Signal_v7.cpp,v 1.1 2009/12/06 Pete Clarke Exp $
/** @class Bs2JpsiPhi_Signal_v7 Bs2JpsiPhi_Signal_v7.cpp
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Peter Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-13
 */

#include "ClassLookUp.h"
#include "Mathematics.h"
#include "PerEventResModel.h"
#include "DoubleResolutionModel.h"
#include "FixedResolutionModel.h"
#include "IResolutionModel.h"
#include "DoubleFixedResModel.h"
#include "TripleFixedResModel.h"
#include "Phis2012ResolutionModel.h"
#include "Bs2JpsiPhi_Angluar_Terms.h"
#include "Bs2JpsiPhi_Signal_v7.h"
#include "SimpleMistagCalib.h"
#include "CombinedMistagCalib.h"
#include "MistagCalib3fb.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
// #include "TF1.h"

using namespace::std;

PDF_CREATOR( Bs2JpsiPhi_Signal_v7 );

//......................................
//Constructor(s)
//New one with configurator
Bs2JpsiPhi_Signal_v7::Bs2JpsiPhi_Signal_v7(PDFConfigurator* configurator) : BasePDF(), 
	// Physics parameters
	gammaName			( configurator->getName("gamma") ), 
	deltaGammaName		( configurator->getName("deltaGamma") ), 
	deltaMName			( configurator->getName("deltaM") ), 
	Azero_sqName			( configurator->getName("Azero_sq") ), 
	Apara_sqName			( configurator->getName("Apara_sq") ), 
	Aperp_sqName			( configurator->getName("Aperp_sq") ), 
	delta_zeroName		( configurator->getName("delta_zero") ), 
	delta_paraName		( configurator->getName("delta_para") ), 
	delta_perpName		( configurator->getName("delta_perp") ), 
	As_sqName			( configurator->getName("F_s") ), 
	delta_sName			( configurator->getName("delta_s") ), 
	CspName			( configurator->getName("Csp") ), 
	cosdparName			( configurator->getName("cosdpar") ), //PELC-COSDPAR Special for fitting cosdpar separately
	phis_Name			( configurator->getName("Phi_s") ), 
	phis_zeroName			( configurator->getName("Phis_zero") ), 
	phis_paraName			( configurator->getName("Phis_para") ), 
	phis_perpName			( configurator->getName("Phis_perp") ), 
	phis_SName			( configurator->getName("Phis_S") ), 
	lambdaName			( configurator->getName("lambda") ), 
	lambda_zeroName		( configurator->getName("lambda_zero") ), 
	lambda_paraName		( configurator->getName("lambda_para") ), 
	lambda_perpName		( configurator->getName("lambda_perp") ), 
	lambda_SName			( configurator->getName("lambda_S") ), 
	// Observables
	timeName			( configurator->getName("time") ), 
	cosThetaName			( configurator->getName("cosTheta") ), 
	cosPsiName			( configurator->getName("cosPsi") ), 
	phiName			( configurator->getName("phi") ), 
	cthetakName 			( configurator->getName("helcosthetaK") ), 
	cthetalName			( configurator->getName("helcosthetaL") ), 
	phihName			( configurator->getName("helphi") ), 
	BetaName			( configurator->getName("beta") ), 
	// Other things
	_useEventResolution(false), 
	_useTimeAcceptance(false), 
	_useHelicityBasis(false), 
	_numericIntegralForce(false), 
	_numericIntegralTimeOnly(false), 
	_useCosAndSin(false), 
	_useCosDpar(false), 
	_usePunziMistag(false), 
	_usePunziSigmat(false), 
	allowNegativeAsSq(false), 
	_usePlotComponents(false), 
	_usePlotAllComponents(false), 
	DebugFlag_v7(true), 
	_offsetToGammaForBetaFactor(), 
	//objects
	t(), ctheta_tr(), phi_tr(), ctheta_1(), ctheta_k(), phi_h(), ctheta_l(),
	_gamma(), dgam(), Aperp_sq(), Apara_sq(), Azero_sq(), As_sq(), delta_para(),
	delta_perp(), delta_zero(), delta_s(), delta1(), delta2(), delta_ms(), phi_s(), _cosphis(), _sinphis(), 
	angAccI1(), angAccI2(), angAccI3(), angAccI4(), angAccI5(), angAccI6(), angAccI7(), angAccI8(), angAccI9(), angAccI10(),
	tlo(), thi(), expL_stored(), expH_stored(), expSin_stored(), expCos_stored(),
	intExpL_stored(), intExpH_stored(), intExpSin_stored(), intExpCos_stored(), timeAcc(NULL),
	CachedA1(), CachedA2(), CachedA3(), CachedA4(), CachedA5(), CachedA6(), CachedA7(), CachedA8(), CachedA9(), CachedA10(),
	_fitDirectlyForApara(false), performingComponentProjection(false), _useDoubleTres(false), _useTripleTres(false), _useNewPhisres(false), resolutionModel(NULL),
	_useBetaParameter(false), _useMultiplePhis(false), RequireInterference(true)
{
	componentIndex = 0;

	bool isCopy = configurator->hasConfigurationValue( "RAPIDFIT_SAYS_THIS_IS_A_COPY", "True" );

	if( !isCopy ) std::cout << "Constructing PDF: Bs2JpsiPhi_Signal_v7 " << endl;

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
	_usePlotComponents = configurator->isTrue( "PlotComponents" ) ;
	_usePlotAllComponents = configurator->isTrue( "PlotAllComponents" ) ; 
	_fitDirectlyForApara = configurator->isTrue( "FitDirectlyForApara" );
	_useNewMistagModel = configurator->isTrue( "useNewMistagModel" );
	DebugFlag_v7 = !configurator->hasConfigurationValue( "DEBUG", "False" );

	_useBetaParameter = configurator->isTrue( "floatBeta" );
	_useMultiplePhis = configurator->isTrue( "useMultiplePhis" );

	if( !_useBetaParameter )
	{
		string offsetToGammaForBetaFactor = configurator->getConfigurationValue( "OffsetToGammaForBetaFactor") ;
		if( offsetToGammaForBetaFactor == "" )
		{
			_offsetToGammaForBetaFactor = 0.0;
		}
		else
		{
			_offsetToGammaForBetaFactor = strtod( offsetToGammaForBetaFactor.c_str(), NULL ) ;
			if( !isCopy ) cout << "Bs2JpsiPhi_Signal_v7:: Adding OffsetToGammaForBetaFactor = " << _offsetToGammaForBetaFactor << endl ;
		}
	}

	//...............................................
	// Configure to use angular acceptance machinery
	string angAccFile = configurator->getConfigurationValue( "AngularAcceptanceFile" ) ;
	_angAccIgnoreNumerator = configurator->isTrue( "AngularAcceptanceIgnoreNumerator" ) ;

	if( !isCopy )
	{
		if( angAccFile == "" ) cout << "Bs2JpsiPhi_Signal_v7:: Using flat angular acceptance " << endl ;
		else cout << "Bs2JpsiPhi_Signal_v7:: Constructing angAcc using file: " << angAccFile << endl ;
		angAcc = new AngularAcceptance( angAccFile, _angAccIgnoreNumerator, _useHelicityBasis ) ;
		angAccI1 = angAcc->af1() ;      cout << "  af1 = " << setprecision(6) << setw(15) << angAccI1 << setw(10) << " ";
		angAccI2 = angAcc->af2() ;	cout << "  af2 = " << setprecision(6) << setw(15) << angAccI2 << setw(10) << " ";
		angAccI3 = angAcc->af3() ;	cout << "  af3 = " << setprecision(6) << setw(15) << angAccI3 << endl ;
		angAccI4 = angAcc->af4() ;	cout << "  af4 = " << setprecision(6) << setw(15) << angAccI4 << setw(10) << " ";
		angAccI5 = angAcc->af5() ;	cout << "  af5 = " << setprecision(6) << setw(15) << angAccI5 << setw(10) << " ";
		angAccI6 = angAcc->af6() ;	cout << "  af6 = " << setprecision(6) << setw(15) << angAccI6 << endl ;
		angAccI7 = angAcc->af7() ;	cout << "  af7 = " << setprecision(6) << setw(15) << angAccI7 << setw(10) << " ";
		angAccI8 = angAcc->af8() ;	cout << "  af8 = " << setprecision(6) << setw(15) << angAccI8 << setw(10) << " ";
		angAccI9 = angAcc->af9() ;	cout << "  af9 = " << setprecision(6) << setw(15) << angAccI9 << setw(10) << " ";
		angAccI10 = angAcc->af10();	cout << "  af10 = " << setprecision(6) << setw(15) << angAccI10 << endl ;
		if( _angAccIgnoreNumerator ) cout << "Bs2JpsiPhi_Signal_v7:: Ignoring angular acceptance numerator " << endl ;
	}
	else
	{
		angAcc = new AngularAcceptance( angAccFile, _angAccIgnoreNumerator, _useHelicityBasis, isCopy ) ;
		angAccI1 = angAcc->af1()  ;
		angAccI2 = angAcc->af2()  ;
		angAccI3 = angAcc->af3()  ;
		angAccI4 = angAcc->af4()  ;
		angAccI5 = angAcc->af5()  ;
		angAccI6 = angAcc->af6()  ;
		angAccI7 = angAcc->af7()  ;
		angAccI8 = angAcc->af8()  ;
		angAccI9 = angAcc->af9()  ;
		angAccI10 = angAcc->af10();
	}
	//...........................................
	// Configure to use time acceptance machinery
	_useTimeAcceptance = configurator->isTrue( "UseTimeAcceptance" ) ;
	if( useTimeAcceptance() )
	{
		if( configurator->hasConfigurationValue( "TimeAcceptanceType", "Upper" ) )
		{
			timeAcc = new SlicedAcceptance( 0., 14.0, 0.00826, isCopy) ;
			//timeAcc = new SlicedAcceptance( 0., 14.0, 0.00, isCopy) ;
			if( !isCopy ) cout << "Bs2JpsiPhi_Signal_v7:: Constructing timeAcc: Upper time acceptance beta=0.00826 [0 < t < 14] " << endl ;
		}
		else if( configurator->getConfigurationValue( "TimeAcceptanceFile" ) != "" )
		{
			timeAcc = new SlicedAcceptance( "File" , configurator->getConfigurationValue( "TimeAcceptanceFile" ), isCopy ) ;
			if( !isCopy ) cout << "Bs2JpsiPhi_Signal_v7:: Constructing timeAcc: using file: " << configurator->getConfigurationValue( "TimeAcceptanceFile" ) << endl ;
		}
	}

	if( timeAcc == NULL )
	{
		timeAcc = new SlicedAcceptance( 0., 20., isCopy ) ;
		if( !isCopy ) cout << "Bs2JpsiPhi_Signal_v7:: Constructing timeAcc: DEFAULT FLAT [0 < t < 20]  " << endl ;
	}

	this->SetNumericalNormalisation( false );

	//..........................................
	// Choose resolution model according to flags
	// For now hard coded.

	//if( _useEventResolution )
		// get model name
		// use class lookup
	//{
		string resolutionModelName = configurator->GetResolutionModel();
		resolutionModel = ClassLookUp::LookUpResName( resolutionModelName, configurator, isCopy );

		if( _useEventResolution ) this->TurnCachingOff();
	//}
	//else
	//{
	//	resolutionModel = ClassLookUp::LookUpResName( "", configurator, isCopy );
	//}

	//if( _useNewMistagModel ) _mistagCalibModel = new CombinedMistagCalib( configurator );
	if( _useNewMistagModel ) _mistagCalibModel = new MistagCalib3fb( configurator );
	else _mistagCalibModel = new SimpleMistagCalib( configurator );

	//........................
	// Now do some actual work
	this->MakePrototypes();

	//PELC  - debug to plot the distribution of PDF values for each event
	//histOfPdfValues = new TH1D( "HistOfPdfValue" ,  "HistOfPdfValue" , 110, -0.00001, 0.00001 ) ;
	//c0  = new TCanvas;
	//histCounter = 0;
	//~PELC

	this->SetCopyConstructorSafe( false );
}

//........................................................
//Destructor
Bs2JpsiPhi_Signal_v7::~Bs2JpsiPhi_Signal_v7()
{
	if( timeAcc != NULL ) delete timeAcc;
	if( angAcc != NULL ) delete angAcc;
	if( resolutionModel != NULL ) delete resolutionModel;
	if( _mistagCalibModel != NULL ) delete _mistagCalibModel;
}


//......................................
//Make the data point and parameter set
void Bs2JpsiPhi_Signal_v7::MakePrototypes()
{
	//............................
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
	//if(useEventResolution()) allObservables.push_back( eventResolutionName );    
	resolutionModel->addObservables( allObservables );
	_mistagCalibModel->addObservables( allObservables );

	//...............................
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

	if( _useCosAndSin )
	{
		parameterNames.push_back( cosphisName );
		parameterNames.push_back( sinphisName );
		parameterNames.push_back( lambdaName );
	}
	else if( _useMultiplePhis )
	{
		parameterNames.push_back( phis_zeroName );
		parameterNames.push_back( phis_paraName );
		parameterNames.push_back( phis_perpName );
		parameterNames.push_back( phis_SName );
		parameterNames.push_back( lambda_zeroName );
		parameterNames.push_back( lambda_paraName );
		parameterNames.push_back( lambda_perpName );
		parameterNames.push_back( lambda_SName );
	}
	else
	{
		parameterNames.push_back( phis_Name );
		parameterNames.push_back( lambdaName );
	}

	if( _fitDirectlyForApara )
	{
		parameterNames.push_back( Apara_sqName );
	}

	//Now add additional paramaters for the resolution model
	resolutionModel->addParameters( parameterNames );
	_mistagCalibModel->addParameters( parameterNames );

	if( _useBetaParameter ) parameterNames.push_back( BetaName );

	allParameters = ParameterSet(parameterNames);
}


//.........................................................
//Return a list of observables not to be integrated
vector<string> Bs2JpsiPhi_Signal_v7::GetDoNotIntegrateList()
{
	vector<string> list;

	//***THIS NEEDS FIXING*** PROBLEM IS ONLY THE RESOULTION MODEL KNOWS THIS NOW - SO HOW DOES ONE ADD TO THE D.I.L CLEANLY ??
	//list.push_back("eventResolution") ; 
	resolutionModel->addObservables( list );

	if( !_usePunziMistag ) _mistagCalibModel->addObservables( list );

	if( _numericIntegralTimeOnly )
	{
		if( _useHelicityBasis )
		{
			list.push_back( cthetakName );
			list.push_back( cthetalName ) ;
			list.push_back( phihName ) ;
		}
		else
		{
			list.push_back( cosThetaName );
			list.push_back( cosPsiName ) ;
			list.push_back( phiName ) ;
		}
	}
	return list;
}


//........................................................
//Set the physics parameters into member variables

bool Bs2JpsiPhi_Signal_v7::SetPhysicsParameters( ParameterSet* NewParameterSet )
{
	// normalisationCacheValid = false;  //This is only used for the untagged events and only if not useing event resolution
	//timeIntegralCacheValid = false;   //This cannot be used if event resolution is used

	bool result = allParameters.SetPhysicsParameters(NewParameterSet);

	//Let the resolution model take its specific parameters out
	resolutionModel->setParameters( allParameters );
	_mistagCalibModel->setParameters( allParameters );

	if( _useBetaParameter ) _offsetToGammaForBetaFactor = allParameters.GetPhysicsParameter( BetaName )->GetValue();

	// Physics parameters.
	_gamma  = allParameters.GetPhysicsParameter( gammaName )->GetValue() + _offsetToGammaForBetaFactor ;
	dgam      = allParameters.GetPhysicsParameter( deltaGammaName )->GetValue();

	Azero_sq = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
	if( (Azero_sq < 0.) || (Azero_sq > 1.)  )
	{
		cout << "Warning in Bs2JpsiPhi_Signal_v7::SetPhysicsParameters: Azero_sq <0 or >1 but left as is" <<  endl ;
	}
	Aperp_sq = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
	if( (Aperp_sq < 0.) || (Aperp_sq > 1.)  )
	{
		cout << "Warning in Bs2JpsiPhi_Signal_v7::SetPhysicsParameters: Aperp_sq <0 or >1 but left as is" <<  endl;
	}

	Apara_sq = 0.;
	if( _fitDirectlyForApara )
	{
		Apara_sq = allParameters.GetPhysicsParameter( Apara_sqName )->GetValue();
	}
	else
	{
		Apara_sq = (1. - Azero_sq - Aperp_sq );
	}
	if( Apara_sq < 0. )
	{
		cout << "Warning in Bs2JpsiPhi_Signal_v7::SetPhysicsParameters: derived parameter Apara_sq <0  and so set to zero" <<  endl ;
		Apara_sq = 0. ;
	}

	double fs = allParameters.GetPhysicsParameter( As_sqName )->GetValue();
	if( (fs < 0.) || (fs >= 1.) )
	{
		cout << "Warning in Bs2JpsiPhi_Signal_v7::SetPhysicsParameters: As_sq <0 or >=1 but left as is" <<  endl ;
	}
	As_sq = fs / (1. - fs );

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

	if( _useMultiplePhis )
	{
		phis_zeroVal = allParameters.GetPhysicsParameter( phis_zeroName )->GetValue();
		phis_paraVal = allParameters.GetPhysicsParameter( phis_paraName )->GetValue();
		phis_perpVal = allParameters.GetPhysicsParameter( phis_perpName )->GetValue();
		phis_SVal = allParameters.GetPhysicsParameter( phis_SName )->GetValue();
		lambda_zeroVal = allParameters.GetPhysicsParameter( lambda_zeroName )->GetValue();
		lambda_paraVal = allParameters.GetPhysicsParameter( lambda_paraName )->GetValue();
		lambda_perpVal = allParameters.GetPhysicsParameter( lambda_perpName )->GetValue();
		lambda_SVal = allParameters.GetPhysicsParameter( lambda_SName )->GetValue();
	}
	else
	{

		if(_useCosAndSin)
		{
			_cosphis = allParameters.GetPhysicsParameter( cosphisName )->GetValue();
			_sinphis = allParameters.GetPhysicsParameter( sinphisName )->GetValue();
			phi_s = 0.5*( acos( _cosphis ) + asin( _sinphis ) );
		}
		else
		{
			phi_s     = allParameters.GetPhysicsParameter( phis_Name )->GetValue();
			_cosphis = cos(phi_s);
			_sinphis = sin(phi_s);
		}

		lambda = allParameters.GetPhysicsParameter( lambdaName )->GetValue();
		if( !_useCosAndSin )
		{
			phis_zeroVal = phi_s;
			phis_paraVal = phi_s;
			phis_perpVal = phi_s;
			phis_SVal = phi_s;
		}
		lambda_zeroVal = lambda;
		lambda_paraVal = lambda;
		lambda_perpVal = lambda;
		lambda_SVal = lambda;
	}

	// New: Prepare the coefficients of all of the time dependent terms (C,D,S etc)
	this->prepareCDS( lambda, phi_s );

	stored_AT = Aperp_sq > 0. ? sqrt(Aperp_sq) : 0.;
	stored_A0 = Azero_sq > 0. ? sqrt(Azero_sq) : 0.;
	stored_AP = Apara_sq > 0. ? sqrt(Apara_sq) : 0.;
	stored_AS = As_sq > 0. ? sqrt(As_sq) : 0.;
	stored_ASint = stored_AS * Csp;

	stored_gammal = (gamma() + ( dgam *0.5 )) > 0. ? (gamma() + ( dgam *0.5 )) : 0.;
	stored_gammah = (gamma() - ( dgam *0.5 )) > 0. ? (gamma() - ( dgam *0.5 )) : 0.;

	double delta_perp_s = delta_perp - delta_s;
	double delta_para_s = delta_para - delta_s;
	double delta_zero_s = delta_zero - delta_s;
	sin_delta_perp_s = sin(delta_perp_s);
	cos_delta_perp_s = cos(delta_perp_s);
	sin_delta_zero_s = sin(delta_zero_s);
	cos_delta_zero_s = cos(delta_zero_s);
	sin_delta_para_s = sin(delta_para_s);
	cos_delta_para_s = cos(delta_para_s);

	sin_delta1 = sin(delta1);
	cos_delta1 = cos(delta1);
	sin_delta2 = sin(delta2);
	cos_delta2 = cos(delta2);
	double delta_2_1 = delta2 - delta1;
	sin_delta_2_1 = sin(delta_2_1);
	cos_delta_2_1 = cos(delta_2_1);

	return result;
}

//.............................................................
//Calculate the PDF value for a given set of observables for use by numeric integral
double Bs2JpsiPhi_Signal_v7::EvaluateForNumericIntegral(DataPoint * measurement)
{
	if( _numericIntegralTimeOnly ) return this->EvaluateTimeOnly(measurement) ;
	else return this->Evaluate(measurement) ;
}


//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_Signal_v7::Evaluate(DataPoint * measurement)
{
	_datapoint = measurement;

	//Let the resolution model pull out its specific obsrvables first.
	//This can only be the case if event resolution is used (so far)
	resolutionModel->setObservables( measurement );
	_mistagCalibModel->setObservables( measurement );

	_eventIsTagged = _mistagCalibModel->eventIsTagged();

	vector<double> angularData;

	if( !_useHelicityBasis )
	{
		angularData.push_back( measurement->GetObservable( cosThetaName )->GetValue() );
		angularData.push_back( measurement->GetObservable( cosPsiName )->GetValue() );
		angularData.push_back( measurement->GetObservable( phiName )->GetValue() );

		A0A0_value = Bs2JpsiPhi_Angular_Terms::TangleFactorA0A0( angularData );
		APAP_value = Bs2JpsiPhi_Angular_Terms::TangleFactorAPAP( angularData );
		ATAT_value = Bs2JpsiPhi_Angular_Terms::TangleFactorATAT( angularData );
		ASAS_value = Bs2JpsiPhi_Angular_Terms::TangleFactorASAS( angularData );
		ImAPAT_value = Bs2JpsiPhi_Angular_Terms::TangleFactorImAPAT( angularData );
		ReA0AP_value = Bs2JpsiPhi_Angular_Terms::TangleFactorReA0AP( angularData );
		ImA0AT_value = Bs2JpsiPhi_Angular_Terms::TangleFactorImA0AT( angularData );
		ReASAP_value = Bs2JpsiPhi_Angular_Terms::TangleFactorReASAP( angularData );
		ImASAT_value = Bs2JpsiPhi_Angular_Terms::TangleFactorImASAT( angularData );
		ReASA0_value = Bs2JpsiPhi_Angular_Terms::TangleFactorReASA0( angularData );
	}
	else
	{
		angularData.push_back( measurement->GetObservable( cthetakName )->GetValue() );
		angularData.push_back( measurement->GetObservable( cthetalName )->GetValue() );
		angularData.push_back( measurement->GetObservable( phihName )->GetValue() );

		A0A0_value = Bs2JpsiPhi_Angular_Terms::HangleFactorA0A0( angularData );
		APAP_value = Bs2JpsiPhi_Angular_Terms::HangleFactorAPAP( angularData );
		ATAT_value = Bs2JpsiPhi_Angular_Terms::HangleFactorATAT( angularData );
		ASAS_value = Bs2JpsiPhi_Angular_Terms::HangleFactorASAS( angularData );
		ImAPAT_value = Bs2JpsiPhi_Angular_Terms::HangleFactorImAPAT( angularData );
		ReA0AP_value = Bs2JpsiPhi_Angular_Terms::HangleFactorReA0AP( angularData );
		ImA0AT_value = Bs2JpsiPhi_Angular_Terms::HangleFactorImA0AT( angularData );
		ReASAP_value = Bs2JpsiPhi_Angular_Terms::HangleFactorReASAP( angularData );
		ImASAT_value = Bs2JpsiPhi_Angular_Terms::HangleFactorImASAT( angularData );
		ReASA0_value = Bs2JpsiPhi_Angular_Terms::HangleFactorReASA0( angularData );
	}

	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() ; // - timeOffset ;

	//Get anglular quantities
	double angAcceptanceFactor = 0 ;
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

	//Cache amplitues and angles terms used in cross section
	this->CacheAmplitudesAndAngles() ;

	double returnValue = this->diffXsec( );

	if( !performingComponentProjection )
	{
		if( DebugFlag_v7 )
		{
			//conditions to throw exception
			bool c1 = std::isnan(returnValue);
			bool c3 = (t>0.) && (returnValue <= 0.);
			if( c1 || c3 )
			{
				this->DebugPrint( " Bs2JpsiPhi_Signal_v7::Evaluate() returns <=0 or nan :" , returnValue ) ;
				if( std::isnan(returnValue) ) throw 10 ;
				if( returnValue <= 0. ) throw 10 ;
			}
		}
	}


	if( ! _angAccIgnoreNumerator ) returnValue*=angAcceptanceFactor;

	return returnValue;
}


//.............................................................
//Calculate the PDF value for a given set of observables

double Bs2JpsiPhi_Signal_v7::EvaluateTimeOnly(DataPoint * measurement)
{
	_datapoint = measurement;

	resolutionModel->setObservables( measurement );
	_mistagCalibModel->setObservables( measurement );

	_eventIsTagged = _mistagCalibModel->eventIsTagged();

	vector<double> angularData;

	if( !_useHelicityBasis )
	{
		angularData.push_back( measurement->GetObservable( cosThetaName )->GetValue() );
		angularData.push_back( measurement->GetObservable( cosPsiName )->GetValue() );
		angularData.push_back( measurement->GetObservable( phiName )->GetValue() );

		A0A0_value = Bs2JpsiPhi_Angular_Terms::TangleFactorA0A0( angularData );
		APAP_value = Bs2JpsiPhi_Angular_Terms::TangleFactorAPAP( angularData );
		ATAT_value = Bs2JpsiPhi_Angular_Terms::TangleFactorATAT( angularData );
		ASAS_value = Bs2JpsiPhi_Angular_Terms::TangleFactorASAS( angularData );
		ImAPAT_value = Bs2JpsiPhi_Angular_Terms::TangleFactorImAPAT( angularData );
		ReA0AP_value = Bs2JpsiPhi_Angular_Terms::TangleFactorReA0AP( angularData );
		ImA0AT_value = Bs2JpsiPhi_Angular_Terms::TangleFactorImA0AT( angularData );
		ReASAP_value = Bs2JpsiPhi_Angular_Terms::TangleFactorReASAP( angularData );
		ImASAT_value = Bs2JpsiPhi_Angular_Terms::TangleFactorImASAT( angularData );
		ReASA0_value = Bs2JpsiPhi_Angular_Terms::TangleFactorReASA0( angularData );
	}
	else
	{
		angularData.push_back( measurement->GetObservable( cthetakName )->GetValue() );
		angularData.push_back( measurement->GetObservable( cthetalName )->GetValue() );
		angularData.push_back( measurement->GetObservable( phihName )->GetValue() );

		A0A0_value = Bs2JpsiPhi_Angular_Terms::HangleFactorA0A0( angularData );
		APAP_value = Bs2JpsiPhi_Angular_Terms::HangleFactorAPAP( angularData );
		ATAT_value = Bs2JpsiPhi_Angular_Terms::HangleFactorATAT( angularData );
		ASAS_value = Bs2JpsiPhi_Angular_Terms::HangleFactorASAS( angularData );
		ImAPAT_value = Bs2JpsiPhi_Angular_Terms::HangleFactorImAPAT( angularData );
		ReA0AP_value = Bs2JpsiPhi_Angular_Terms::HangleFactorReA0AP( angularData );
		ImA0AT_value = Bs2JpsiPhi_Angular_Terms::HangleFactorImA0AT( angularData );
		ReASAP_value = Bs2JpsiPhi_Angular_Terms::HangleFactorReASAP( angularData );
		ImASAT_value = Bs2JpsiPhi_Angular_Terms::HangleFactorImASAT( angularData );
		ReASA0_value = Bs2JpsiPhi_Angular_Terms::HangleFactorReASA0( angularData );
	}

	// Get observables into member variables
	t = measurement->GetObservable( timeName )->GetValue() ; // - timeOffset ;

	double returnValue = this->diffXsecTimeOnly( );

	if( DebugFlag_v7 )
	{
		//conditions to throw exception
		bool c1 = std::isnan(returnValue) ;
		bool c3 = (t>0.) && (returnValue <= 0.)  ;
		if( (c1 || c3)  )
		{
			this->DebugPrintEvaluate( " Bs2JpsiPhi_Signal_v7::EvaluateTimeOnly() returns <=0 or nan :" , returnValue ) ;
			if( std::isnan(returnValue) ) throw 10 ;
			if( returnValue <= 0. ) throw 10 ;
		}
	}


	return returnValue ;

}


//...............................................................
//Calculate the normalisation for a given set of physics parameters and boundary

double Bs2JpsiPhi_Signal_v7::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
	_datapoint = measurement;

	resolutionModel->setObservables( measurement );
	_mistagCalibModel->setObservables( measurement );

	_eventIsTagged = _mistagCalibModel->eventIsTagged();

	if( _numericIntegralForce ) return -1.;

	// Get time boundaries into member variables
	IConstraint * timeBound = boundary->GetConstraint( timeName );
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		return 0;
	}
	else
	{
		tlo = timeBound->GetMinimum();
		thi = timeBound->GetMaximum();
	}

	double returnValue = this->diffXsecCompositeNorm1( );

	if( !performingComponentProjection )
	{
		if( DebugFlag_v7 )
		{
			// Conditions to throw exception
			bool c1 = std::isnan(returnValue);
			bool c2 = (returnValue <= 0.);
			if( c1 || c2 )
			{
				this->DebugPrintNormalisation( " Bs2JpsiPhi_Signal_v7::Normalisation() returns <=0 or nan :" , returnValue ) ;
				if( std::isnan(returnValue) ) throw 10 ;
				if( returnValue <= 0. ) throw 10 ;
			}
		}
	}

	return returnValue;
}



//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_Signal_v7::preCalculateTimeFactors()
{
	expL_stored = resolutionModel->Exp( t, gamma_l() );
	expH_stored = resolutionModel->Exp( t, gamma_h() );

	if( _eventIsTagged && RequireInterference )
	{
		expSin_stored = resolutionModel->ExpSin( t, gamma(), delta_ms );
		expCos_stored = resolutionModel->ExpCos( t, gamma(), delta_ms );
	}
	else
	{
		expSin_stored = 0.;
		expCos_stored = 0.;
	}
	return;
}


//.......................................................
// Pre calculate the time integrals : this is becaue these functions are called many times for each event due to the 10 angular terms
void Bs2JpsiPhi_Signal_v7::preCalculateTimeIntegrals()
{
	intExpL_stored = resolutionModel->ExpInt( tlo, thi, gamma_l() );
	intExpH_stored = resolutionModel->ExpInt( tlo, thi, gamma_h() );
}

void Bs2JpsiPhi_Signal_v7::preCalculateSinusoidIntegrals()
{
	if( _eventIsTagged && RequireInterference )
	{
		intExpSin_stored = resolutionModel->ExpSinInt( tlo, thi, gamma(), delta_ms );
		intExpCos_stored = resolutionModel->ExpCosInt( tlo, thi, gamma(), delta_ms );
	}
	else
	{
		intExpSin_stored = 0.;
		intExpCos_stored = 0.;
	}
	return;
}

void Bs2JpsiPhi_Signal_v7::preCalculateTimeIntegralsDebug() const
{
	cout << "Timey-Wimey stuff:" << endl;
	cout << "expL_stored: " << expL_stored << endl;
	cout << "expH_stored: " << expH_stored << endl;
	cout << "expSin_stored: " << expSin_stored << endl;
	cout << "expCos_stored: " << expCos_stored << endl;
	cout << "intExpL_stored: " << intExpL_stored << endl;
	cout << "intExpH_stored: " << intExpH_stored << endl;
	cout << "intExpSin_stored: " << intExpSin_stored << endl;
	cout << "intExpCos_stored: " << intExpCos_stored << endl;
	cout << endl;
}

vector<string> Bs2JpsiPhi_Signal_v7::PDFComponents()
{
	vector<string> this_component_list;
	if( _usePlotComponents && !_usePlotAllComponents ) {
		if( allParameters.GetPhysicsParameter(Azero_sqName)->GetValue() > 1E-10 )
		{
			this_component_list.push_back( "CP-Even" );
		}
		if( allParameters.GetPhysicsParameter(Aperp_sqName)->GetValue() > 1E-10 )
		{
			this_component_list.push_back( "CP-Odd" );
		}
		if( allParameters.GetPhysicsParameter(As_sqName)->GetValue() > 1E-10 )
		{
			this_component_list.push_back( "As" );
		}
		this_component_list.push_back( "0" );
	}
	else if( _usePlotAllComponents && !_usePlotComponents )
	{
		this_component_list.push_back( "0" );
		this_component_list.push_back( "A0A0" );
		this_component_list.push_back( "AparaApara" );
		this_component_list.push_back( "AperpAperp" );
		this_component_list.push_back( "ImAparaAperp" );
		this_component_list.push_back( "ReA0Apara" );
		this_component_list.push_back( "ImA0Aperp" );
		this_component_list.push_back( "AsAs" );
		this_component_list.push_back( "ImAsApara" );
		this_component_list.push_back( "ReAsAperp" );
		this_component_list.push_back( "ImAsA0" );
	}

	return this_component_list;
}

double Bs2JpsiPhi_Signal_v7::EvaluateComponent( DataPoint* input, ComponentRef* Component )
{
	performingComponentProjection = true;
	componentIndex = Component->getComponentNumber();

	RequireInterference = false;

	if( _usePlotComponents && !_usePlotAllComponents )
	{
		if( componentIndex == -1 )
		{
			string ComponentName = Component->getComponentName();
			if( ComponentName.compare( "CP-Even" ) == 0 )
			{
				Component->setComponentNumber( 1 );
				componentIndex = 1;
				RequireInterference = true;
			}
			else if( ComponentName.compare( "CP-Odd" ) == 0 )
			{
				Component->setComponentNumber( 2 );
				componentIndex = 2;
				RequireInterference = false;
			}
			else if( ComponentName.compare( "As" ) == 0 )
			{
				Component->setComponentNumber( 3 );
				componentIndex = 3;
				RequireInterference = false;
			}
			else
			{
				Component->setComponentNumber( 0 );
				componentIndex = 0;
				RequireInterference = true;
			}
		}
	}
	else if( _usePlotAllComponents && !_usePlotComponents )
	{
		if( componentIndex == -1 )
		{
			string ComponentName = Component->getComponentName();
			if( ComponentName.compare( "A0A0" ) == 0 )
			{
				Component->setComponentNumber( 1 );
				componentIndex = 1;
				RequireInterference = false;
			}
			else if( ComponentName.compare( "AparaApara" ) == 0 )
			{
				Component->setComponentNumber( 2 );
				componentIndex = 2;
				RequireInterference = false;
			}
			else if( ComponentName.compare( "AperpAperp" ) == 0 )
			{
				Component->setComponentNumber( 3 );
				componentIndex = 3;
				RequireInterference = false;
			}
			else if( ComponentName.compare( "ImAparaAperp" ) == 0 )
			{
				Component->setComponentNumber( 4 );
				componentIndex = 4;
				RequireInterference = true;
			}
			else if( ComponentName.compare( "ReA0Apara" ) == 0 )
			{
				Component->setComponentNumber( 6 );
				componentIndex = 5;
				RequireInterference = true;
			}
			else if( ComponentName.compare( "ImA0Aperp" ) == 0 )
			{
				Component->setComponentNumber( 6 );
				componentIndex = 6;
				RequireInterference = true;
			}
			else if( ComponentName.compare( "AsAs" ) == 0 )
			{
				Component->setComponentNumber( 7 );
				componentIndex = 7;
				RequireInterference = false;
			}
			else if( ComponentName.compare( "ImAsApara" ) == 0 )
			{
				Component->setComponentNumber( 8 );
				componentIndex = 8;
				RequireInterference = true;
			}
			else if( ComponentName.compare( "ReAsAperp" ) == 0 )
			{
				Component->setComponentNumber( 9 );
				componentIndex = 9;
				RequireInterference = true;
			}
			else if( ComponentName.compare( "ImAsA0" ) == 0 )
			{
				Component->setComponentNumber( 10 );
				componentIndex = 10;
				RequireInterference = true;
			}
			else
			{
				Component->setComponentNumber( 0 );
				componentIndex = 0;
			}
		}
	}
	else
	{
		RequireInterference = true;
		return this->Evaluate( input );
	}

	double return_value = this->Evaluate( input );
	componentIndex = 0;

	performingComponentProjection = false;
	RequireInterference = true;
	return return_value;
}

//...................................
// Main Diff cross section

double Bs2JpsiPhi_Signal_v7::diffXsec()
{
	preCalculateTimeFactors();

	double xsec=-1.;
	if( _usePlotAllComponents && !_usePlotComponents )
	{
		switch( componentIndex )
		{
			case 1:
				xsec = CachedA1 * timeFactorA0A0(  );
				break;
			case 2:
				xsec = CachedA2 * timeFactorAPAP(  );
				break;
			case 3:
				xsec = CachedA3 * timeFactorATAT(  );
				break;
			case 4:
				xsec = CachedA4 * timeFactorImAPAT(  );
				break;
			case 5:
				xsec = CachedA5 * timeFactorReA0AP(  );
				break;
			case 6:
				xsec = CachedA6 * timeFactorImA0AT(  );
				break;
			case 7:
				xsec = CachedA7 * timeFactorASAS(  );
				break;
			case 8:
				xsec = CachedA8 * timeFactorReASAP(  );
				break;
			case 9:
				xsec = CachedA9 * timeFactorImASAT(  );
				break;
			case 10:
				xsec = CachedA10 * timeFactorReASA0(  );
				break;
			default:
				xsec = CachedA1 * timeFactorA0A0(  );
				xsec+= CachedA2 * timeFactorAPAP(  );
				xsec+= CachedA3 * timeFactorATAT(  );
				xsec+= CachedA4 * timeFactorImAPAT(  );
				xsec+= CachedA5 * timeFactorReA0AP(  );
				xsec+= CachedA6 * timeFactorImA0AT(  );
				xsec+= CachedA7 * timeFactorASAS(  );
				xsec+= CachedA8 * timeFactorReASAP(  );
				xsec+= CachedA9 * timeFactorImASAT(  );
				xsec+= CachedA10 * timeFactorReASA0(  );
				break;
		}
	}
	else if( !_usePlotAllComponents && _usePlotComponents )
	{
		switch( componentIndex )
		{
			case 1:         //      CP-Even         CP-Odd=0 && S-Wave=0
				xsec = CachedA1 * timeFactorA0A0(  );
				xsec += CachedA2 * timeFactorAPAP(  );
				xsec += CachedA5 * timeFactorReA0AP(  );
				break;
			case 2:         //      CP-Odd          CP-Even=0 && S-Wave=0
				xsec = CachedA3 * timeFactorATAT(  );
				break;
			case 3:         //      S-Wave          CP-Even=0 && CP-Odd=0
				xsec = CachedA7 * timeFactorASAS(  );
				break;
			default:        //      Everything
				xsec =

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

				//PELC - This turned out to be an important debugging tool
				//switch it on to see the values of PDF being returend.  If ANY go negative, it means there is a sign wrong in one or more of the terms
				//You need to enable in the .h file as well
				//histOfPdfValues->Fill(xsec) ;
				//histCounter++ ;
				//if( histCounter > 10000 ) {
				//      histOfPdfValues->Draw() ;
				//      c0->Update() ;
				//      c0->SaveAs( "histOfPdfValues-from-Evaluate.eps" ) ;
				//      histCounter = 0 ;
				//}
				break;
		}
	}
	else if( _usePlotAllComponents && _usePlotComponents )
	{
		return 0.;
	}
	else
	{
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
	}

	Observable* timeObs = _datapoint->GetObservable( timeName );
	//if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue( timeObs, timeOffset );
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue( timeObs, 0.0 );
	if( DebugFlag_v7 )
	{
		if( xsec < 0)
		{
			this->DebugPrintXsec( " Bs2JpsiPhi_Signal_v7_v1::diffXsec( ) : return value < 0 = ", xsec );
		}
	}

	return xsec;
}


//...................................
// Integral over angles only for a fixed time.

double Bs2JpsiPhi_Signal_v7::diffXsecTimeOnly()
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

	Observable* timeObs = _datapoint->GetObservable( timeName );
	if( useTimeAcceptance() ) xsec = xsec * timeAcc->getValue( timeObs, 0.0 );

	if( DebugFlag_v7 )
	{
		if( xsec < 0 )
		{
			this->DebugPrintXsec( " Bs2JpsiPhi_Signal_v7_v1::diffXsecTimeOnly( ) : return value < 0 = ", xsec );
		}
	}

	return xsec;
}




//...................................
// Integral over all variables: t + angles

double Bs2JpsiPhi_Signal_v7::diffXsecNorm1()
{
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

	if( DebugFlag_v7 )
	{
		if( norm < 0 )
		{
			this->DebugPrintNorm( " Bs2JpsiPhi_Signal_v7_v1::diffXsecNorm1( ) : return value < 0 = ", norm );
		}
	}
	return norm ;
}

void Bs2JpsiPhi_Signal_v7::generateTimeIntegrals()
{
	double stored_gammaLInt=0.;
	double stored_gammaHInt=0.;

	double tlo_boundary = tlo;
	double thi_boundary = thi;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		timeBinNum = islice;
		AcceptanceSlice* thisSlice = timeAcc->getSlice(islice);

		const double slice_lo = thisSlice->tlow();
		const double slice_hi = thisSlice->thigh();

		tlo = tlo_boundary > slice_lo ? tlo_boundary : slice_lo;
		thi = thi_boundary < slice_hi ? thi_boundary : slice_hi;
		if( thi > tlo )
		{
			this->preCalculateTimeIntegrals();
			stored_gammaLInt += intExpL_stored * thisSlice->height();
			stored_gammaHInt += intExpH_stored * thisSlice->height();
		}
	}
	tlo = tlo_boundary;
	thi = thi_boundary;

	intExpL_stored = stored_gammaLInt;
	intExpH_stored = stored_gammaHInt;
}

void Bs2JpsiPhi_Signal_v7::generateSinusoidIntegrals()
{
	double stored_ExpSinInt=0.;
	double stored_ExpCosInt=0.;

	double tlo_boundary = tlo;
	double thi_boundary = thi;

	for( unsigned int islice = 0; islice < (unsigned) timeAcc->numberOfSlices(); ++islice )
	{
		timeBinNum = islice;
		AcceptanceSlice* thisSlice = timeAcc->getSlice(islice);

		const double slice_lo = thisSlice->tlow();
		const double slice_hi = thisSlice->thigh();

		tlo = tlo_boundary > slice_lo ? tlo_boundary : slice_lo;
		thi = thi_boundary < slice_hi ? thi_boundary : slice_hi;
		if( thi > tlo )
		{
			this->preCalculateSinusoidIntegrals();
			stored_ExpSinInt += intExpSin_stored * thisSlice->height();
			stored_ExpCosInt += intExpCos_stored * thisSlice->height();
		}
	}
	tlo = tlo_boundary;
	thi = thi_boundary;

	intExpSin_stored = stored_ExpSinInt;
	intExpCos_stored = stored_ExpCosInt;
}

void Bs2JpsiPhi_Signal_v7::ConstructTimeIntegrals()
{
	//this->generateIntegrals();

	vector<double> perEventData = _datapoint->GetPerEventData();

	if( perEventData.empty() || !resolutionModel->CacheValid() )
	{
		this->generateTimeIntegrals();
		this->generateSinusoidIntegrals();
		perEventData.clear();
		perEventData.push_back( _gamma );
		perEventData.push_back( dgam );
		perEventData.push_back( delta_ms );
		perEventData.push_back( intExpL_stored );
		perEventData.push_back( intExpH_stored );
		perEventData.push_back( intExpSin_stored );
		perEventData.push_back( intExpCos_stored );
		_datapoint->SetPerEventData( perEventData );
	}
	else if( _gamma != perEventData[0] )
	{
		this->generateTimeIntegrals();
		this->generateSinusoidIntegrals();
		perEventData.clear();
		perEventData.push_back( _gamma );
		perEventData.push_back( dgam );
		perEventData.push_back( delta_ms );
		perEventData.push_back( intExpL_stored );
		perEventData.push_back( intExpH_stored );
		perEventData.push_back( intExpSin_stored );
		perEventData.push_back( intExpCos_stored );
		_datapoint->SetPerEventData( perEventData );
	}
	else if( dgam != perEventData[1] )
	{
		intExpSin_stored = perEventData[5];
		intExpCos_stored = perEventData[6];
		this->generateTimeIntegrals();
		perEventData.clear();
		perEventData.push_back( _gamma );
		perEventData.push_back( dgam );
		perEventData.push_back( delta_ms );
		perEventData.push_back( intExpL_stored );
		perEventData.push_back( intExpH_stored );
		perEventData.push_back( intExpSin_stored );
		perEventData.push_back( intExpCos_stored );
		_datapoint->SetPerEventData( perEventData );
	}
	else if( delta_ms != perEventData[2] )
	{
		intExpL_stored = perEventData[3];
		intExpH_stored = perEventData[4];
		this->generateSinusoidIntegrals();
		perEventData.clear();
		perEventData.push_back( _gamma );
		perEventData.push_back( dgam );
		perEventData.push_back( delta_ms );
		perEventData.push_back( intExpL_stored );
		perEventData.push_back( intExpH_stored );
		perEventData.push_back( intExpSin_stored );
		perEventData.push_back( intExpCos_stored );
		_datapoint->SetPerEventData( perEventData );
	}
	else
	{
		intExpL_stored = perEventData[3];
		intExpH_stored = perEventData[4];
		intExpSin_stored = perEventData[5];
		intExpCos_stored = perEventData[6];
	}
}

//....................................................
// New method to calculate normalisation using a histogrammed "low-end" time acceptance function
// The acceptance function information is all contained in the timeAcceptance member object,
double Bs2JpsiPhi_Signal_v7::diffXsecCompositeNorm1(  )
{
	double returnValue = 0;

	this->ConstructTimeIntegrals();

	returnValue = this->diffXsecNorm1();

	return returnValue;
}


//.......................................................
// New speed up method to Cache time integrals
void Bs2JpsiPhi_Signal_v7::CacheAmplitudesAndAngles()
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

}

//....................................................
// New to prepare all of the coeefficients needed in the time dependen terms
void Bs2JpsiPhi_Signal_v7::prepareCDS( double lambdaVal, double PhisVal )
{
	double lambda_sq = lambdaVal*lambdaVal;
	double inv_lambda = 1./(1.0 + lambda_sq);

	double F1 = 2.0*lambdaVal *inv_lambda;
	double F2 = (1.0 - lambda_sq) *inv_lambda;

	if( !_useCosAndSin )
	{
		_sinphis = sin(PhisVal);
		_cosphis = cos(PhisVal);
	}
	_SS = _sinphis * F1;
	_DD = _cosphis * F1;
	_CC = F2;

}


void Bs2JpsiPhi_Signal_v7::DebugPrintEvaluate( string message, double value )
{
	PDF_THREAD_LOCK

		cout << "1: " << " A: " << CachedA1 << " * " << " T: " << timeFactorA0A0(  ) << " = " << CachedA1*timeFactorA0A0(  ) << endl;
	cout << "2: " << " A: " << CachedA2 << " * " << " T: " << timeFactorAPAP(  ) << " = " << CachedA2*timeFactorAPAP(  ) << endl;
	cout << "3: " << " A: " << CachedA3 << " * " << " T: " << timeFactorATAT(  ) << " = " << CachedA3*timeFactorATAT(  ) << endl;

	cout << "4: " << " A: " << CachedA4 << " * " << " T: " << timeFactorImAPAT(  ) << " = " << CachedA4*timeFactorImAPAT(  ) << endl;
	cout << "5: " << " A: " << CachedA5 << " * " << " T: " << timeFactorReA0AP(  ) << " = " << CachedA5*timeFactorReA0AP(  ) << endl;
	cout << "6: " << " A: " << CachedA6 << " * " << " T: " << timeFactorImA0AT(  ) << " = " << CachedA6*timeFactorImA0AT(  ) << endl;

	cout << "7: " << " A: " << CachedA7 << " * " << " T: " << timeFactorASAS(  ) << " = " << CachedA7*timeFactorASAS(  ) << endl;

	cout << "8: " << " A: " << CachedA8 << " * " << " T: " << timeFactorReASAP(  ) << " = " << CachedA8*timeFactorReASAP(  ) << endl;
	cout << "9: " << " A: " << CachedA9 << " * " << " T: " << timeFactorImASAT(  ) << " = " << CachedA9*timeFactorImASAT(  ) << endl;
	cout << "10:" << " A: " << CachedA10 << " * " << " T: " << timeFactorReASA0(  ) << " = " << CachedA10*timeFactorReASA0(  ) << endl;

	cout << endl;
	this->preCalculateTimeIntegralsDebug();

	this->DebugPrint( message, value );
	PDF_THREAD_UNLOCK
}

void Bs2JpsiPhi_Signal_v7::DebugPrintNormalisation( string message, double value )
{
	PDF_THREAD_LOCK

		cout << "1: " << " A: " << A0()*A0() << " * " << " T: " << timeFactorA0A0Int(  ) << " * " << " Acc: " <<  angAccI1 << " = " << A0()*A0()*timeFactorA0A0Int(  ) * angAccI1 << endl;
	cout << "2: " << " A: " << AP()*AP() << " * " << " T: " << timeFactorAPAPInt(  ) << " * " << " Acc: " <<  angAccI2 << " = " << AP()*AP()*timeFactorAPAPInt(  ) * angAccI2 << endl;
	cout << "3: " << " A: " << AT()*AT() << " * " << " T: " << timeFactorATATInt(  ) << " * " << " Acc: " <<  angAccI3 << " = " << AT()*AT()*timeFactorATATInt(  ) * angAccI3 << endl;

	cout << "4: " << " A: " << AP()*AT() << " * " << " T: " << timeFactorImAPATInt(  ) << " * " << " Acc: " <<  angAccI4 << " = " << AP()*AT()*timeFactorImAPATInt(  ) * angAccI4 << endl;
	cout << "5: " << " A: " << A0()*AP() << " * " << " T: " << timeFactorReA0APInt(  ) << " * " << " Acc: " <<  angAccI5 << " = " << A0()*AP()*timeFactorReA0APInt(  ) * angAccI5 << endl;
	cout << "6: " << " A: " << A0()*AT() << " * " << " T: " << timeFactorImA0ATInt(  ) << " * " << " Acc: " <<  angAccI6 << " = " << A0()*AT()*timeFactorImA0ATInt(  ) * angAccI6 << endl;

	cout << "7: " << " A: " << AS()*AS() << " * " << " T: " << timeFactorASASInt(  ) << " * " << " Acc: " <<  angAccI7 << " = " << AS()*AS()*timeFactorASASInt(  ) * angAccI7 << endl;

	cout << "8: " << " A: " << ASint()*AP() << " * " << " T: " << timeFactorReASAPInt(  ) << " * " << " Acc: " <<  angAccI8 << " = " << ASint()*AP()*timeFactorReASAPInt(  ) * angAccI8 << endl;
	cout << "9: " << " A: " << ASint()*AT() << " * " << " T: " << timeFactorImASATInt(  ) << " * " << " Acc: " <<  angAccI9 << " = " << ASint()*AT()*timeFactorImASATInt(  ) * angAccI9 << endl;
	cout << "10:" << " A: " << ASint()*A0() << " * " << " T: " << timeFactorReASA0Int(  ) << " * " << " Acc: " <<  angAccI10 << " = " << ASint()*A0()*timeFactorReASA0Int(  ) * angAccI10 << endl;

	cout << endl;
	this->preCalculateTimeIntegralsDebug();

	this->DebugPrint( message, value );
	PDF_THREAD_UNLOCK
}

//===========================================================================================
// Debug printout
//===========================================================================================
void Bs2JpsiPhi_Signal_v7::DebugPrint( string message, double value )  const
{
	if( !performingComponentProjection )
	{
		(void) message; (void) value;
		cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_Signal_v7::DebugPrint ***************************" << endl ;
		cout << message << value << endl <<endl;

		cout << endl ;
		cout << "   gamma " << gamma() << endl;
		cout << "   gl    " << gamma_l() << endl;
		cout << "   gh    " << gamma_h()  << endl;
		cout << "   AT^2    " << AT()*AT() << endl;
		cout << "   AP^2    " << AP()*AP() << endl;
		cout << "   A0^2    " << A0()*A0() << endl;
		cout << "   AS^2    " << AS()*AS() << endl;
		cout << "   P-Wave  " << A0()*A0()+AP()*AP()+AT()*AT() << endl;
		cout << "   ATOTAL  " << AS()*AS()+A0()*A0()+AP()*AP()+AT()*AT() << endl;
		cout << "   delta_ms       " << delta_ms << endl;
		//cout << "   mistag         " << mistag() << endl;
		//cout << "   mistagP1       " << _mistagP1 << endl;
		//cout << "   mistagP0       " << _mistagP0 << endl;
		//cout << "   mistagSetPoint " << _mistagSetPoint << endl;
		cout << " For event with:" << endl;
		//cout << "   time      " << t << endl ;
		//cout << "   ctheta_tr " << ctheta_tr << endl ;
		//cout << "   ctheta_1 " << ctheta_1 << endl ;
		//cout << "   phi_tr " << phi_tr << endl ;
		if( _datapoint ) _datapoint->Print();
		cout << endl << " For ParameterSet:" << endl;
		allParameters.Print();
	}
}


void Bs2JpsiPhi_Signal_v7::DebugPrintXsec( string message, double value )
{
	PDF_THREAD_LOCK

		if( !performingComponentProjection )
		{
			(void) message; (void) value;
			cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_Signal_v7::DebugPrintXsec ***************************" << endl;
			cout << message << value << endl << endl;
			cout << "   A0()*A0() term: " <<  A0()*A0() * timeFactorA0A0(  ) * A0A0_value << endl;
			cout << "   AP()*AP() term: " <<AP()*AP() * timeFactorAPAP(  ) * APAP_value << endl;
			cout << "   AT()*AT() term: " <<AT()*AT() * timeFactorATAT(  ) * ATAT_value << endl << endl;

			cout << "   AP()*AT() term: " <<AP()*AT() * timeFactorImAPAT(  ) * ImAPAT_value << endl;
			cout << "                 : " <<AP()*AT() <<" / "<<  timeFactorImAPAT( )  <<" / "<<  ImAPAT_value << endl;
			cout << "   A0()*AP() term: " <<A0()*AP() * timeFactorReA0AP(  ) * ReA0AP_value << endl;
			cout << "                 : " <<A0()*AP() <<" / "<<  timeFactorReA0AP(  ) <<" / "<<  ReA0AP_value << endl;
			cout << "   A0()*AT() term: " <<A0()*AT() * timeFactorImA0AT(  ) * ImA0AT_value << endl << endl;
			cout << "                 : " <<A0()*AT() <<" / "<<  timeFactorImA0AT(  ) <<" / "<<  ImA0AT_value << endl << endl;

			cout << "   AS()*AS() term: " <<AS()*AS() * timeFactorASAS(  ) * ASAS_value << endl << endl;

			cout << "   AS()*AP() term: " <<AS()*AP() * timeFactorReASAP(  ) * ReASAP_value << endl;
			cout << "                 : " <<AS()*AP() <<" / "<<   timeFactorReASAP(  ) <<" / "<<  ReASAP_value << endl;
			cout << "   AS()*AT() term: " <<AS()*AT() * timeFactorImASAT(  ) * ImASAT_value << endl;
			cout << "                 : " <<AS()*AT() <<" / "<<   timeFactorImASAT(  ) <<" / "<<   ImASAT_value << endl;
			cout << "   AS()*A0() term: " <<AS()*A0() * timeFactorReASA0(  ) * ReASA0_value<< endl;
			cout << "                 : " <<AS()*A0() <<" / "<<   timeFactorReASA0(  ) <<" / "<<  ReASA0_value << endl << endl;

			double PwaveTot =
				A0()*A0() * timeFactorA0A0(  ) * A0A0_value +
				AP()*AP() * timeFactorAPAP(  ) * APAP_value +
				AT()*AT() * timeFactorATAT(  ) * ATAT_value +
				AP()*AT() * timeFactorImAPAT(  ) * ImAPAT_value +
				A0()*AP() * timeFactorReA0AP(  ) * ReA0AP_value +
				A0()*AT() * timeFactorImA0AT(  ) * ImA0AT_value;

			double SwaveAdditions =
				AS()*AS() * timeFactorASAS(  ) * ASAS_value +
				AS()*AP() * timeFactorReASAP(  ) * ReASAP_value +
				AS()*AT() * timeFactorImASAT(  ) * ImASAT_value +
				AS()*A0() * timeFactorReASA0(  ) * ReASA0_value;

			cout << "   Pwave Only : " << PwaveTot << endl;
			cout << "   Swave add : " <<  SwaveAdditions << endl;
			if( _datapoint ) _datapoint->Print();
		}

	PDF_THREAD_UNLOCK
}

void Bs2JpsiPhi_Signal_v7::DebugPrintNorm( string message, double value )
{
	PDF_THREAD_LOCK

		if( !performingComponentProjection )
		{
			(void) message; (void) value;
			cout << "*************DEBUG OUTPUT FROM Bs2JpsiPhi_Signal_v7::DebugPrintNorm ***************************" << endl;
			cout << message << value << endl << endl;

			cout << endl;
			cout <<  A0()*A0() * timeFactorA0A0Int(  )* angAccI1 << endl;
			cout <<  AP()*AP() * timeFactorAPAPInt(  )* angAccI2 << endl;
			cout <<  AT()*AT() * timeFactorATATInt(  )* angAccI3 << endl << endl;
			cout <<  AP()*AT() * timeFactorImAPATInt(  ) * angAccI4 << endl;
			cout <<  A0()*AP() * timeFactorReA0APInt(  ) * angAccI5 << endl;
			cout <<  A0()*AT() * timeFactorImA0ATInt(  ) * angAccI6 << endl << endl;
			cout <<  AS()*AS() * timeFactorASASInt(  ) * angAccI7 << endl;
			cout <<  AS()*AP() * timeFactorReASAPInt(  ) * angAccI8<< endl;
			cout <<  AS()*AT() * timeFactorImASATInt(  ) * angAccI9<< endl;
			cout <<  AS()*A0() * timeFactorReASA0Int(  ) * angAccI10<< endl;
		}

	PDF_THREAD_UNLOCK
}

double Bs2JpsiPhi_Signal_v7::timeFactorEven()  const
{
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( 1.0 + cosphis() ) * expL()
					+ ( 1.0 - cosphis() ) * expH()
					) +
			_mistagCalibModel->D2() * (
					( 2.0 * sinphis() ) * expSin( )
					+ ( 2.0 * CC()      ) * expCos( )
					);
	}
	else
	{
		return ( 1.0 + cosphis() ) * expL() + ( 1.0 - cosphis() ) * expH();
	}
}

void Bs2JpsiPhi_Signal_v7::timeFactorEvenDebug() const
{
	cout << "D1: " << _mistagCalibModel->D1() << " * ( " << ( 1.0 + cosphis() ) * expL( ) << " + " << ( 1.0 + cosphis() ) * expH( ) << " )" << endl;
	cout << "D2: " << _mistagCalibModel->D2() << " * ( " << ( 2.0 * sinphis() ) * expSin( ) << " + " << + ( 2.0 * CC()      ) * expCos( ) << " )" << endl;
}

double Bs2JpsiPhi_Signal_v7::timeFactorEvenInt()  const
{
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( 1.0 + cosphis() )  * intExpL()
					+ ( 1.0 - cosphis() )  * intExpH()
					) +
			_mistagCalibModel->D2() * (
					( 2.0 * sinphis() ) * intExpSin( )
					+  ( 2.0 * CC()      ) * intExpCos( )
					);
	}
	else
	{
		return ( 1.0 + cosphis() ) * intExpL() + ( 1.0 - cosphis() ) * intExpH();
	}
}


void Bs2JpsiPhi_Signal_v7::timeFactorEvenIntDebug() const
{
	cout << "D1: " << _mistagCalibModel->D1() << "  intExpL(): " << intExpL() << "  intExpH(): " << intExpH() << endl;
	cout << "D2: " << _mistagCalibModel->D2() << "  intExpSin(): " << intExpSin( ) << "  intExpCos(): " << intExpSin( ) << endl;
	cout << "   " << _mistagCalibModel->D1() *( ( 1.0 + cosphis() )  * intExpL() + ( 1.0 - cosphis() )  * intExpH() ) << endl;
	cout << " + " << _mistagCalibModel->D2() *( ( 2.0 * sinphis() ) * intExpSin( ) + ( 2.0 * CC()      ) * intExpCos( ) )<< endl;
}


double Bs2JpsiPhi_Signal_v7::timeFactorOdd(  )   const
{
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( 1.0 - cosphis() ) * expL()
					+ ( 1.0 + cosphis() ) * expH()
					) +
			_mistagCalibModel->D2() * (
					-  ( 2.0 * sinphis() ) * expSin( )
					+  ( 2.0 * CC()      ) * expCos( )
					);
	}
	else
	{
		return ( 1.0 - cosphis() ) * expL() + ( 1.0 + cosphis() ) * expH();
	}
}

double Bs2JpsiPhi_Signal_v7::timeFactorOddInt(  )  const
{
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( 1.0 - cosphis() ) * intExpL()
					+ ( 1.0 + cosphis() ) * intExpH()
					) +
			_mistagCalibModel->D2() * (
					-  ( 2.0 * sinphis() ) * intExpSin( )
					+  ( 2.0 * CC()      ) * intExpCos( )
					);
	}
	else
	{
		return ( 1.0 - cosphis() ) * intExpL() + ( 1.0 + cosphis() ) * intExpH();
	}
}

//...........................
double Bs2JpsiPhi_Signal_v7::timeFactorA0A0( )
{
	this->prepareCDS( lambda_zeroVal, phis_zeroVal );
	return timeFactorEven( );
}

double Bs2JpsiPhi_Signal_v7::timeFactorA0A0Int( )
{
	this->prepareCDS( lambda_zeroVal, phis_zeroVal );
	return timeFactorEvenInt( );
}

//...........................
double Bs2JpsiPhi_Signal_v7::timeFactorAPAP( )
{
	this->prepareCDS( lambda_paraVal, phis_paraVal );
	return timeFactorEven( );
}
double Bs2JpsiPhi_Signal_v7::timeFactorAPAPInt( )
{
	this->prepareCDS( lambda_paraVal, phis_paraVal );
	return timeFactorEvenInt( );
}

//...........................
double Bs2JpsiPhi_Signal_v7::timeFactorATAT( )
{
	this->prepareCDS( lambda_perpVal, phis_perpVal );
	return timeFactorOdd( );
}

double Bs2JpsiPhi_Signal_v7::timeFactorATATInt( )
{
	this->prepareCDS( lambda_perpVal, phis_perpVal );
	return timeFactorOddInt( );
}

double Bs2JpsiPhi_Signal_v7::timeFactorImAPAT( )
{
	this->prepareCDS( sqrt(lambda_paraVal*lambda_perpVal), 0.5*(phis_paraVal+phis_perpVal) );

	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( expL( ) - expH( ) ) * cos_delta1 * sinphis()
					+ ( expL( ) + expH( ) ) * sin_delta1 * CC()
					) +
			_mistagCalibModel->D2() * (
					2.0  * ( sin_delta1*expCos( ) - cos_delta1*cosphis()*expSin( ) )
					);
	}
	else
	{
		return ( expL( ) - expH( ) ) * cos_delta1 * sinphis() + ( expL( ) + expH( ) ) * sin_delta1 * CC();
	}
}

double Bs2JpsiPhi_Signal_v7::timeFactorImAPATInt( )
{
	this->prepareCDS( sqrt(lambda_paraVal*lambda_perpVal), 0.5*(phis_paraVal+phis_perpVal) );
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( intExpL() - intExpH() ) * cos_delta1 * sinphis()
					+ ( intExpL() + intExpH() ) * sin_delta1 * CC()
					) +
			_mistagCalibModel->D2() * (
					2.0  * ( sin_delta1*intExpCos() - cos_delta1*cosphis()*intExpSin() )
					);
	}
	else
	{
		return ( intExpL() - intExpH() ) * cos_delta1 * sinphis() +  ( intExpL() + intExpH() ) * sin_delta1 * CC();
	}
}

double Bs2JpsiPhi_Signal_v7::timeFactorReA0AP( )
{
	this->prepareCDS( sqrt(lambda_zeroVal*lambda_perpVal), 0.5*(phis_zeroVal+phis_perpVal) );
	if( _useCosDpar ) return cosdpar * this->timeFactorEven(  );//PELC-COSDPAR Special for fitting cosdpar separately
	else return cos_delta_2_1 * this->timeFactorEven(  );
}

double Bs2JpsiPhi_Signal_v7::timeFactorReA0APInt( )
{
	this->prepareCDS( sqrt(lambda_zeroVal*lambda_perpVal), 0.5*(phis_zeroVal+phis_perpVal) );
	if( _useCosDpar ) return cosdpar * this->timeFactorEvenInt( );//PELC-COSDPAR Special for fitting cosdpar separately
	else return cos_delta_2_1 * this->timeFactorEvenInt( );
}

double Bs2JpsiPhi_Signal_v7::timeFactorImA0AT(  )
{
	this->prepareCDS( sqrt(lambda_zeroVal*lambda_paraVal), 0.5*(phis_zeroVal+phis_paraVal) );
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( expL( ) - expH( ) ) * cos_delta2 * sinphis()
					+ ( expL( ) + expH( ) ) * sin_delta2 * CC()
					) +
			_mistagCalibModel->D2() * (
					2.0  * ( sin_delta2*expCos( ) - cos_delta2*cosphis()*expSin( ) )
					);
	}
	else
	{
		return ( expL( ) - expH( ) ) * cos_delta2 * sinphis() + ( expL( ) + expH( ) ) * sin_delta2 * CC();
	}
}

double Bs2JpsiPhi_Signal_v7::timeFactorImA0ATInt( )
{
	this->prepareCDS( sqrt(lambda_zeroVal*lambda_paraVal), 0.5*(phis_zeroVal+phis_paraVal) );
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( intExpL() - intExpH()  ) * cos_delta2 * sinphis()
					+ ( intExpL() + intExpH()  ) * sin_delta2 * CC()
					) +
			_mistagCalibModel->D2() * (
					2.0  * ( sin_delta2*intExpCos() - cos_delta2*cosphis()*intExpSin()  )
					);
	}
	else
	{
		return ( intExpL() - intExpH()  ) * cos_delta2 * sinphis() + ( intExpL() + intExpH()  ) * sin_delta2 * CC();
	}
}


double Bs2JpsiPhi_Signal_v7::timeFactorASAS( )
{
	this->prepareCDS( lambda_SVal, phis_SVal );
	return timeFactorOdd( );
}
double Bs2JpsiPhi_Signal_v7::timeFactorASASInt( )
{
	this->prepareCDS( lambda_SVal, phis_SVal );
	return timeFactorOddInt( );
}


double Bs2JpsiPhi_Signal_v7::timeFactorReASAP( )
{
	this->prepareCDS( sqrt(lambda_SVal*lambda_perpVal), 0.5*(phis_SVal+phis_perpVal) );
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( expL( ) - expH( ) ) * sin_delta_para_s * sinphis()
					+ ( expL( ) + expH( ) ) * cos_delta_para_s * CC()
					) +
			_mistagCalibModel->D2() * (
					2.0  * ( cos_delta_para_s*expCos( ) - sin_delta_para_s*cosphis()*expSin( ) )
					);
	}
	else
	{
		return ( expL( ) - expH( ) ) * sin_delta_para_s * sinphis() + ( expL( ) + expH( ) ) * cos_delta_para_s * CC();
	}
}

double Bs2JpsiPhi_Signal_v7::timeFactorReASAPInt( )
{
	this->prepareCDS( sqrt(lambda_SVal*lambda_perpVal), 0.5*(phis_SVal+phis_perpVal) );
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( intExpL() - intExpH() ) * sin_delta_para_s * sinphis()
					+ ( intExpL() + intExpH() ) * cos_delta_para_s * CC()
					) +
			_mistagCalibModel->D2() * (
					2.0  * ( cos_delta_para_s*intExpCos() - sin_delta_para_s*cosphis()*intExpSin() )
					);
	}
	else
	{
		return ( intExpL() - intExpH() ) * sin_delta_para_s * sinphis() + ( intExpL() + intExpH() ) * cos_delta_para_s * CC();
	}
}


double Bs2JpsiPhi_Signal_v7::timeFactorImASAT( )
{
	this->prepareCDS( sqrt(lambda_SVal*lambda_paraVal), 0.5*(phis_SVal+phis_paraVal) );
	return sin_delta_perp_s * this->timeFactorOdd(  );
}

double Bs2JpsiPhi_Signal_v7::timeFactorImASATInt( )
{
	this->prepareCDS( sqrt(lambda_SVal*lambda_paraVal), 0.5*(phis_SVal+phis_paraVal) );
	return sin_delta_perp_s * this->timeFactorOddInt( );
}

double Bs2JpsiPhi_Signal_v7::timeFactorReASA0( )
{
	this->prepareCDS( sqrt(lambda_SVal*lambda_zeroVal), 0.5*(phis_SVal+phis_zeroVal) );
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( expL( ) - expH( ) ) * sin_delta_zero_s * sinphis()
					+ ( expL( ) + expH( ) ) * cos_delta_zero_s * CC()
					) +
			_mistagCalibModel->D2() * (
					2.0  * ( cos_delta_zero_s*expCos( ) - sin_delta_zero_s*cosphis()*expSin( ) )
					);
	}
	else
	{
		return ( expL( ) - expH( ) ) * sin_delta_zero_s * sinphis() + ( expL( ) + expH( ) ) * cos_delta_zero_s * CC();
	}
}

double Bs2JpsiPhi_Signal_v7::timeFactorReASA0Int( )
{
	this->prepareCDS( sqrt(lambda_SVal*lambda_zeroVal), 0.5*(phis_SVal+phis_zeroVal) );
	if( _eventIsTagged )
	{
		return
			_mistagCalibModel->D1() * (
					( intExpL() - intExpH() ) * sin_delta_zero_s * sinphis()
					+ ( intExpL() + intExpH() ) * cos_delta_zero_s * CC()
					) +
			_mistagCalibModel->D2() * (
					2.0  * ( cos_delta_zero_s*intExpCos() - sin_delta_zero_s*cosphis()*intExpSin() )
					);
	}
	else
	{
		return ( intExpL() - intExpH() ) * sin_delta_zero_s * sinphis() + ( intExpL() + intExpH() ) * cos_delta_zero_s * CC();
	}
}

