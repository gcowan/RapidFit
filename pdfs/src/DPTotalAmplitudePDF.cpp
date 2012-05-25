// $Id: DPTotalAmplitudePDF.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class DPTotalAmplitudePDF DPTotalAmplitudePDF.cpp
 *
 *  PDF for Bd2JpsiKpi spin-0
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#include "DPTotalAmplitudePDF.h"
#include "DPJpsiKaon.hh"
#include "DPZplusK.hh"
#include "DPHelpers.hh"
#include "DPComponent.hh"

#include <iostream>
#include "math.h"
#include "TMath.h"
#include "TComplex.h"
#include "RooMath.h"

PDF_CREATOR( DPTotalAmplitudePDF );

//Constructor
DPTotalAmplitudePDF::DPTotalAmplitudePDF( PDFConfigurator* configurator) :

	// Physics parameters
	  fracA0sqZplusName	( configurator->getName("fracA0sqZplus") )
	, fracApsqZplusName	( configurator->getName("fracApsqZplus") )
	, fracZplusName		( configurator->getName("fracZplus") )
	, phaseA0ZplusName	( configurator->getName("phaseA0Zplus") )
	, phaseApZplusName	( configurator->getName("phaseApZplus") )
	, phaseAmZplusName	( configurator->getName("phaseAmZplus") )
	
	, fracA0sqKst892Name	( configurator->getName("fracA0sqKst892") )
	, fracApsqKst892Name	( configurator->getName("fracApsqKst892") )
	//, fracKst892Name	( configurator->getName("fracKst892") )
	, phaseA0Kst892Name	( configurator->getName("phaseA0Kst892") )
	, phaseApKst892Name	( configurator->getName("phaseApKst892") )
	, phaseAmKst892Name	( configurator->getName("phaseAmKst892") )
	
	, fracA0sqKst1410Name	( configurator->getName("fracA0sqKst1410") )
	, fracApsqKst1410Name	( configurator->getName("fracApsqKst1410") )
	, fracKst1410Name	( configurator->getName("fracKst1410") )
	, phaseA0Kst1410Name	( configurator->getName("phaseA0Kst1410") )
	, phaseApKst1410Name	( configurator->getName("phaseApKst1410") )
	, phaseAmKst1410Name	( configurator->getName("phaseAmKst1410") )
	
	, fracA0sqKst1680Name	( configurator->getName("fracA0sqKst1680") )
	, fracApsqKst1680Name	( configurator->getName("fracApsqKst1680") )
	, fracKst1680Name	( configurator->getName("fracKst1680") )
	, phaseA0Kst1680Name	( configurator->getName("phaseA0Kst1680") )
	, phaseApKst1680Name	( configurator->getName("phaseApKst1680") )
	, phaseAmKst1680Name	( configurator->getName("phaseAmKst1680") )
	
	, fracK01430Name	( configurator->getName("fracK01430") )
	, phaseA0K01430Name	( configurator->getName("phaseA0K01430") )
	
	, fracK21430Name	( configurator->getName("fracK21430") )
	, phaseA0K21430Name	( configurator->getName("phaseA0K21430") )
	
	, massZplusName		( configurator->getName("massZplus") )
	, widthZplusName	( configurator->getName("widthZplus") )
	, massKst892Name	( configurator->getName("massKst892") )
	, widthKst892Name	( configurator->getName("widthKst892") )
	, massKst1410Name	( configurator->getName("massKst1410") )
	, widthKst1410Name	( configurator->getName("widthKst1410") )
	, massKst1680Name	( configurator->getName("massKst1680") )
	, widthKst1680Name	( configurator->getName("widthKst1680") )
	, massK01430Name	( configurator->getName("massK01430") )
	, widthK01430Name	( configurator->getName("widthK01430") )
	, massK21430Name	( configurator->getName("massK21430") )
	, widthK21430Name	( configurator->getName("widthK21430") )
	// Observables
	, m23Name	( configurator->getName("m23") )
	, cosTheta1Name	( configurator->getName("cosTheta1") )
	, cosTheta2Name	( configurator->getName("cosTheta2") )
	, phiName	( configurator->getName("phi") )
	// The actual values of the parameters and observables
	, fracA0sqZplus(),  fracApsqZplus(),  fracZplus(),  phaseA0Zplus(),  phaseApZplus(),  phaseAmZplus()
	, fracA0sqKst892(),  fracApsqKst892(),  fracKst892(),  phaseA0Kst892(),  phaseApKst892(),  phaseAmKst892()
	, fracA0sqKst1410(), fracApsqKst1410(), fracKst1410(), phaseA0Kst1410(), phaseApKst1410(), phaseAmKst1410()
	, fracA0sqKst1680(), fracApsqKst1680(), fracKst1680(), phaseA0Kst1680(), phaseApKst1680(), phaseAmKst1680()
	, fracK01430(),  				  phaseA0K01430()
	, fracK21430(),  				  phaseA0K21430()
	, massZplus(), widthZplus()
	, massKst892(), widthKst892()
	, massKst1410(), widthKst1410()
	, massKst1680(), widthKst1680()
	, massK01430(), widthK01430()
	, massK21430(), widthK21430()
	, m23(), cosTheta1(), cosTheta2(), phi()
	, pMuPlus(0., 0., 0., 0.), pMuMinus(0., 0., 0., 0.), pPi(0., 0., 0., 0.), pK(0., 0., 0., 0.), pB(0., 0., 0., 5.279)
	, cosARefs()
, histogramFile(), angularAccHistCosTheta1(), angularAccHistPhi(), angularAccHistMassCosTheta2()
{
	MakePrototypes();

	componentIndex = 0;

	mJpsi=3.096916;
	// Construct all components we need
	DPComponent * tmp;
	// B0 --> J/psi K*
	tmp=new DPJpsiKaon(0, 1, 5.279, 0.89594, 0.0487, 0.493677,
			0.13957018, 5.0, 1.5, 3.096916,1);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K*(1410)
	tmp=new DPJpsiKaon(0, 1, 5.279, 1.414, 0.232, 0.493677,
			0.13957018, 5.0, 1.5, 3.096916,1);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K*(1680)
	tmp=new DPJpsiKaon(0, 1, 5.279, 1.717, 0.322, 0.493677,
			0.13957018, 5.0, 1.5, 3.096916,1);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K0(1430)
	tmp=new DPJpsiKaon(0, 0, 5.279, 1.425, 0.270, 0.493677,
			0.13957018, 5.0, 1.5, 3.096916,0);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K2(1430)
	tmp=new DPJpsiKaon(0, 0, 5.279, 1.4324, 0.109, 0.493677,
			0.13957018, 5.0, 1.5, 3.096916,2);
	KpiComponents.push_back(tmp);
  
	// Kpi s-wave using LASS
  	//tmp=new DPJpsiKaon(0, 0, 5.279, 1.425, 0.270, 0.493677,
        //             0.13957018, 5.0, 1.5, 3.096916,0,
        //             "LASS", 0.00415, 0.00136);
   	//KpiComponents.push_back(tmp);

	// B0 --> Z+ K-
	tmp=new DPZplusK(0,0,5.279,4.430,0.100,0.493677,
			0.13957018, 5.0, 1.5, 3.096916,0);
	//ZComponents.push_back(tmp);

	this->SetNumericalNormalisation( true );
	this->TurnCachingOff();
	useAngularAcceptance = false;
	if ( configurator->isTrue( "UseAngularAcceptance" ) )
	{
		useAngularAcceptance = true;
		histogramFile = TFile::Open("~gcowan/public/RapidFit/pdfs/configdata/Bd2JpsiKpi_TotalEff.root");
		angularAccHistCosTheta1 = (TH1D*)histogramFile->Get("cosmu_effTot");
		angularAccHistPhi 	= (TH1D*)histogramFile->Get("delta_phi_effTot");
		angularAccHistMassCosTheta2 = (TH2D*)histogramFile->Get("mass_cos_effTot");
	}
}

DPTotalAmplitudePDF::DPTotalAmplitudePDF( const DPTotalAmplitudePDF &copy ) :
	BasePDF( (BasePDF)copy )
	,m23Name(copy.m23Name)
	,cosTheta1Name(copy.cosTheta1Name)
	,cosTheta2Name(copy.cosTheta2Name)
	,phiName(copy.phiName)
	,m23(copy.m23)
	,cosTheta1(copy.cosTheta1)
	,cosTheta2(copy.cosTheta2)
	,phi(copy.phi)
	,fracA0sqZplusName(copy.fracA0sqZplusName)
	,fracApsqZplusName(copy.fracApsqZplusName)
	,fracZplusName(copy.fracZplusName)
	,phaseA0ZplusName(copy.phaseA0ZplusName)
	,phaseApZplusName(copy.phaseApZplusName)
	,phaseAmZplusName(copy.phaseAmZplusName)
	
	,fracA0sqKst892Name(copy.fracA0sqKst892Name)
	,fracApsqKst892Name(copy.fracApsqKst892Name)
	//,fracKst892Name(copy.fracKst892Name)
	,phaseA0Kst892Name(copy.phaseA0Kst892Name)
	,phaseApKst892Name(copy.phaseApKst892Name)
	,phaseAmKst892Name(copy.phaseAmKst892Name)
	
	,fracA0sqKst1410Name(copy.fracA0sqKst1410Name)
	,fracApsqKst1410Name(copy.fracApsqKst1410Name)
	,fracKst1410Name(copy.fracKst1410Name)
	,phaseA0Kst1410Name(copy.phaseA0Kst1410Name)
	,phaseApKst1410Name(copy.phaseApKst1410Name)
	,phaseAmKst1410Name(copy.phaseAmKst1410Name)
	
	,fracA0sqKst1680Name(copy.fracA0sqKst1680Name)
	,fracApsqKst1680Name(copy.fracApsqKst1680Name)
	,fracKst1680Name(copy.fracKst1680Name)
	,phaseA0Kst1680Name(copy.phaseA0Kst1680Name)
	,phaseApKst1680Name(copy.phaseApKst1680Name)
	,phaseAmKst1680Name(copy.phaseAmKst1680Name)
	
	,fracK01430Name(copy.fracK01430Name)
	,phaseA0K01430Name(copy.phaseA0K01430Name)
	,fracK21430Name(copy.fracK21430Name)
	,phaseA0K21430Name(copy.phaseA0K21430Name)
	
	,massZplusName(copy.massZplusName)
	,widthZplusName(copy.widthZplusName)
	,massKst892Name(copy.massKst892Name)
	,widthKst892Name(copy.widthKst892Name)
	,massKst1410Name(copy.massKst1410Name)
	,widthKst1410Name(copy.widthKst1410Name)
	,massKst1680Name(copy.massKst1680Name)
	,widthKst1680Name(copy.widthKst1680Name)
	,massK01430Name(copy.massK01430Name)
	,widthK01430Name(copy.widthK01430Name)
	,massK21430Name(copy.massK21430Name)
	,widthK21430Name(copy.widthK21430Name)
	
	,fracA0sqZplus(copy.fracA0sqZplus)
	,fracApsqZplus(copy.fracApsqZplus)
	,fracZplus(copy.fracZplus)
	,phaseA0Zplus(copy.phaseA0Zplus)
	,phaseApZplus(copy.phaseApZplus)
	,phaseAmZplus(copy.phaseAmZplus)
	
	,fracA0sqKst892(copy.fracA0sqKst892)
	,fracApsqKst892(copy.fracApsqKst892)
	,fracKst892(copy.fracKst892)
	,phaseA0Kst892(copy.phaseA0Kst892)
	,phaseApKst892(copy.phaseApKst892)
	,phaseAmKst892(copy.phaseAmKst892)
	
	,fracA0sqKst1410(copy.fracA0sqKst1410)
	,fracApsqKst1410(copy.fracApsqKst1410)
	,fracKst1410(copy.fracKst1410)
	,phaseA0Kst1410(copy.phaseA0Kst1410)
	,phaseApKst1410(copy.phaseApKst1410)
	,phaseAmKst1410(copy.phaseAmKst1410)
	
	,fracA0sqKst1680(copy.fracA0sqKst1680)
	,fracApsqKst1680(copy.fracApsqKst1680)
	,fracKst1680(copy.fracKst1680)
	,phaseA0Kst1680(copy.phaseA0Kst1680)
	,phaseApKst1680(copy.phaseApKst1680)
	,phaseAmKst1680(copy.phaseAmKst1680)
	
	,fracK01430(copy.fracK01430)
	,phaseA0K01430(copy.phaseA0K01430)
	,fracK21430(copy.fracK21430)
	,phaseA0K21430(copy.phaseA0K21430)
	
	,massZplus(copy.massZplus)
	,widthZplus(copy.widthZplus)
	,massKst892(copy.massKst892)
	,widthKst892(copy.widthKst892)
	,massKst1410(copy.massKst1410)
	,widthKst1410(copy.widthKst1410)
	,massKst1680(copy.massKst1680)
	,widthKst1680(copy.widthKst1680)
	,massK01430(copy.massK01430)
	,widthK01430(copy.widthK01430)
	,massK21430(copy.massK21430)
	,widthK21430(copy.widthK21430)
	,pMuPlus(copy.pMuPlus)
	,pMuMinus(copy.pMuMinus)
	,pPi(copy.pPi)
	,pK(copy.pK)
	,KpiComponents()//copy.KpiComponents)
	,ZComponents()//copy.ZComponents)
	,wigner(copy.wigner)
	,mJpsi(copy.mJpsi)
	,useAngularAcceptance(copy.useAngularAcceptance)
	,pB(copy.pB)
	,cosARefs(copy.cosARefs)
	,histogramFile(copy.histogramFile)
	,angularAccHistCosTheta1(copy.angularAccHistCosTheta1)
	,angularAccHistPhi(copy.angularAccHistPhi)
	,angularAccHistMassCosTheta2(copy.angularAccHistMassCosTheta2)
{
	this->SetNumericalNormalisation(true);
	this->TurnCachingOff();
	componentIndex = 0;
	vector<DPComponent*>::const_iterator component_i = copy.KpiComponents.begin();
	for( ; component_i != copy.KpiComponents.end(); ++component_i )
	{
		KpiComponents.push_back( new DPJpsiKaon( *( dynamic_cast<DPJpsiKaon*>( *(component_i) ) ) ) );
	}
	vector<DPComponent*>::const_iterator component_j = copy.ZComponents.begin();
	for( ; component_j != copy.ZComponents.end(); ++component_j )
	{
		ZComponents.push_back( new DPZplusK( *( dynamic_cast<DPZplusK*>( *(component_j) ) ) ) );
	}
}


//Make the data point and parameter set
void DPTotalAmplitudePDF::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( m23Name );
	allObservables.push_back( cosTheta1Name );
	allObservables.push_back( cosTheta2Name );
	allObservables.push_back( phiName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( fracA0sqZplusName );
	parameterNames.push_back( fracApsqZplusName );
	parameterNames.push_back( fracZplusName );
	parameterNames.push_back( phaseA0ZplusName );
	parameterNames.push_back( phaseApZplusName );
	parameterNames.push_back( phaseAmZplusName );
	
	parameterNames.push_back( fracA0sqKst892Name );
	parameterNames.push_back( fracApsqKst892Name );
	//parameterNames.push_back( fracKst892Name );
	parameterNames.push_back( phaseA0Kst892Name );
	parameterNames.push_back( phaseApKst892Name );
	parameterNames.push_back( phaseAmKst892Name );
	
	parameterNames.push_back( fracA0sqKst1410Name );
	parameterNames.push_back( fracApsqKst1410Name );
	parameterNames.push_back( fracKst1410Name );
	parameterNames.push_back( phaseA0Kst1410Name );
	parameterNames.push_back( phaseApKst1410Name );
	parameterNames.push_back( phaseAmKst1410Name );
	
	parameterNames.push_back( fracA0sqKst1680Name );
	parameterNames.push_back( fracApsqKst1680Name );
	parameterNames.push_back( fracKst1680Name );
	parameterNames.push_back( phaseA0Kst1680Name );
	parameterNames.push_back( phaseApKst1680Name );
	parameterNames.push_back( phaseAmKst1680Name );
	
	parameterNames.push_back( fracK01430Name );
	parameterNames.push_back( phaseA0K01430Name );
	parameterNames.push_back( fracK21430Name );
	parameterNames.push_back( phaseA0K21430Name );
	
	parameterNames.push_back( massZplusName );
	parameterNames.push_back( widthZplusName );
	parameterNames.push_back( massKst892Name );
	parameterNames.push_back( widthKst892Name );
	parameterNames.push_back( massKst1410Name );
	parameterNames.push_back( widthKst1410Name );
	parameterNames.push_back( massKst1680Name );
	parameterNames.push_back( widthKst1680Name );
	parameterNames.push_back( massK01430Name );
	parameterNames.push_back( widthK01430Name );
	parameterNames.push_back( massK21430Name );
	parameterNames.push_back( widthK21430Name );
	allParameters = ParameterSet(parameterNames);
}

//Destructor
DPTotalAmplitudePDF::~DPTotalAmplitudePDF()
{
	// destroy components
	for (unsigned int i = 0; i < KpiComponents.size(); ++i)
	{
		delete KpiComponents[i];
		KpiComponents[i]=0;
	}
	for (unsigned int i=0;i<ZComponents.size();++i)
	{
		delete ZComponents[i];
		ZComponents[i]=0;
	}
	if ( useAngularAcceptance ) histogramFile->Close();
}

bool DPTotalAmplitudePDF::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	fracA0sqZplus    = allParameters.GetPhysicsParameter( fracA0sqZplusName )->GetValue();
	fracApsqZplus    = allParameters.GetPhysicsParameter( fracApsqZplusName )->GetValue();
	fracZplus    = allParameters.GetPhysicsParameter( fracZplusName )->GetValue();
	phaseA0Zplus   = allParameters.GetPhysicsParameter( phaseA0ZplusName )->GetValue();
	phaseApZplus   = allParameters.GetPhysicsParameter( phaseApZplusName )->GetValue();
	phaseAmZplus   = allParameters.GetPhysicsParameter( phaseAmZplusName )->GetValue();
	
	fracA0sqKst892    = allParameters.GetPhysicsParameter( fracA0sqKst892Name )->GetValue();
	fracApsqKst892    = allParameters.GetPhysicsParameter( fracApsqKst892Name )->GetValue();
	//fracKst892    = allParameters.GetPhysicsParameter( fracKst892Name )->GetValue();
	phaseA0Kst892  = allParameters.GetPhysicsParameter( phaseA0Kst892Name )->GetValue();
	phaseApKst892  = allParameters.GetPhysicsParameter( phaseApKst892Name )->GetValue();
	phaseAmKst892  = allParameters.GetPhysicsParameter( phaseAmKst892Name )->GetValue();
	
	fracA0sqKst1410   = allParameters.GetPhysicsParameter( fracA0sqKst1410Name )->GetValue();
	fracApsqKst1410   = allParameters.GetPhysicsParameter( fracApsqKst1410Name )->GetValue();
	fracKst1410   = allParameters.GetPhysicsParameter( fracKst1410Name )->GetValue();
	phaseA0Kst1410 = allParameters.GetPhysicsParameter( phaseA0Kst1410Name )->GetValue();
	phaseApKst1410 = allParameters.GetPhysicsParameter( phaseApKst1410Name )->GetValue();
	phaseAmKst1410 = allParameters.GetPhysicsParameter( phaseAmKst1410Name )->GetValue();
	
	fracA0sqKst1680   = allParameters.GetPhysicsParameter( fracA0sqKst1680Name )->GetValue();
	fracApsqKst1680   = allParameters.GetPhysicsParameter( fracApsqKst1680Name )->GetValue();
	fracKst1680   = allParameters.GetPhysicsParameter( fracKst1680Name )->GetValue();
	phaseA0Kst1680 = allParameters.GetPhysicsParameter( phaseA0Kst1680Name )->GetValue();
	phaseApKst1680 = allParameters.GetPhysicsParameter( phaseApKst1680Name )->GetValue();
	phaseAmKst1680 = allParameters.GetPhysicsParameter( phaseAmKst1680Name )->GetValue();
	
	fracK01430    = allParameters.GetPhysicsParameter( fracK01430Name )->GetValue();
	phaseA0K01430  = allParameters.GetPhysicsParameter( phaseA0K01430Name )->GetValue();
	fracK21430    = allParameters.GetPhysicsParameter( fracK21430Name )->GetValue();
	phaseA0K21430  = allParameters.GetPhysicsParameter( phaseA0K21430Name )->GetValue();

	// Sum of all amplitudes must equal 1
	fracKst892 = ((1. - fracZplus - fracKst1410 - fracKst1680 - fracK01430 - fracK21430) < 0.) ? 0. : (1. - fracZplus - fracKst1410 - fracKst1680 - fracK01430 - fracK21430);

	//cout << fracKst892 << " " << fracZplus << " " << fracKst1410 << " " << fracKst1680 << " " << fracK01430 << " " << fracK21430 << endl;
	
	double magA0Zplus   = sqrt(fracA0sqZplus*fracZplus);
	double magApZplus   = sqrt(fracApsqZplus*fracZplus);
	double magAmZplus   = ( (fracZplus - magA0Zplus*magA0Zplus - magApZplus*magApZplus) < 0.) ? 0. : sqrt(fracZplus - magA0Zplus*magA0Zplus - magApZplus*magApZplus);
	double magA0Kst892  = sqrt(fracA0sqKst892*fracKst892);
	double magApKst892  = sqrt(fracApsqKst892*fracKst892);
	double magAmKst892  = ( (fracKst892 - magA0Kst892*magA0Kst892 - magApKst892*magApKst892) < 0.) ? 0. : sqrt(fracKst892 - magA0Kst892*magA0Kst892 - magApKst892*magApKst892);
	double magA0Kst1410 = sqrt(fracA0sqKst1410*fracKst1410);
	double magApKst1410 = sqrt(fracApsqKst1410*fracKst1410);
	double magAmKst1410 = ( (fracKst1410 - magA0Kst1410*magA0Kst1410 - magApKst1410*magApKst1410) < 0.) ? 0. : sqrt(fracKst1410 - magA0Kst1410*magA0Kst1410 - magApKst1410*magApKst1410);
	double magA0Kst1680 = sqrt(fracA0sqKst1680*fracKst1680);
	double magApKst1680 = sqrt(fracApsqKst1680*fracKst1680);
	double magAmKst1680 = ( (fracKst1680 - magA0Kst1680*magA0Kst1680 - magApKst1680*magApKst1680) < 0.) ? 0. : sqrt(fracKst1680 - magA0Kst1680*magA0Kst1680 - magApKst1680*magApKst1680);

	//cout << fracKst892 << " " << magA0Kst892 << " " << magApKst892 << " " << magAmKst892 << " " << (1.- (fracKst892 - magA0Kst892*magA0Kst892 - magApKst892*magApKst892)) << endl;
	
	massZplus  = allParameters.GetPhysicsParameter( massZplusName )->GetValue();
	widthZplus = allParameters.GetPhysicsParameter( widthZplusName )->GetValue();
	massKst892  = allParameters.GetPhysicsParameter( massKst892Name )->GetValue();
	widthKst892 = allParameters.GetPhysicsParameter( widthKst892Name )->GetValue();
	massKst1410  = allParameters.GetPhysicsParameter( massKst1410Name )->GetValue();
	widthKst1410 = allParameters.GetPhysicsParameter( widthKst1410Name )->GetValue();
	massKst1680  = allParameters.GetPhysicsParameter( massKst1680Name )->GetValue();
	widthKst1680 = allParameters.GetPhysicsParameter( widthKst1680Name )->GetValue();
	massK01430  = allParameters.GetPhysicsParameter( massK01430Name )->GetValue();
	widthK01430 = allParameters.GetPhysicsParameter( widthK01430Name )->GetValue();
	massK21430  = allParameters.GetPhysicsParameter( massK21430Name )->GetValue();
	widthK21430 = allParameters.GetPhysicsParameter( widthK21430Name )->GetValue();

	// No checks performed here to ensure that parameters are set correctly		
	//ZComponents[0]  ->setResonanceParameters( massZplus, widthZplus );
	KpiComponents[0]->setResonanceParameters( massKst892, widthKst892 );
	KpiComponents[1]->setResonanceParameters( massKst1410, widthKst1410 );
	KpiComponents[2]->setResonanceParameters( massKst1680, widthKst1680 );
	KpiComponents[3]->setResonanceParameters( massK01430, widthK01430 );
	KpiComponents[4]->setResonanceParameters( massK21430, widthK21430 );
	//ZComponents[0]  ->setHelicityAmplitudes(magA0Zplus, magApZplus, magAmZplus, phaseA0Zplus, phaseApZplus, phaseAmZplus);
	KpiComponents[0]->setHelicityAmplitudes(magA0Kst892,  magApKst892, magAmKst892, phaseA0Kst892, phaseApKst892, phaseAmKst892);
	KpiComponents[1]->setHelicityAmplitudes(magA0Kst1410, magApKst1410, magAmKst1410, phaseA0Kst1410, phaseApKst1410, phaseAmKst1410);
	KpiComponents[2]->setHelicityAmplitudes(magA0Kst1680, magApKst1680, magAmKst1680, phaseA0Kst1680, phaseApKst1680, phaseAmKst1680);
	KpiComponents[3]->setHelicityAmplitudes(sqrt(fracK01430), 0., 0., phaseA0K01430, 0., 0.);
	KpiComponents[4]->setHelicityAmplitudes(sqrt(fracK21430), 0., 0., phaseA0K21430, 0., 0.);
	
	return isOK;
}

//Calculate the function value
double DPTotalAmplitudePDF::Evaluate(DataPoint * measurement)
{
	// Observables
	m23       = measurement->GetObservable( m23Name )->GetValue();
	cosTheta1 = measurement->GetObservable( cosTheta1Name )->GetValue();
	cosTheta2 = measurement->GetObservable( cosTheta2Name )->GetValue();
	phi       = measurement->GetObservable( phiName )->GetValue();

	double angularAccCosTheta1 = 1.;
	double angularAccPhi = 1.;
	double angularAccMassCosTheta2 = 1.;
	if ( useAngularAcceptance )
	{
		angularAccCosTheta1     = angularAccHistCosTheta1	->GetBinContent( angularAccHistCosTheta1	->FindBin(cosTheta1) );
		angularAccPhi           = angularAccHistPhi		->GetBinContent( angularAccHistPhi		->FindBin(phi) );
		angularAccMassCosTheta2 = angularAccHistMassCosTheta2	->GetBinContent( angularAccHistMassCosTheta2	->FindBin(m23, cosTheta2) );
	}
	/*
	// Need angle between reference axis
	DPHelpers::calculateFinalStateMomenta(5.279, m23, mJpsi,
	cosTheta1,  cosTheta2, phi, 0.105, 0.105, 0.13957018, 0.493677,
	pMuPlus, pMuMinus, pPi, pK);
	// Cos of the angle between psi reference axis
	cosARefs = DPHelpers::referenceAxisCosAngle(pB, pMuPlus, pMuMinus, pPi, pK);
	double cosThetaZ;
	double cosThetaPsi;
	double dphi;
	DPHelpers::calculateZplusAngles(pB, pMuPlus, pMuMinus, pPi, pK,
	&cosThetaZ, &cosThetaPsi, &dphi);
	double m13 = (pMuPlus + pMuMinus + pPi).M();

	cout << cosARefs << " " << m13 << " " << cosThetaZ << " " << cosThetaPsi << " " << dphi << endl;
	 */
	double result = 0.;
	TComplex tmp(0,0);

	unsigned int lower = componentIndex - 1;
	unsigned int upper = componentIndex;

	//cout << "componentIndex " << componentIndex << endl;

	if ( componentIndex == 0 ) {
		lower = 0;
		upper = KpiComponents.size();
	}
	// Now sum over final state helicities (this is not general code, but
	// knows about internals of components
	for (int twoLambda = -2; twoLambda <= 2; twoLambda += 4) // Sum over +-1
	{
		tmp = TComplex(0,0);
		for (int twoLambdaPsi = -2; twoLambdaPsi <= 2; twoLambdaPsi += 2) // Sum over -1,0,+1
		{
			for (unsigned int i = lower; i < upper; ++i) // sum over all components
			{
				tmp += KpiComponents[i]->amplitude(m23, cosTheta1, cosTheta2, phi,
						twoLambda, twoLambdaPsi);
				//cout << "m23: " << m23 << " " << cosTheta1 << " " << cosTheta2 << " " << phi << " " << tmp.Re() << " " << tmp.Im() << endl;
			}

			/*
			// Now comes sum over Z+ components and lambdaPsiPrime
			for (unsigned int i = 0; i < ZComponents.size(); ++i)
			{
			// Sum over lambdaPsiPrime
			for (int twoLambdaPsiPrime = -2; twoLambdaPsiPrime <= 2; twoLambdaPsiPrime += 2)
			{
			tmp += wigner.function(cosARefs, twoLambdaPsiPrime/2, twoLambdaPsi/2)*
			ZComponents[i]->amplitude(m13, cosThetaZ, cosThetaPsi, dphi,
			twoLambda,twoLambdaPsiPrime);
			}
			}
			 */
		}
		result += tmp.Rho2();
	}
	//cout << angularAccCosTheta1*angularAccPhi*angularAccMassCosTheta2 << endl;
	return result * angularAccCosTheta1*angularAccPhi*angularAccMassCosTheta2;
}

vector<string> DPTotalAmplitudePDF::PDFComponents()
{
        vector<string> component_list;
        component_list.push_back( "892" );
        //component_list.push_back( "1410" );
        //component_list.push_back( "1680" );
        //component_list.push_back( "1430" );
        //component_list.push_back( "1430_2" );
        //component_list.push_back( "LASS" );
        component_list.push_back( "0" );
        return component_list;
}

//Calculate the function value
double DPTotalAmplitudePDF::EvaluateComponent(DataPoint * measurement, ComponentRef* Component)
{
        componentIndex = Component->getComponentNumber();
        if( componentIndex == -1 )
        {
                string ComponentName = Component->getComponentName();
                if( ComponentName.compare( "892" ) == 0 )
                {
                        Component->setComponentNumber( 1 );
                        componentIndex = 1;
                }
                else if( ComponentName.compare( "1410" ) == 0 )
                {
                        Component->setComponentNumber( 2 );
                        componentIndex = 2;
                }
                else if( ComponentName.compare( "1680" ) == 0 )
                {
                        Component->setComponentNumber( 3 );
                        componentIndex = 3;
                }
                else if( ComponentName.compare( "1430" ) == 0 )
                {
                        Component->setComponentNumber( 4 );
                        componentIndex = 4;
                }
                else if( ComponentName.compare( "1430_2" ) == 0 )
                {
                        Component->setComponentNumber( 5 );
                        componentIndex = 5;
                }
                else if( ComponentName.compare( "LASS" ) == 0 )
                {
                        Component->setComponentNumber( 6 );
                        componentIndex = 6;
                }
                else
                {
                        Component->setComponentNumber( 0 );
                        componentIndex = 0;
                }
        }

        double return_value = this->Evaluate( measurement );

        return return_value;
}
