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
	  magA0ZplusName		( configurator->getName("magA0Zplus") )
	, magApZplusName		( configurator->getName("magApZplus") )
	, magAmZplusName		( configurator->getName("magAmZplus") )
	, phaseA0ZplusName	( configurator->getName("phaseA0Zplus") )
	, phaseApZplusName	( configurator->getName("phaseApZplus") )
	, phaseAmZplusName	( configurator->getName("phaseAmZplus") )
	
	, magA0Kst892Name	( configurator->getName("magA0Kst892") )
	, magApKst892Name	( configurator->getName("magApKst892") )
	, magAmKst892Name	( configurator->getName("magAmKst892") )
	, phaseA0Kst892Name	( configurator->getName("phaseA0Kst892") )
	, phaseApKst892Name	( configurator->getName("phaseApKst892") )
	, phaseAmKst892Name	( configurator->getName("phaseAmKst892") )
	
	, magA0Kst1410Name	( configurator->getName("magA0Kst1410") )
	, magApKst1410Name	( configurator->getName("magApKst1410") )
	, magAmKst1410Name	( configurator->getName("magAmKst1410") )
	, phaseA0Kst1410Name	( configurator->getName("phaseA0Kst1410") )
	, phaseApKst1410Name	( configurator->getName("phaseApKst1410") )
	, phaseAmKst1410Name	( configurator->getName("phaseAmKst1410") )
	
	, magA0Kst1680Name	( configurator->getName("magA0Kst1680") )
	, magApKst1680Name	( configurator->getName("magApKst1680") )
	, magAmKst1680Name	( configurator->getName("magAmKst1680") )
	, phaseA0Kst1680Name	( configurator->getName("phaseA0Kst1680") )
	, phaseApKst1680Name	( configurator->getName("phaseApKst1680") )
	, phaseAmKst1680Name	( configurator->getName("phaseAmKst1680") )
	
	, magA0K01430Name	( configurator->getName("magA0K01430") )
	, phaseA0K01430Name	( configurator->getName("phaseA0K01430") )
	
	, magA0K21430Name	( configurator->getName("magA0K21430") )
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

  	, mag_LASSName	( configurator->getName("mag_LASS") )
  	, phase_LASSName	( configurator->getName("phase_LASS") )
	, a_LASSName	( configurator->getName("a_LASS") )
	, r_LASSName	( configurator->getName("r_LASS") )
	// The actual values of the parameters and observables
	, magA0Zplus(),  magApZplus(),  magAmZplus(),  phaseA0Zplus(),  phaseApZplus(),  phaseAmZplus()
	, magA0Kst892(),  magApKst892(), magAmKst892(), phaseA0Kst892(),  phaseApKst892(),  phaseAmKst892()
	, magA0Kst1410(), magApKst1410(), magAmKst1410(), phaseA0Kst1410(), phaseApKst1410(), phaseAmKst1410()
	, magA0Kst1680(), magApKst1680(), magAmKst1680(), phaseA0Kst1680(), phaseApKst1680(), phaseAmKst1680()
	, magA0K01430(),  				  phaseA0K01430()
	, magA0K21430(),  				  phaseA0K21430()
	, massZplus(), widthZplus()
	, massKst892(), widthKst892()
	, massKst1410(), widthKst1410()
	, massKst1680(), widthKst1680()
	, massK01430(), widthK01430()
	, massK21430(), widthK21430()
	, m23(), cosTheta1(), cosTheta2(), phi()

	//LASS parameters
	  , a_LASS(), r_LASS(), mag_LASS(), phase_LASS() 

	, pMuPlus(0., 0., 0., 0.), pMuMinus(0., 0., 0., 0.), pPi(0., 0., 0., 0.), pK(0., 0., 0., 0.), pB(0., 0., 0., 5.279)
	, cosARefs()
	, useFourDHistogram(false)
	, fullFileName(), histogramFile(), histo(), angularAccHistCosTheta1(), angularAccHistPhi(), angularAccHistMassCosTheta2()
        , xaxis(), yaxis(), zaxis(), maxis()
        , nxbins(), nybins(), nzbins(), nmbins()
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

  	tmp=new DPJpsiKaon(0, 0, 5.279, 1.435, 0.279, 0.493677,
                     0.13957018, 5.0, 1.5, 3.096916,0,
                     "LASS", 1.94, 1.76);

	KpiComponents.push_back(tmp);



	// B0 --> Z+ K-
	tmp=new DPZplusK(0,0,5.279,4.430,0.100,0.493677,
			0.13957018, 5.0, 1.5, 3.096916,0);
	//ZComponents.push_back(tmp);

	this->SetNumericalNormalisation( true );
	this->TurnCachingOff();
	useAngularAcceptance = false;
	string fileName = configurator->getConfigurationValue( "AngularAcceptanceHistogram" ) ;
	if ( configurator->isTrue( "UseAngularAcceptance" ) )
	{
		useAngularAcceptance = true;
		//File location
                ifstream input_file;
                input_file.open( fileName.c_str(), ifstream::in );
                input_file.close();

                bool local_fail = input_file.fail();

                if( !getenv("RAPIDFITROOT") && local_fail )
                {
                        cerr << "\n" << endl;
                        //cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                        cerr << "$RAPIDFITROOT NOT DEFINED, PLEASE DEFINE IT SO I CAN USE ACCEPTANCE DATA" << endl;
                        //cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                        cerr << "\n" << endl;
                        exit(-987);
                }

                if( getenv("RAPIDFITROOT") )
                {
                        string path( getenv("RAPIDFITROOT") ) ;

                        cout << "RAPIDFITROOT defined as: " << path << endl;

                        fullFileName = path+"/pdfs/configdata/"+fileName ;
                        input_file.open( fullFileName.c_str(), ifstream::in );
                        input_file.close();
                }
                bool elsewhere_fail = input_file.fail();

                if( elsewhere_fail && local_fail )
                {
                        cerr << "\n\tFileName:\t" << fullFileName << "\t NOT FOUND PLEASE CHECK YOUR RAPIDFITROOT" << endl;
                        cerr << "\t\tAlternativley make sure your XML points to the correct file, or that the file is in the current working directory\n" << endl;
                        exit(-89);
                }

                if( fullFileName.empty() || !local_fail )
                {
                        fullFileName = fileName;
                }

		histogramFile = TFile::Open(fullFileName.c_str());
                
		if (useFourDHistogram) {
		histo = (THnSparse*)histogramFile->Get("histo_4var_eff"); //(fileName.c_str())));
 
                // cos mu
                xaxis = histo->GetAxis(0);
                nxbins = xaxis->GetNbins();

                // cos k
                yaxis = histo->GetAxis(1);
                nybins = yaxis->GetNbins();

                // delta phi
                zaxis = histo->GetAxis(2);
                nzbins = zaxis->GetNbins();

                // m Kpi
                maxis = histo->GetAxis(3);
                nmbins = maxis->GetNbins();

                int total_num_bins = nxbins * nybins * nzbins * nmbins;
		int sum = 0;
		
		vector<int> zero_bins;
                //loop over each bin in histogram and print out how many zero bins there are
                int idx[4] = {0,0,0,0};
                for (int i=1; i < nxbins+1; ++i)
                {
                        for (int j=1; j < nybins+1; ++j)
                        {
                                for (int k=1; k < nzbins+1; ++k)
                                {
                                        for (int l=1; l < nmbins+1; ++l)
                                        {
                                                idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
                                                double bin_content = 0.;//histo->GetBinContent(idx); // why does this not work anymore?
                                                //cout << "Bin content: " << bin_content << endl;
                                                if(bin_content<=0)
                                                {
                                                        zero_bins.push_back(1);
                                                }
                                                //cout << " Zero bins " << zero_bins.size() << endl;}
                                                else if (bin_content>0)
                                                {
                                                        sum += (int) bin_content;
                                                }
                                        }
                                }
                        }
                }

                int average_bin_content = sum / total_num_bins;

                cout << "\n\n\t\t" << "****" << "For total number of bins " << total_num_bins << " there are " << zero_bins.size() << " bins containing zero events " << "****" << endl;
                cout <<  "\t\t\t" << "***" << "Average number of entries of non-zero bins: " << average_bin_content << "***" << endl;
                cout << endl;
		}

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
	,magA0ZplusName(copy.magA0ZplusName)
	,magApZplusName(copy.magApZplusName)
	,magAmZplusName(copy.magAmZplusName)
	,phaseA0ZplusName(copy.phaseA0ZplusName)
	,phaseApZplusName(copy.phaseApZplusName)
	,phaseAmZplusName(copy.phaseAmZplusName)
	
	,magA0Kst892Name(copy.magA0Kst892Name)
	,magApKst892Name(copy.magApKst892Name)
	,magAmKst892Name(copy.magAmKst892Name)
	,phaseA0Kst892Name(copy.phaseA0Kst892Name)
	,phaseApKst892Name(copy.phaseApKst892Name)
	,phaseAmKst892Name(copy.phaseAmKst892Name)
	
	,magA0Kst1410Name(copy.magA0Kst1410Name)
	,magApKst1410Name(copy.magApKst1410Name)
	,magAmKst1410Name(copy.magAmKst1410Name)
	,phaseA0Kst1410Name(copy.phaseA0Kst1410Name)
	,phaseApKst1410Name(copy.phaseApKst1410Name)
	,phaseAmKst1410Name(copy.phaseAmKst1410Name)
	
	,magA0Kst1680Name(copy.magA0Kst1680Name)
	,magApKst1680Name(copy.magApKst1680Name)
	,magAmKst1680Name(copy.magAmKst1680Name)
	,phaseA0Kst1680Name(copy.phaseA0Kst1680Name)
	,phaseApKst1680Name(copy.phaseApKst1680Name)
	,phaseAmKst1680Name(copy.phaseAmKst1680Name)
	
	,magA0K01430Name(copy.magA0K01430Name)
	,phaseA0K01430Name(copy.phaseA0K01430Name)
	,magA0K21430Name(copy.magA0K21430Name)
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


  	,mag_LASSName(copy.mag_LASSName)
  	,phase_LASSName(copy.phase_LASSName)
	,a_LASSName(copy.a_LASSName)
	,r_LASSName(copy.r_LASSName)
	
	,magA0Zplus(copy.magA0Zplus)
	,magApZplus(copy.magApZplus)
	,magAmZplus(copy.magAmZplus)
	,phaseA0Zplus(copy.phaseA0Zplus)
	,phaseApZplus(copy.phaseApZplus)
	,phaseAmZplus(copy.phaseAmZplus)
	
	,magA0Kst892(copy.magA0Kst892)
	,magApKst892(copy.magApKst892)
	,magAmKst892(copy.magAmKst892)
	,phaseA0Kst892(copy.phaseA0Kst892)
	,phaseApKst892(copy.phaseApKst892)
	,phaseAmKst892(copy.phaseAmKst892)
	
	,magA0Kst1410(copy.magA0Kst1410)
	,magApKst1410(copy.magApKst1410)
	,magAmKst1410(copy.magAmKst1410)
	,phaseA0Kst1410(copy.phaseA0Kst1410)
	,phaseApKst1410(copy.phaseApKst1410)
	,phaseAmKst1410(copy.phaseAmKst1410)
	
	,magA0Kst1680(copy.magA0Kst1680)
	,magApKst1680(copy.magApKst1680)
	,magAmKst1680(copy.magAmKst1680)
	,phaseA0Kst1680(copy.phaseA0Kst1680)
	,phaseApKst1680(copy.phaseApKst1680)
	,phaseAmKst1680(copy.phaseAmKst1680)
	
	,magA0K01430(copy.magA0K01430)
	,phaseA0K01430(copy.phaseA0K01430)
	,magA0K21430(copy.magA0K21430)
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
	,useFourDHistogram(copy.useFourDHistogram)
	,fullFileName(copy.fullFileName)
	,histogramFile()//copy.histogramFile)
	,angularAccHistCosTheta1()//copy.angularAccHistCosTheta1)
	,angularAccHistPhi()//copy.angularAccHistPhi)
	,angularAccHistMassCosTheta2()//copy.angularAccHistMassCosTheta2)
	,histo()
	,xaxis(copy.xaxis), yaxis(copy.yaxis), zaxis(copy.zaxis), maxis(copy.maxis)
	,nxbins(copy.nxbins), nybins(copy.nybins), nzbins(copy.nzbins), nmbins(copy.nmbins)
	,mag_LASS(copy.mag_LASS)
	,phase_LASS(copy.phase_LASS)
        ,a_LASS(copy.a_LASS)
        ,r_LASS(copy.r_LASS)
{
	this->SetNumericalNormalisation(true);
	this->TurnCachingOff();
	componentIndex = 0;

	cout << "Making copy of DPTotalAmplitudePDF. Acceptance: " << useAngularAcceptance << endl;

	for( unsigned int i=0; i < copy.KpiComponents.size(); ++i )
	{
		KpiComponents.push_back( new DPJpsiKaon( *((DPJpsiKaon*)copy.KpiComponents[i]) ) );
	}

	for( unsigned int i=0; i < copy.ZComponents.size(); ++i )
	{
		ZComponents.push_back( new DPZplusK( *((DPZplusK*)copy.ZComponents[i]) ) );
	}
        
	if ( useAngularAcceptance )
        {
                histogramFile = TFile::Open(fullFileName.c_str());
                
		if (useFourDHistogram) histo = (THnSparse*)histogramFile->Get("histo_4var_eff"); //(fileName.c_str())));
		else {
	        	angularAccHistCosTheta1 = (TH1D*)histogramFile->Get("cosmu_effTot");
                	angularAccHistPhi       = (TH1D*)histogramFile->Get("delta_phi_effTot");
                	angularAccHistMassCosTheta2 = (TH2D*)histogramFile->Get("mass_cos_effTot");
		}
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
	parameterNames.push_back( magA0ZplusName );
	parameterNames.push_back( magApZplusName );
	parameterNames.push_back( magAmZplusName );
	parameterNames.push_back( phaseA0ZplusName );
	parameterNames.push_back( phaseApZplusName );
	parameterNames.push_back( phaseAmZplusName );
	
	parameterNames.push_back( magA0Kst892Name );
	parameterNames.push_back( magApKst892Name );
	parameterNames.push_back( magAmKst892Name );
	parameterNames.push_back( phaseA0Kst892Name );
	parameterNames.push_back( phaseApKst892Name );
	parameterNames.push_back( phaseAmKst892Name );
	
	parameterNames.push_back( magA0Kst1410Name );
	parameterNames.push_back( magApKst1410Name );
	parameterNames.push_back( magAmKst1410Name );
	parameterNames.push_back( phaseA0Kst1410Name );
	parameterNames.push_back( phaseApKst1410Name );
	parameterNames.push_back( phaseAmKst1410Name );
	
	parameterNames.push_back( magA0Kst1680Name );
	parameterNames.push_back( magApKst1680Name );
	parameterNames.push_back( magAmKst1680Name );
	parameterNames.push_back( phaseA0Kst1680Name );
	parameterNames.push_back( phaseApKst1680Name );
	parameterNames.push_back( phaseAmKst1680Name );
	
	parameterNames.push_back( magA0K01430Name );
	parameterNames.push_back( phaseA0K01430Name );
	parameterNames.push_back( magA0K21430Name );
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

	parameterNames.push_back( mag_LASSName );
	parameterNames.push_back( phase_LASSName );
	parameterNames.push_back( a_LASSName );
	parameterNames.push_back( r_LASSName );


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
	magA0Zplus    = allParameters.GetPhysicsParameter( magA0ZplusName )->GetValue();
	magApZplus    = allParameters.GetPhysicsParameter( magApZplusName )->GetValue();
	magAmZplus    = allParameters.GetPhysicsParameter( magAmZplusName )->GetValue();
	phaseA0Zplus   = allParameters.GetPhysicsParameter( phaseA0ZplusName )->GetValue();
	phaseApZplus   = allParameters.GetPhysicsParameter( phaseApZplusName )->GetValue();
	phaseAmZplus   = allParameters.GetPhysicsParameter( phaseAmZplusName )->GetValue();
	
	magA0Kst892    = allParameters.GetPhysicsParameter( magA0Kst892Name )->GetValue();
	magApKst892    = allParameters.GetPhysicsParameter( magApKst892Name )->GetValue();
	magAmKst892        = allParameters.GetPhysicsParameter( magAmKst892Name )->GetValue();
	phaseA0Kst892  = allParameters.GetPhysicsParameter( phaseA0Kst892Name )->GetValue();
	phaseApKst892  = allParameters.GetPhysicsParameter( phaseApKst892Name )->GetValue();
	phaseAmKst892  = allParameters.GetPhysicsParameter( phaseAmKst892Name )->GetValue();
	
	magA0Kst1410   = allParameters.GetPhysicsParameter( magA0Kst1410Name )->GetValue();
	magApKst1410   = allParameters.GetPhysicsParameter( magApKst1410Name )->GetValue();
	magAmKst1410   = allParameters.GetPhysicsParameter( magAmKst1410Name )->GetValue();
	phaseA0Kst1410 = allParameters.GetPhysicsParameter( phaseA0Kst1410Name )->GetValue();
	phaseApKst1410 = allParameters.GetPhysicsParameter( phaseApKst1410Name )->GetValue();
	phaseAmKst1410 = allParameters.GetPhysicsParameter( phaseAmKst1410Name )->GetValue();
	
	magA0Kst1680   = allParameters.GetPhysicsParameter( magA0Kst1680Name )->GetValue();
	magApKst1680   = allParameters.GetPhysicsParameter( magApKst1680Name )->GetValue();
	magAmKst1680   = allParameters.GetPhysicsParameter( magAmKst1680Name )->GetValue();
	phaseA0Kst1680 = allParameters.GetPhysicsParameter( phaseA0Kst1680Name )->GetValue();
	phaseApKst1680 = allParameters.GetPhysicsParameter( phaseApKst1680Name )->GetValue();
	phaseAmKst1680 = allParameters.GetPhysicsParameter( phaseAmKst1680Name )->GetValue();
	
	magA0K01430    = allParameters.GetPhysicsParameter( magA0K01430Name )->GetValue();
	phaseA0K01430  = allParameters.GetPhysicsParameter( phaseA0K01430Name )->GetValue();
	magA0K21430    = allParameters.GetPhysicsParameter( magA0K21430Name )->GetValue();
	phaseA0K21430  = allParameters.GetPhysicsParameter( phaseA0K21430Name )->GetValue();

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

	mag_LASS = allParameters.GetPhysicsParameter( mag_LASSName )->GetValue();
	phase_LASS = allParameters.GetPhysicsParameter( phase_LASSName )->GetValue();
	a_LASS = allParameters.GetPhysicsParameter( a_LASSName )->GetValue();
	r_LASS = allParameters.GetPhysicsParameter( r_LASSName )->GetValue();


	// No checks performed here to ensure that parameters are set correctly		
	//ZComponents[0]  ->setResonanceParameters( massZplus, widthZplus );
	KpiComponents[0]->setResonanceParameters( massKst892, widthKst892 );
	KpiComponents[1]->setResonanceParameters( massKst1410, widthKst1410 );
	KpiComponents[2]->setResonanceParameters( massKst1680, widthKst1680 );
	KpiComponents[3]->setResonanceParameters( massK01430, widthK01430 );
	KpiComponents[4]->setResonanceParameters( massK21430, widthK21430 );
	KpiComponents[5]->setResonanceParameters( a_LASS, r_LASS );
	//ZComponents[0]  ->setHelicityAmplitudes(magA0Zplus, magApZplus, magAmZplus, phaseA0Zplus, phaseApZplus, phaseAmZplus);
	KpiComponents[0]->setHelicityAmplitudes(magA0Kst892,  magApKst892, magAmKst892, phaseA0Kst892, phaseApKst892, phaseAmKst892);
	KpiComponents[1]->setHelicityAmplitudes(magA0Kst1410, magApKst1410, magAmKst1410, phaseA0Kst1410, phaseApKst1410, phaseAmKst1410);
	KpiComponents[2]->setHelicityAmplitudes(magA0Kst1680, magApKst1680, magAmKst1680, phaseA0Kst1680, phaseApKst1680, phaseAmKst1680);
	KpiComponents[3]->setHelicityAmplitudes(magA0K01430, 0., 0., phaseA0K01430, 0., 0.);
	KpiComponents[4]->setHelicityAmplitudes(magA0K21430, 0., 0., phaseA0K21430, 0., 0.);
	KpiComponents[5]->setHelicityAmplitudes(mag_LASS, 0., 0., phase_LASS, 0., 0.);

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

	int globalbin = -1;
        int xbin = -1, ybin = -1, zbin = -1, mbin = -1;

	double angularAcc = 1.;
	double angularAccCosTheta1 = 1.;
	double angularAccPhi = 1.;
	double angularAccMassCosTheta2 = 1.;
	if ( useAngularAcceptance )
	{
		if ( useFourDHistogram ) {
			//Find global bin number for values of angles, find number of entries per bin, divide by volume per bin and normalise with total number of entries in the histogram
                	xbin = xaxis->FindFixBin( cosTheta1 ); if( xbin > nxbins ) xbin = nxbins;
                	ybin = yaxis->FindFixBin( cosTheta2 ); if( ybin > nybins ) ybin = nybins;
                	zbin = zaxis->FindFixBin( phi  	    ); if( zbin > nzbins ) zbin = nzbins;
                	mbin = maxis->FindFixBin( m23 	    ); if( mbin > nmbins ) mbin = nmbins;

                	int idx[4] = { xbin, ybin, zbin, mbin };
                	globalbin = (int)histo->GetBin( idx );
                	angularAcc = histo->GetBinContent(globalbin);
		}
		else {
			angularAccCosTheta1     = angularAccHistCosTheta1	->GetBinContent( angularAccHistCosTheta1	->FindBin(cosTheta1) );
			angularAccPhi           = angularAccHistPhi		->GetBinContent( angularAccHistPhi		->FindBin(phi) );
			angularAccMassCosTheta2 = angularAccHistMassCosTheta2	->GetBinContent( angularAccHistMassCosTheta2	->FindBin(m23, cosTheta2) );
			angularAcc = angularAccCosTheta1*angularAccPhi*angularAccMassCosTheta2; // factor of 81 = 3^4 to get acceptance on same scale as other quantities
		}
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

	unsigned int lower = (unsigned)(componentIndex - 1);
	unsigned int upper = (unsigned)componentIndex;

	//std::cout << "componentIndex " << componentIndex << std::endl;

	if ( componentIndex == 0 ) {
		lower = 0;
		upper = (unsigned)KpiComponents.size();
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
				//cout << "m23: " << m23 << " " << cosTheta1 << " " << cosTheta2 << " " << phi << " " << tmp.Re() << " " << tmp.Im() << " " << i << endl;
			}

			/*
			// Now comes sum over Z+ components and lambdaPsiPrime
			for (unsigned int i = 0; i < ZComponents.size(); ++i)
			{
			{
			tmp += ZComponents[i]->amplitude(m13, cosThetaZ, cosThetaPsi, dphi,
			twoLambda,twoLambdaPsi); // need to check that we pass right helicities
			}
			}
			 */
		}
		result += tmp.Rho2();
	}
	//cout << angularAccCosTheta1*angularAccPhi*angularAccMassCosTheta2 << endl;

	//momenta are defined on eq 39.20a/b of the 2010 PDG
	const double m1 = 0.493677;    // kaon mass
	const double m2 = 0.13957018; // pion mass
	const double m3 = 3.096916;   // psi mass
	const double MB0= 5.2795; // B0 mass

	double t1 = m23*m23-(m1+m2)*(m1+m2);
	double t2 = m23*m23-(m1-m2)*(m1-m2);
	
	double t31 = MB0*MB0 - (m23 + m3)*(m23 + m3);
	double t32 = MB0*MB0 - (m23 - m3)*(m23 - m3);

	double p1_st = sqrt(t1*t2)/m23/2.;
	double p3    = sqrt(t31*t32)/MB0/2.;

	//std::cout << result << " " << angularAcc << " " << p1_st << " " << p3 << std::endl;

	double returnable_value = result * angularAcc * p1_st * p3;

//  std::cout<<"DEBUG: "<<m23<<" "<<cosTheta1<<" "<<cosTheta2<<" "<<phi<<" "<<result<<" "<<p1_st*p3<<std::endl;

	if( isnan(returnable_value) || returnable_value < 0 ) return 0.;
	else return returnable_value;
}

vector<string> DPTotalAmplitudePDF::PDFComponents()
{
        vector<string> components_list;
        components_list.push_back( "892" );
        components_list.push_back( "1410" );
        components_list.push_back( "1680" );
        components_list.push_back( "1430" );
        components_list.push_back( "1430_2" );
        components_list.push_back( "LASS" );
        components_list.push_back( "0" );
        return components_list;
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

double DPTotalAmplitudePDF::Normalisation(PhaseSpaceBoundary * boundary)
{
        (void) boundary;
	return -1.;
}

