
#include <string>
#include <iostream>

#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TMath.h"

#include "ROOT_File_Processing.h"
#include "TTree_Processing.h"
#include "Histo_Processing.h"
#include "StringOperations.h"

using namespace::std;

void AngDist( int argc, char* argv[] );

int main( int argc, char* argv[] )
{
	if( argc != 2 )
	{
		cout << "Usage: " << argv[0] << " someInputFile.root"<< endl;
		exit(0);
	}
	AngDist(argc,argv);
	return 0;
}

void AngDist( int argc, char* argv[] )
{
	TString fileName( argv[1] );

	TTree* inputTuple = ROOT_File_Processing::GetFirstTree( fileName.Data() );

	TString lowerSideBand="5205<mass&&mass<5325";
	TString upperSideBand="5400<mass&&mass<5520";

	TString MassCut="("+lowerSideBand+")||("+upperSideBand+")";

	//TString CutString=TString("((selA==1)||(selB==1))*")+MassCut;

	//TString CutString="sWeight*(time>0.3&&time<20.)";

	TString CutString="(time>0.3&&time<14.)&&"+MassCut;

	int binning_costheta = 7;	double costheta_hi=1., costheta_lo=-1.;
	int binning_cospsi = 5;		double cospsi_hi=1., cospsi_lo=-1.;
	int binning_phi = 9;		double phi_hi=TMath::Pi(), phi_lo=-TMath::Pi();

	int binning_costhetaK = 16;	double costhetaK_lo=-1., costhetaK_hi=1.;
	int binning_costhetaL = 24;	double costhetaL_lo=-1., costhetaL_hi=1.;
	int binning_phih = 5;		double helphi_lo=-TMath::Pi(), helphi_hi=TMath::Pi();

	TString cosPsiName="trcospsi";
	TString cosThetaName="trcostheta";
	TString phiName="trphi";

	TString helcosthetaKName="helcosthetaK";
	TString helcosthetaLName="helcosthetaL";
	TString helphiName="helphi";


	TString helphiTransform="("+helphiName+"+TMath::Pi()-("+helphiName+">0.)*(2.*TMath::Pi()))";

	//	1D Projections Transversity
	TH1D* psihisto     = new TH1D(cosPsiName+"_proj", cosPsiName, binning_costheta, cospsi_lo, cospsi_hi);
	TH1D* thetahisto   = new TH1D(cosThetaName+"_proj", cosThetaName, binning_cospsi, costheta_lo, costheta_hi);
	TH1D* phihisto     = new TH1D(phiName+"_proj", phiName, binning_phi, phi_lo, phi_hi);
	//	1D Projections Helicity
	TH1D* helkhisto    = new TH1D(helcosthetaKName+"_proj", helcosthetaKName, binning_costhetaK, costhetaK_lo, costhetaK_hi);
	TH1D* hellhisto    = new TH1D(helcosthetaLName+"_proj", helcosthetaLName, binning_costhetaL, costhetaL_lo, costhetaL_hi);
	TH1D* helphihisto  = new TH1D(helphiName+"_proj", helphiName, binning_phih, helphi_lo, helphi_hi);

	//      2D Projections Transversity
	TH2D* psiTheta     = new TH2D(cosPsiName+"_"+cosThetaName,cosPsiName+":"+cosThetaName, binning_cospsi, cospsi_lo, cospsi_hi, binning_costheta, costheta_lo,costheta_hi);
	TH2D* psiPhi       = new TH2D(cosPsiName+"_"+phiName,cosPsiName+":"+phiName, binning_cospsi, cospsi_lo,cospsi_hi, binning_phi, phi_lo, phi_hi);
	TH2D* thetaPhi     = new TH2D(cosThetaName+"_"+phiName,cosThetaName+":"+phiName, binning_costheta, costheta_lo,costheta_hi, binning_phi, phi_lo, phi_hi);
	//      2D Projections Helicity
	TH2D* thetaKthataL = new TH2D(helcosthetaKName+"_"+helcosthetaLName,helcosthetaKName+":"+helcosthetaLName, binning_costhetaK, costhetaK_lo,costhetaK_hi, binning_costhetaL, costhetaL_lo,costhetaL_hi );
	TH2D* thetaKphi    = new TH2D(helcosthetaKName+"_"+helphiName,helcosthetaKName+":"+helphiName, binning_costhetaK, costhetaK_lo,costhetaK_hi, binning_phih, helphi_lo,helphi_hi );
	TH2D* thetaLphi    = new TH2D(helcosthetaLName+"_"+helphiName,helcosthetaLName+":"+helphiName, binning_costhetaL, costhetaL_lo,costhetaL_hi, binning_phih, helphi_lo,helphi_hi );

	//      3D Projections Transversity
	TH3D* threeD  = new TH3D("histo", "histo", binning_cospsi,cospsi_lo,cospsi_hi, binning_costheta, costheta_lo, costheta_hi, binning_phi, phi_lo, phi_hi);
	TH3D* threeDH = new TH3D("histoHel", "histoHel", binning_costhetaK,costhetaK_lo,costhetaK_hi, binning_costhetaL,costhetaL_lo,costhetaL_hi, binning_phih, helphi_lo,helphi_hi);


	vector<double>* costheta_data = TTree_Processing::Buffer_Branch( inputTuple, cosThetaName, CutString );
	vector<double>* cospsi_data = TTree_Processing::Buffer_Branch( inputTuple, cosPsiName, CutString );
	vector<double>* phi_data = TTree_Processing::Buffer_Branch( inputTuple, phiName, CutString );

	vector<double>* costhetaK_data = TTree_Processing::Buffer_Branch( inputTuple, helcosthetaKName, CutString );
	vector<double>* costhetaL_data = TTree_Processing::Buffer_Branch( inputTuple, helcosthetaLName, CutString );
	vector<double>* helphi_data = TTree_Processing::Buffer_Branch( inputTuple, helphiTransform, CutString );

	for( unsigned int i=0; i< costheta_data->size(); ++i )
	{
		thetahisto->Fill( (*costheta_data)[i] );
		psihisto->Fill( (*cospsi_data)[i] );
		phihisto->Fill( (*phi_data)[i] );
		helkhisto->Fill( (*costhetaK_data)[i] );
		hellhisto->Fill( (*costhetaL_data)[i] );
		helphihisto->Fill( (*helphi_data)[i] );

		thetaPhi->Fill( (*costheta_data)[i], (*phi_data)[i] );
		psiTheta->Fill( (*cospsi_data)[i], (*costheta_data)[i] );
		psiPhi->Fill( (*cospsi_data)[i], (*phi_data)[i] );
		thetaKthataL->Fill( (*costhetaK_data)[i], (*costhetaL_data)[i] );
		thetaKphi->Fill( (*costhetaK_data)[i], (*helphi_data)[i] );
		thetaLphi->Fill( (*costhetaL_data)[i], (*helphi_data)[i] );

		threeD->Fill( (*cospsi_data)[i], (*costheta_data)[i], (*phi_data)[i] );
		threeDH->Fill( (*costhetaK_data)[i], (*costhetaL_data)[i], (*helphi_data)[i] );
	}

	TString AngDist="AngularDist_"; AngDist.Append(StringOperations::TimeString());

	TFile* outFile = new TFile(AngDist+".root","RECREATE");

	threeD->Write( "", TObject::kOverwrite );
	threeDH->Write( "", TObject::kOverwrite );
	psihisto->Write( "", TObject::kOverwrite );
	thetahisto->Write( "", TObject::kOverwrite );
	phihisto->Write( "", TObject::kOverwrite );
	helkhisto->Write( "", TObject::kOverwrite );
	hellhisto->Write( "", TObject::kOverwrite );
	helphihisto->Write( "", TObject::kOverwrite );
	psiTheta->Write( "", TObject::kOverwrite );
	psiPhi->Write( "", TObject::kOverwrite );
	thetaPhi->Write( "", TObject::kOverwrite );
	thetaKthataL->Write( "", TObject::kOverwrite );
	thetaKphi->Write( "", TObject::kOverwrite );
	thetaLphi->Write( "", TObject::kOverwrite );

	outFile->Close();

}

