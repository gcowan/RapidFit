
#include "TMath.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"

#include "ROOT_File_Processing.h"
#include "TTree_Processing.h"
#include "Template_Functions.h"

#include <string>
#include <map>
#include <vector>

using namespace::std;


int main( int argc, char* argv[] )
{
	if( argc != 3 )
	{
		cout << "Usage: " << argv[0] << " nTuple_withOut_sWave.root  nTuple_with_sWave.root" << endl;
		exit(0);
	}
	string woTupleName=argv[1];//"Scale_woSWave.root";
	string wTupleName=argv[2];//"Scale_wSWave.root";

	string costhetaLName="helcosthetaL";
	string costhetaKName="helcosthetaK";
	string phiName="helphi";

	TTree* woTree = ROOT_File_Processing::GetFirstTree( woTupleName );

	vector<double>* wocosthetaLData = TTree_Processing::Buffer_Branch( woTree, costhetaLName );
	vector<double>* wocosthetaKData = TTree_Processing::Buffer_Branch( woTree, costhetaKName );
	vector<double>* wophiData = TTree_Processing::Buffer_Branch( woTree, phiName );

	TTree* wTree = ROOT_File_Processing::GetFirstTree( wTupleName );

	vector<double>* wcosthetaLData = TTree_Processing::Buffer_Branch( wTree, costhetaLName );
	vector<double>* wcosthetaKData = TTree_Processing::Buffer_Branch( wTree, costhetaKName );
	vector<double>* wphiData = TTree_Processing::Buffer_Branch( wTree, phiName );

	TFile* output = new TFile("weights.root", "RECREATE");

	TH1D* wo_costhetaL = new TH1D( "wo_costhetaL", "costhetaL" , 100, -1., 1. );
	TH1D* wo_costhetaK = new TH1D( "wo_costhetaK", "costhetaK" , 100, -1., 1. );
	TH1D* wo_phi = new TH1D( "wo_phi", "phi" , 100, -TMath::Pi(), TMath::Pi() );

	TH1D* w_costhetaL = new TH1D( "w_costhetaL", "costhetaL" , 100, -1., 1. );
	TH1D* w_costhetaK = new TH1D( "w_costhetaK", "costhetaK" , 100, -1., 1. );
	TH1D* w_phi = new TH1D( "w_phi", "phi" , 100, -TMath::Pi(), TMath::Pi() );

	for( unsigned int i=0; i< wcosthetaLData->size(); ++i )
	{
		w_costhetaL->Fill( (*wcosthetaLData)[i] );
		w_costhetaK->Fill( (*wcosthetaKData)[i] );
		w_phi->Fill( (*wphiData)[i] );

		wo_costhetaL->Fill( (*wocosthetaLData)[i] );
		wo_costhetaK->Fill( (*wocosthetaKData)[i] );
		wo_phi->Fill( (*wophiData)[i] );
	}

	TH1D* comp_helcosthetaL = new TH1D( "comp_costhetaL", "costhetaL" , 100, -1., 1. );
	TH1D* comp_helcosthetaK = new TH1D( "comp_costhetaK", "costhetaL" , 100, -1., 1. );
	TH1D* comp_helphi = new TH1D( "comp_helphi", "phi", 100, -TMath::Pi(), TMath::Pi() );

	comp_helcosthetaL->Divide( w_costhetaL, wo_costhetaL );
	comp_helcosthetaK->Divide( w_costhetaK, wo_costhetaK );
	comp_helphi->Divide( w_phi, wo_phi );

	comp_helcosthetaL->Write("", TObject::kOverwrite );
	comp_helcosthetaK->Write("", TObject::kOverwrite );
	comp_helphi->Write("", TObject::kOverwrite );

	wo_costhetaL->Write("",TObject::kOverwrite);
	wo_costhetaK->Write("",TObject::kOverwrite);
	wo_phi->Write("",TObject::kOverwrite);

	w_costhetaL->Write("",TObject::kOverwrite);
	w_costhetaK->Write("",TObject::kOverwrite);
	w_phi->Write("",TObject::kOverwrite);

	output->Write("",TObject::kOverwrite);
	output->Close();
	return 0;
}


