
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
		cout << "Usage: " << argv[0] << " nTuple2Weight.root  weights.root" << endl;
		exit(0);
	}

	TFile* output = new TFile( argv[2], "READ");

	TH1D* comp_helcosthetaL = (TH1D*) gDirectory->Get("comp_costhetaL");
	TH1D* comp_helcosthetaK = (TH1D*) gDirectory->Get("comp_costhetaK");
	TH1D* comp_helphi = (TH1D*) gDirectory->Get("comp_helphi");

	TFile* MC = new TFile( argv[1], "UPDATE" );
	(void) MC;

	TTree* MasterTree = ROOT_File_Processing::GetFirstTree( argv[1] );

	string masterCut="";//"selA&&(5200<=mass&&mass<=5550)&&(0.3<=time&&time<=14.)&&sigmat<=0.12&&(990<=mdau2&&mdau2<=1050)";

	vector<double>* timeData      = TTree_Processing::Buffer_Branch( MasterTree, "time", masterCut );
	vector<double>* sigmatData    = TTree_Processing::Buffer_Branch( MasterTree, "sigmat", masterCut );
	vector<double>* costhetaKData = TTree_Processing::Buffer_Branch( MasterTree, "helcosthetaK", masterCut );
	vector<double>* costhetaLData = TTree_Processing::Buffer_Branch( MasterTree, "helcosthetaL", masterCut );
	vector<double>* helphiData    = TTree_Processing::Buffer_Branch( MasterTree, "helphi", masterCut );
	vector<double>* costhetaData = TTree_Processing::Buffer_Branch( MasterTree, "trcostheta", masterCut );
	vector<double>* cospsiData = TTree_Processing::Buffer_Branch( MasterTree, "trcospsi", masterCut );
	vector<double>* trphiData    = TTree_Processing::Buffer_Branch( MasterTree, "trphi", masterCut );

	vector<double> weighted_time, weighted_sigmat, weighted_costhetaL, weighted_costhetaK, weighted_helphi;
	vector<double> weighted_trphi, weighted_costheta, weighted_cospsi;
	vector<double> total_weights;

	for( unsigned int i=0; i< costhetaKData->size(); ++i )
	{
		double weight_costhetaL = comp_helcosthetaL->GetBinContent( comp_helcosthetaL->FindBin( (*costhetaLData)[i] ) );
		double weight_costhetaK = comp_helcosthetaK->GetBinContent( comp_helcosthetaK->FindBin( (*costhetaKData)[i] ) );
		double weight_phi = comp_helphi->GetBinContent( comp_helphi->FindBin( (*helphiData)[i] ) );
		double total_weight = weight_costhetaL*weight_costhetaK*weight_phi;

		weighted_time.push_back( (*timeData)[i]*total_weight );
		weighted_sigmat.push_back( (*sigmatData)[i]*total_weight );
		weighted_costhetaL.push_back( (*costhetaLData)[i]*total_weight );
		weighted_costhetaK.push_back( (*costhetaKData)[i]*total_weight );
		weighted_helphi.push_back( (*helphiData)[i]*total_weight );
		weighted_costheta.push_back( (*costhetaData)[i]*total_weight );
		weighted_cospsi.push_back( (*cospsiData)[i]*total_weight );
		weighted_trphi.push_back( (*trphiData)[i]*total_weight );

		total_weights.push_back(total_weight);
	}

	TTree_Processing::AddBranch( MasterTree, "time_reWeighted", weighted_time );
	TTree_Processing::AddBranch( MasterTree, "sigmat_reWeighted", weighted_sigmat );
	TTree_Processing::AddBranch( MasterTree, "helcosthetaL_reWeighted", weighted_costhetaL );
	TTree_Processing::AddBranch( MasterTree, "helcosthetaK_reWeighted", weighted_costhetaK );
	TTree_Processing::AddBranch( MasterTree, "helphi_reWeighted", weighted_helphi );
	TTree_Processing::AddBranch( MasterTree, "trcostheta_reWeighted", weighted_costheta );
	TTree_Processing::AddBranch( MasterTree, "trcospsi_reWeighted", weighted_cospsi );
	TTree_Processing::AddBranch( MasterTree, "trphi_reWeighted", weighted_trphi );
	TTree_Processing::AddBranch( MasterTree, "SWave_Weight", total_weights );

	MasterTree->Write("",TObject::kOverwrite);

	return 0;
}

