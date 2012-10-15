
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

struct eventStruct
{
	double mass;
	double costhetaL;
	double costhetaK;
	double phi;
	double time;
	double sWeight;
};

int main( int argc, char* argv[] )
{
	string MasterTupleName="Bs2JpsiPhi_ntupleB_for_fitting_20120821_MagDownMagUp.root";
	string EdinburghTupleName="1fbUB-6binSig.root";
	string HDTupleName="hd_ubOnly.root";

	string evtNum="100*eventNumber+runNumber";
	string massName="mass";
	string costhetaKName="helcosthetaK";
	string costhetaLName="helcosthetaL";
	string phiName="helphi";
	string timeName="time";

	string masterCut="selA&&(5200<=mass&&mass<=5550)&&(0.3<=time&&time<=14.)&&sigmat<=0.12&&(990<=mdau2&&mdau2<=1050)";
	string EdinburghCut="(5200<=mass&&mass<=5550)&&(0.3<=time&&time<=14.)";
	string HDCut="";

	TTree* MasterTree = ROOT_File_Processing::GetFirstTree( MasterTupleName );

	vector<double>* MasterKey = TTree_Processing::Buffer_Branch( MasterTree, evtNum, masterCut );
	vector<double>* MastermassData = TTree_Processing::Buffer_Branch( MasterTree, massName, masterCut );
	vector<double>* MastercosthetaLData = TTree_Processing::Buffer_Branch( MasterTree, costhetaLName, masterCut );
	vector<double>* MastercosthetaKData = TTree_Processing::Buffer_Branch( MasterTree, costhetaKName, masterCut );
	vector<double>* MasterphiData = TTree_Processing::Buffer_Branch( MasterTree, phiName, masterCut );
	vector<double>* MastertimeData = TTree_Processing::Buffer_Branch( MasterTree, timeName, masterCut );

	TTree* EdinburghTree = ROOT_File_Processing::GetFirstTree( EdinburghTupleName );

	vector<double>* EdiKey = TTree_Processing::Buffer_Branch( EdinburghTree, evtNum, EdinburghCut );
	vector<double>* EdimassData = TTree_Processing::Buffer_Branch( EdinburghTree, massName, EdinburghCut );
	vector<double>* EdicosthetaLData = TTree_Processing::Buffer_Branch( EdinburghTree, costhetaLName, EdinburghCut );
	vector<double>* EdicosthetaKData = TTree_Processing::Buffer_Branch( EdinburghTree, costhetaKName, EdinburghCut );
	vector<double>* EdiphiData = TTree_Processing::Buffer_Branch( EdinburghTree, phiName, EdinburghCut );
	vector<double>* EditimeData = TTree_Processing::Buffer_Branch( EdinburghTree, timeName, EdinburghCut );

	TTree* HDTree = ROOT_File_Processing::GetFirstTree( HDTupleName );

	vector<double>* HDKey = TTree_Processing::Buffer_Branch( HDTree, evtNum, HDCut );
	vector<double>* HDmassData = TTree_Processing::Buffer_Branch( HDTree, massName, HDCut );
	vector<double>* HDcosthetaLData = TTree_Processing::Buffer_Branch( HDTree, costhetaLName, HDCut );
	vector<double>* HDcosthetaKData = TTree_Processing::Buffer_Branch( HDTree, costhetaKName, HDCut );
	vector<double>* HDphiData = TTree_Processing::Buffer_Branch( HDTree, phiName, HDCut );
	vector<double>* HDtimeData = TTree_Processing::Buffer_Branch( HDTree, timeName, HDCut );

	cout << "Num Events:" << endl;
	cout << "Edi:\t" << EdimassData->size() << endl;
	cout << "Master:\t" << MastermassData->size() << endl;
	cout << "HD:\t" << HDmassData->size() << endl;

	map< double, eventStruct > MasterMap;
	map< double, eventStruct > EdiMap;
	map< double, eventStruct > HDMap;

	for( unsigned int i=0; i< MasterKey->size(); ++i )
	{
		struct eventStruct thisEvent;
		thisEvent.mass = (*MastermassData)[i];
		thisEvent.costhetaL = (*MastercosthetaLData)[i];
		thisEvent.costhetaK = (*MastercosthetaKData)[i];
		thisEvent.phi = (*MasterphiData)[i];
		thisEvent.time = (*MastertimeData)[i];
		MasterMap[(*MasterKey)[i]] = thisEvent;
	}

	for( unsigned int i=0; i< EdiKey->size(); ++i )
	{
		struct eventStruct thisEvent;
		thisEvent.mass = (*EdimassData)[i];
		thisEvent.costhetaL = (*EdicosthetaLData)[i];
		thisEvent.costhetaK = (*EdicosthetaKData)[i];
		thisEvent.phi = (*EdiphiData)[i];
		thisEvent.time = (*EditimeData)[i];
		EdiMap[(*EdiKey)[i]] = thisEvent;
	}

	for( unsigned int i=0; i< HDKey->size(); ++i )
	{
		struct eventStruct thisEvent;
		thisEvent.mass = (*HDmassData)[i];
		thisEvent.costhetaL = (*HDcosthetaLData)[i];
		thisEvent.costhetaK = (*HDcosthetaKData)[i];
		thisEvent.phi = (*HDphiData)[i];
		thisEvent.time = (*HDtimeData)[i];
		HDMap[(*HDKey)[i]] = thisEvent;
	}

	cout << "Num Events:" << endl;
	cout << "Edi:\t" << EdiMap.size() << endl;    
	cout << "Master:\t" << MasterMap.size() << endl;  
	cout << "HD:\t" << HDMap.size() << endl;

	vector<double> Edires_costhetaL, Edires_costhetaK, Edires_phi, Edires_time, Edires_mass;
	vector<double> HDres_costhetaL, HDres_costhetaK, HDres_phi, HDres_time, HDres_mass;

	vector<double> sWeightDiff;

	for( map<double, eventStruct>::iterator Master = MasterMap.begin(); Master != MasterMap.end(); ++Master )
	{
		double thisMass = Master->first;

		Edires_costhetaL.push_back( MasterMap[thisMass].costhetaL - EdiMap[thisMass].costhetaL );
		HDres_costhetaL.push_back( MasterMap[thisMass].costhetaK - HDMap[thisMass].costhetaL );

		Edires_costhetaK.push_back( MasterMap[thisMass].costhetaK - EdiMap[thisMass].costhetaK );
		HDres_costhetaK.push_back( MasterMap[thisMass].costhetaL - HDMap[thisMass].costhetaK );

		Edires_phi.push_back( MasterMap[thisMass].phi - EdiMap[thisMass].phi );
		HDres_phi.push_back( MasterMap[thisMass].phi - HDMap[thisMass].phi );

		Edires_time.push_back( MasterMap[thisMass].time - EdiMap[thisMass].time );
		HDres_time.push_back( MasterMap[thisMass].time - HDMap[thisMass].time );

		Edires_mass.push_back( MasterMap[thisMass].mass - EdiMap[thisMass].mass );
		HDres_mass.push_back( MasterMap[thisMass].mass - HDMap[thisMass].mass );

		sWeightDiff.push_back( EdiMap[thisMass].sWeight - HDMap[thisMass].sWeight );
	}


	TFile* residuals = new TFile("Residuals.root","RECREATE");

	TH1D* Edi_costhetaL = new TH1D( "Edires_costhetaL", "costhetaL" , 100, get_minimum(Edires_costhetaL), get_maximum(Edires_costhetaL) );
	TH1D* Edi_costhetaK = new TH1D( "Edires_costhetaK", "costhetaK" , 100, get_minimum(Edires_costhetaK), get_maximum(Edires_costhetaK) );
	TH1D* Edi_phi = new TH1D( "Edires_phi", "phi" , 100, get_minimum(Edires_phi), get_maximum(Edires_phi) );
	TH1D* Edi_mass = new TH1D( "Edires_mass", "mass" , 100, get_minimum(Edires_mass), get_maximum(Edires_mass) );
	TH1D* Edi_time = new TH1D( "Edires_time", "time" , 100, get_minimum(Edires_time), get_maximum(Edires_time) );

	TH1D* HD_costhetaL = new TH1D( "HDres_costhetaL", "costhetaL" , 100, get_minimum(HDres_costhetaL), get_maximum(HDres_costhetaL) );
	TH1D* HD_costhetaK = new TH1D( "HDres_costhetaK", "costhetaK" , 100, get_minimum(HDres_costhetaK), get_maximum(HDres_costhetaK) );
	TH1D* HD_phi = new TH1D( "HDres_phi", "phi" , 100, get_minimum(HDres_phi), get_maximum(HDres_phi) );
	TH1D* HD_mass = new TH1D( "HDres_mass", "mass" , 100, get_minimum(HDres_mass), get_maximum(HDres_mass) );
	TH1D* HD_time = new TH1D( "HDres_time", "time" , 100, get_minimum(HDres_time), get_maximum(HDres_time) );

	TH2D* mass_diff = new TH2D( "sWeightDiff", "diff", 100, get_minimum(*MastermassData), get_maximum(*MastermassData), 100, get_minimum(sWeightDiff), get_maximum(sWeightDiff) );

	for( unsigned int i=0; i< Edires_costhetaL.size(); ++i )
	{
		Edi_costhetaL->Fill( Edires_costhetaL[i] );
		Edi_costhetaK->Fill( Edires_costhetaK[i] );
		Edi_phi->Fill( Edires_phi[i] );
		Edi_time->Fill( Edires_time[i] );
		Edi_mass->Fill( Edires_mass[i] );

		HD_costhetaL->Fill( HDres_costhetaL[i] );
		HD_costhetaK->Fill( HDres_costhetaK[i] );
		HD_phi->Fill( HDres_phi[i] );
		HD_time->Fill( HDres_time[i] );
		HD_mass->Fill( HDres_mass[i] );

		mass_diff->Fill( (*MastermassData)[i] , sWeightDiff[i] );
	}

	Edi_costhetaL->Write("",TObject::kOverwrite);
	Edi_costhetaK->Write("",TObject::kOverwrite);
	Edi_phi->Write("",TObject::kOverwrite);
	Edi_mass->Write("",TObject::kOverwrite);
	Edi_time->Write("",TObject::kOverwrite);

	HD_costhetaL->Write("",TObject::kOverwrite);
	HD_costhetaK->Write("",TObject::kOverwrite);
	HD_phi->Write("",TObject::kOverwrite); 
	HD_mass->Write("",TObject::kOverwrite);
	HD_time->Write("",TObject::kOverwrite);

	mass_diff->Write("",TObject::kOverwrite);

	residuals->Close();
	return 0;
}


