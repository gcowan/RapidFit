#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TMath.h"

/*
*  Add the rapid fit variables to the standard ntuple. To use, type then root and at the prompt
*   .L AddRapidFitInfo.C++
*   AddRapidFitInfo("TheFile.root", "DecayTree")
*   It will write out TheFile_RapidFit.root
*/


std::string AddRapidFitInfo(std::string fileName = "theFile.root", std::string treeName = "DecayTree"){

  // get the input
  TChain* decaytree = new TChain(treeName.c_str());
  decaytree->Add(fileName.c_str());


  // make the output
  std::string outputName = fileName.substr(0,fileName.size() - 5);
  outputName += "_RapidFit.root";
  TFile* outFile  =new TFile(outputName.c_str(),"RECREATE");

  std::cout << "Reading: " << treeName << " from " << fileName  << " to " << outputName << std::endl;

  // clone the tree..
  TTree*  newtree = decaytree->CloneTree(-1);

  //add the new branches
  float time; TBranch* branch_time =  newtree->Branch("time",&time, "time/F");
  float cosTheta; TBranch* branch_cosTheta =  newtree->Branch("cosTheta",&cosTheta, "cosTheta/F");
  float cosPsi; TBranch* branch_cosPsi =  newtree->Branch("cosPsi",&cosPsi, "cosPsi/F");
  float phi; TBranch* branch_phi =  newtree->Branch("phi",&phi, "phi/F");
  float tag; TBranch* branch_tag =  newtree->Branch("tag",&tag, "tag/F");
  float mistag; TBranch* branch_mistag =  newtree->Branch("mistag",&mistag, "mistag/F");
  float mass; TBranch* branch_mass =  newtree->Branch("mass",&mass, "mass/F");

  // address of the old variables
  newtree->SetBranchAddress("B_s0_LOKI_DTF_CTAU",&time);
  newtree->SetBranchAddress("B_s0_ThetaTr",&cosTheta);
  newtree->SetBranchAddress("B_s0_ThetaVtr",&cosPsi);
  newtree->SetBranchAddress("B_s0_PhiTr",&phi);
  int itag; newtree->SetBranchAddress("B_s0_TAGDECISION",&itag);
  newtree->SetBranchAddress("B_s0_TAGOMEGA",&mistag);
  newtree->SetBranchAddress("B_s0_LOKI_MASS_JpsiConstr",&mass);

  // loop and make the conversion
  int num_entries  = newtree->GetEntries();
  for (int i = 0; i <num_entries; ++i) {
    newtree->GetEntry(i);
    time /= 0.299792458;
	tag = (float) itag;
    cosTheta = TMath::Cos(cosTheta);
    cosPsi = TMath::Cos(cosPsi);
    branch_time->Fill();
    branch_cosTheta->Fill();
    branch_cosPsi->Fill();
    branch_phi->Fill();
    branch_tag->Fill();
    branch_mistag->Fill();
    branch_mass->Fill();
  } // loop i

  newtree->Write();
  outFile->Close();

  return outputName;
}
