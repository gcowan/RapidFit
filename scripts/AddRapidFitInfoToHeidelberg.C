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


std::string AddRapidFitInfoToHeidelberg(std::string fileName = "heidelberg-data.root", std::string treeName = "Bs2JpsiPhi", bool setTagToZero = false){

  // get the input
  TChain* decaytree = new TChain(treeName.c_str());
  decaytree->Add(fileName.c_str());

	
  // make the output
  std::string outputName = fileName.substr(0,fileName.size() - 5);
 
  if (setTagToZero == false){
     outputName += "_RapidFit.root";
  }
  else {
     outputName += "_RapidFit_Untagged.root";
  }
 
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
  double dtime; newtree->SetBranchAddress("t",&dtime); 
  double dcosTheta; newtree->SetBranchAddress("trcostheta",&dcosTheta); 
  double dcosPsi; newtree->SetBranchAddress("trcospsi",&dcosPsi);
  double dphi; newtree->SetBranchAddress("trphi",&dphi);
  int itag; newtree->SetBranchAddress("tagdecision",&itag);
  double dmistag; newtree->SetBranchAddress("tagomega",&dmistag);
  double dmass; newtree->SetBranchAddress("m",&dmass);

  // loop and make the conversion
  int num_entries  = newtree->GetEntries();
  for (int i = 0; i <num_entries; ++i) {
    newtree->GetEntry(i);
    time = dtime;
    cosTheta = dcosTheta;
    cosPsi = dcosPsi;
    phi = dphi;
    if (setTagToZero == false){ 
      tag = (float) itag;
    }
    else {
      tag = 0.0;
    }    
    mass = dmass;
    mistag = dmistag;
  
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
