#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"

/*
*  Example of how to write a reduced ntuple
*  Two selections, Matts (add on to roadmap):
*  Use: .L ApplySelection.C+
*        ApplySelection("theFile.root")
*
* Yuehongs implementation of the roadmap
*  Use .L ApplySelection.C+
*      ApplyRoadMap("theFile.root")
*
*/

/** First some helpers to avoid duplicating code */

std::string createOutputName(std::string& name, std::string& trailer){
  // helper for making the output name
  std::string outputName = name.substr(0,name.size() - 5);
  outputName += trailer;
  return outputName;
}

TFile* openOutput(std::string& tree, std::string& input, std::string& output) {
  // helper for opening the file
  TFile* outFile  =new TFile(output.c_str(),"RECREATE");
  std::cout << "Reading: " << tree << " from " << input  << " to " << output << std::endl;
  return outFile;
}

void finalize(TTree* cutTree, TFile* output){
  // helper for finalizing
  TTree*  newtree = cutTree->CloneTree(-1);
  newtree->Write();
  output->Close();
}

/* Now the main business */

std::string ApplySelection(std::string fileName = "theFile.root", std::string treeName = "DecayTree", std::string trailer = "_Selected.root"){

  // Matts selection, generally applied on top of the road map

  // get the input
  TChain* decaytree = new TChain(treeName.c_str());
  decaytree->Add(fileName.c_str());

  // make the output file name
  std::string outputName = createOutputName(fileName, trailer);

  // make the output file
  TFile* outFile = openOutput(treeName,fileName,outputName);

  // the cut strings to apply
  TCut cutJpsi = "TMath::Abs(J_psi_1S_MM - 3096.7)/J_psi_1S_MMERR < 4.5";
  TCut cutK = "Kplus_TRACK_CHI2NDOF<4&&Kminus_TRACK_CHI2NDOF<4";
  TCut qualityCut = "B_s0_LOKI_DTF_CHI2NDOF > 0 &&B_s0_LOKI_DTF_VCHI2NDOF > 0";
  TCut oneCand = "hasBestVtxChi2 > 0";
  TCut cleanTail = "(B_s0_MINIPCHI2NEXTBEST > 50 || B_s0_MINIPCHI2NEXTBEST < 0 )";
  //  TCut cut0 = "B0_LOKI_MASS_JpsiConstr>5275&&B0_LOKI_MASS_JpsiConstr<5450";
  //  TCut cut0 = "B_s0_LOKI_MASS_JpsiConstr>5300&&B_s0_LOKI_MASS_JpsiConstr<5450";

  // make the reduced tree
  TTree* smalltree = decaytree->CopyTree(oneCand&&cutK&&qualityCut&&cutJpsi&&cleanTail);
  finalize(smalltree,outFile);

  return outputName;
}

std::string ApplyRoadMap(std::string fileName = "theFile.root", std::string treeName = "DecayTree", std::string trailer = "_Roadmap.root"){

  // Matts selection, generally applied on top of the road map

  // get the input
  TChain* decaytree = new TChain(treeName.c_str());
  decaytree->Add(fileName.c_str());

  // make the output file name
  std::string outputName = createOutputName(fileName, trailer);

  // make the output file
  TFile* outFile = openOutput(treeName,fileName,outputName);

  // the cut strings to apply
  TCut cutMu = "min(muplus_PIDmu,muminus_PIDmu)>0&&max(muplus_TRACK_CHI2NDOF,muminus_TRACK_CHI2NDOF)<4&&min(muplus_PT,muminus_PT)>500";
  TCut cutJ = "J_psi_1S_ENDVERTEX_CHI2/J_psi_1S_ENDVERTEX_NDOF<11";
  TCut cutK = "Kplus_PIDK>0&&Kplus_TRACK_CHI2NDOF<10&&Kminus_PIDK>0&&Kminus_TRACK_CHI2NDOF<10";
  TCut cutPhi = "phi_1020_PT>1000&&phi_1020_ENDVERTEX_CHI2/phi_1020_ENDVERTEX_NDOF<20&&abs(phi_1020_MM-1020)<12";
  TCut cutB = "B_s0_IPCHI2_OWNPV<25&&B_s0_LOKI_MASS_JpsiConstr>5200&&B_s0_LOKI_MASS_JpsiConstr<5550";

  // make the reduced tree
  TTree* smalltree = decaytree->CopyTree(cutMu&&cutJ&&cutK&&cutPhi&&cutB);
  finalize(smalltree,outFile);

  return outputName;
}
