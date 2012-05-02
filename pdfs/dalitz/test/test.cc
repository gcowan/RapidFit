/*
 * This is small program to run some tests of amplitude calculator
 */

#include "DPTotalAmplitude.hh"

#include "TTree.h"
#include "TFile.h"

struct dataStruct
{
// Angles
  float theta1;
  float theta2;
  float phi;
// Invariand masses
  float mKpi;
  float mJpsipi;
  float mJpsiK;
// Store also invariant masses of final state particles (debugging reasons)
  float mK;
  float mpi;
  float mmp;
  float mmm;
  float mJpsi1;
  float mJpsi2;
// Weight
  float wt;
// Do we have K* or K*bar
  int pid;
} ;

char dataStructString[]="cosTheta1/F:cosTheta2/F:phi/F:mKpi/F:mJpsipi/F:mJpsiK/F:mK/F:mpi/F:mmp/F:mmm/F:mJpsi1/F:mJpsi2/F:w/F:id/I";

int main(int argc, char** argv)
{

  DPTotalAmplitude calculator;

  // Open file with phase space decays and get tree from it
  TFile* input=TFile::Open("../PhaseSpacePreparation/phsp_B0JpsiKspin1.root"); 
  TTree* inputTree=(TTree*)input->Get("tree"); 
  dataStruct infoToStore; 
  inputTree->SetBranchAddress("angles", (void*)&infoToStore);

  // Output tree
  TFile* output=TFile::Open("test/testJpsiK892.root","recreate");
  TTree *tree=new TTree("treeOut","");
  tree->Branch("angles",&infoToStore,dataStructString);

  // Event loop
//  for (int i=0; i<inputTree->GetEntries();++i)
  for (int i=0; i<100000;++i)
  {
    inputTree->GetEntry(i);

    infoToStore.wt=calculator.matrixElement(infoToStore.mKpi,
                   infoToStore.theta1, infoToStore.theta2, infoToStore.phi);

    tree->Fill();
  }

  tree->Write("tree");
  output->Close();
  
}
