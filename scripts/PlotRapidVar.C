#include <string>
#include <vector>
#include <iostream>

#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"

/* 
* Plot the rapid fit variables. 
* .L PlotRapidVar.C+ 
*  PlotRapidVar("theFile.root")
*   
*/


void plotHisto(TCanvas* tCan, int pad, TH1F* aHisto, const std::string& ax, bool logy = false, bool setMin = true){

  tCan->cd(pad);
  TVirtualPad* vpad = tCan->GetPad(pad);  
  if (logy) vpad->SetLogy();

  aHisto->SetTitle("");
  aHisto->Draw("HISTO");
  TAxis* xachse = aHisto->GetXaxis();

  xachse->SetTitleFont (132);
  xachse->SetLabelFont (132);
  xachse->SetTitle(ax.c_str());
  aHisto->GetYaxis()->SetLabelFont(132); 
  aHisto->GetYaxis()->SetTitleFont(132); 
  aHisto->GetYaxis()->SetTitle("# events"); 
  aHisto->GetYaxis()->SetTitleOffset(1.2);
  if (setMin) aHisto->SetMinimum(0);
  aHisto->Draw("HISTO");
  
}

void PlotRapidVar(std::string fileName = "theName.root", std::string treeName = "DecayTree"){

  // book the histograms
  TH1F* thisto = new TH1F("time", "time/ps", 100, -10., 10);
  TH1F* mhisto = new TH1F("mass", "mass/MeV", 100, 5100, 5600);
  TH1F* psihisto = new TH1F("cospsi", "cos(#psi)", 100, -1, 1);
  TH1F* thetahisto = new TH1F("costheta", "cos(#theta)", 100, -1, 1);
  TH1F* phihisto = new TH1F("phi", "#phi", 100, 0., 3.14);
  TH1F* taghisto = new TH1F("tag", "tag", 11, -5.5, 5.5);
  TH1F* mistaghisto = new TH1F("mistag", "mistag", 100, 0, 1);


  // get the input
  TChain* decaytree = new TChain(treeName.c_str());
  decaytree->Add(fileName.c_str());

  // address of the variables
  float time;decaytree->SetBranchAddress("time",&time); 
  float mass;decaytree->SetBranchAddress("mass",&mass);  
  float costheta;decaytree->SetBranchAddress("cosTheta",&costheta);   
  float cospsi;decaytree->SetBranchAddress("cosPsi",&cospsi);   
  float phi;decaytree->SetBranchAddress("phi",&phi);   
  float tag;decaytree->SetBranchAddress("tag",&tag);   
  float mistag;decaytree->SetBranchAddress("mistag",&mistag);   
 

  // loop and fill histograms
  int num_entries  = decaytree->GetEntries();
  for (int i = 0; i <num_entries; ++i) {
    decaytree->GetEntry(i);
    thisto->Fill(time);
    mhisto->Fill(mass);
    psihisto->Fill(cospsi);
    thetahisto->Fill(costheta);
    phihisto->Fill(phi);
    taghisto->Fill(tag);
    mistaghisto->Fill(mistag);
  } // loop i
 
  TCanvas *can = new TCanvas("TCan", "Tcan",10,44,800,800);
  can->Divide(2,4);

  // the t should be on log scale
  plotHisto(can, 1 , thisto, "t/ps", true, false);
 
  can->cd(2);
  plotHisto(can, 2, mhisto, "M/MeV");

  can->cd(3);
  plotHisto(can, 3, thetahisto, "cos(#theta)");
  
  can->cd(4);
  plotHisto(can, 4, psihisto, "cos(#psi)");

  can->cd(5);
  plotHisto(can, 5, phihisto, "#psi");

  can->cd(6);
  plotHisto(can, 6, taghisto, "tag");

  can->cd(7);
  plotHisto(can, 7, mistaghisto, "mistag");


}
