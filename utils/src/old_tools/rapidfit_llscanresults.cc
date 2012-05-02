/* stackergen: Part of the simpletools package
 * (c) Conor Fitzpatrick, 2008
 *
 * If you find this program useful in whole or in part 
 * please cite this paper: 
 *
 * Feel free to send bugreports, feature requests, patches etc to:
 * conor.fitzpatrick@cern.ch
 *
 */

#include <stdlib.h>
#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCint.h>
#include <TMath.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <fstream>
#include "stdio.h"
#include "string"
#include "TStyle.h"
#include "Riostream.h"
#include "TAxis.h"
#include <TRandom3.h>
#include <cctype>
#include <cmath>
#include <vector>
#include <TArrow.h>
#include <TEntryList.h>
#include <TLegend.h>
#include<THStack.h>
#include<TColor.h>
#include<Rtypes.h>
#include<TObjArray.h>
#include<TKey.h>
#include<TVirtualHistPainter.h>
#include<TH2F.h>
#include<TGraph2D.h>
#include<TPaveText.h>

using std::cout;
using std::endl;

typedef vector<Double_t> array1d;
typedef vector<TString> array1s;

inline TString prettyPrint(Double_t value){
	char pretty[20];
	TString prettyString;
	sprintf (pretty, "%1.3g",value);
	prettyString = pretty;
	return prettyString;
}


inline TString prettyPrint(Double_t val, Double_t err){
	TString outstr = "$";
	char errstr[20];
	char valstr[20];
	Int_t n = sprintf (errstr, "%1.1g",err);
	sprintf(valstr,"%1.*f",n-2,val);
	outstr += valstr;
	outstr += "\\pm";
	outstr += errstr;
	outstr+= "$";
	return outstr;


}

TPaveText* addLHCbLabel(TString footer){
	TPaveText * label = new TPaveText(0.14, 0.88, 0.14, 0.73,"BRNDC");
	//TPaveText * label = new TPaveText(0.12, 0.58, 0.12, 0.43,"BRNDC");
	label->SetFillColor(0);
	label->SetBorderSize(0);     
	label->SetTextAlign(11);          
	label->SetTextSize(Float_t(0.04));
	TText * labeltext = 0;
	labeltext = label->AddText("LHC#font[12]{b} 2011 Data");
	labeltext = label->AddText("#sqrt{s} = 7TeV");
	labeltext = label->AddText(footer);
	return label;
}


int main(int argc, char *argv[]){

	if(argc !=5 ){
		cout << "Plots 2D LLScans from rapidfit flat ntuples. Usage:" << endl;
		cout <<  argv[0] << " <flatntuple.root> <param1> <param2> <outputdir>"	<< endl;
		exit(1);
	}
	gStyle->SetPalette(1);
	gROOT->SetStyle("Plain");
	gStyle->SetFrameBorderMode(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetStatColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetFillColor(0);
	gStyle->SetFrameFillColor(0);


	TFile * input = TFile::Open( argv[1] );
	TString outputdir = argv[4];

	TString param1string = argv[2];
	TString param2string = argv[3];
	gSystem->mkdir( outputdir );
	TFile * output = new TFile(outputdir+"/llscanresults.root","RECREATE");
	input->cd();
	TNtuple* allresults;
	input->GetObject("RapidFitResult", allresults);
	if(!allresults){
		input->GetObject("RapidFitResult/RapidFitResult",allresults);
		if(!allresults){
			cout << "Couldn't find ntuple RapidFitResult in TFile" << endl;
			exit(1);
		}
	}


	Long64_t entries = allresults->GetEntries();
	TString Cut_String(param1string);
	Cut_String.Append("_gen>-9999");
	TNtuple *intermediate_results = (TNtuple*)allresults->CopyTree("");//"Fit_Status==3");
	TNtuple *results = (TNtuple*)intermediate_results->CopyTree(Cut_String);
	Long64_t goodentries = results->GetEntries();
	cout << "Warning: Did not use " << entries - goodentries << " results as these had a bad fit status" << endl;

	Float_t p1maximum = Float_t(results->GetMaximum(param1string+"_value"));
	Float_t p1minimum = Float_t(results->GetMinimum(param1string+"_value"));
	Float_t p2maximum = Float_t(results->GetMaximum(param2string+"_value"));
	Float_t p2minimum = Float_t(results->GetMinimum(param2string+"_value"));
	Float_t minNLL = Float_t(intermediate_results->GetMinimum("NLL"));
	TString minNLLString = "";
	minNLLString += minNLL;
	cout << "plot ranges are: " << param1string << ": " << p1minimum << " " << p1maximum << " " << param2string << ": " << p2minimum << " " << p2maximum << endl;
	cout << "Minimum NLL is: " << minNLL << endl;	

	TString bins = "(200,0,0,200,0,0)";

	results->Draw(param1string+"_value:"+param2string+"_value>>nllhist"+bins,"(NLL-"+minNLLString+")");
	//TGraph2D *nllgraph = (TGraph2D*)gDirectory->Get("nllhist");
	TGraph2D *nllgraph = new TGraph2D((TH2*)gDirectory->Get("nllhist"));
	nllgraph->SetNpx(100);
	nllgraph->SetNpy(100);
	//Double_t * x = nllgraph->GetX();
	//Double_t * y = nllgraph->GetY();
	//Double_t * z = nllgraph->GetZ();
	TH2D* nllhist = nllgraph->GetHistogram();


	output->cd();
	nllhist->SetTitle("");
	nllhist->GetXaxis()->SetTitle(param2string);
	nllhist->GetYaxis()->SetTitle(param1string);
	nllhist->Write();

	TCanvas *tempCanvas = new TCanvas("temperature plot", "temperature plot",2048,1536);
	nllhist->Draw("colz");

	addLHCbLabel("NLL Scan")->Draw();
	tempCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_temp.pdf");
	tempCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_temp.eps");
	tempCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_temp.png");

	tempCanvas->Write();

	TCanvas *contCanvas = new TCanvas("contour plot", "contour plot",2048,1536);
	nllhist->Draw("cont1z");
	addLHCbLabel("NLL Contours")->Draw();
	contCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_cont.pdf");
	contCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_cont.eps");
	contCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_cont.png");
	contCanvas->Write();

	double conts[4] = {1.15,2.31,3.0,4.61};
	double confs[4] = {68.0,90.0,95.0,99.0};
	nllhist->SetContour(4,conts);

	TCanvas *confCanvas = new TCanvas("onfidence interval plot","confidence interval plot",2048,1536);
	nllhist->Draw("cont LIST");
	confCanvas->Update();
	TObjArray *contObjArr = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TList* contLevel = NULL;
	TGraph* curv     = NULL;
	TGraph* gc    = NULL;
	int TotalConts = contObjArr->GetSize();
	TLegend *leg = new TLegend(0.80,0.89,0.95,0.7);
	leg->SetHeader("Conf. Levels");
	leg->SetBorderSize(0);
        leg->SetFillStyle(0);	
	for(int i = 0; i < TotalConts; i++){
		TString confname = "";
		double cl = confs[i];
		confname +=cl;
		confname += "\% C.L.";
		contLevel = (TList*)contObjArr->At(i);
		curv = (TGraph*)contLevel->First();
		gc = (TGraph*)curv->Clone();
		leg->AddEntry(gc,confname,"F");
		//	gc->Draw("L");
	}
	confCanvas->Clear();
	nllhist->Draw("cont1");
	addLHCbLabel("NLL Scan")->Draw();
	leg->Draw();
	confCanvas->Update();

	confCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_conf.pdf");
	confCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_conf.eps");
	confCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_conf.png");
	confCanvas->Write();

	TCanvas *pubCanvas = new TCanvas("onfidence interval plot (pub)","confidence interval plot (pub)",2048,1536);
	nllhist->Draw("cont2 LIST");
	pubCanvas->Update();
	contObjArr = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	pubCanvas->Clear();
	contLevel = NULL;
	curv     = NULL;
	gc = NULL;
	TotalConts = contObjArr->GetSize();
	leg = new TLegend(0.80,0.89,0.95,0.7);
	leg->SetHeader("Conf. Levels");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	nllhist->Draw("AXIS");
	for(int i = 0; i < TotalConts; i++){
		TString confname = ""; 
		double cl = confs[i];
		confname +=cl;
		confname += "\% C.L.";
		contLevel = (TList*)contObjArr->At(i);
		for(int j =0; j<contLevel->GetSize(); j++){
		curv = (TGraph*)contLevel->At(j);
		gc = (TGraph*)curv->Clone();
		gc->SetLineStyle(Style_t(i+1));
		gc->Draw("L");
		}
		leg->AddEntry(gc,confname, "L");
	}


	addLHCbLabel("NLL Scan")->Draw();
	leg->Draw();
	pubCanvas->Update();
	pubCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_pub.pdf");
	pubCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_pub.eps");
	pubCanvas->Print(outputdir+"/"+param1string+"_"+param2string+"_pub.png");
	pubCanvas->Write();


	output->Write();
	output->Close();


} 
