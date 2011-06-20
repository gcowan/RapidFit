#include "TFile.h"
#include "TPaveText.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TString.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TGraph2D.h"
#include "TMarker.h"
#include <cstdlib>
#include <stdlib.h>
#include <iostream>

using namespace::std;



TPaveText* addLHCbLabel(TString footer){
	TPaveText * label = new TPaveText(0.70,0.12,0.9,0.32,"BRNDC");
	//TPaveText * label = new TPaveText(0.12, 0.58, 0.12, 0.43,"BRNDC");
	//label->SetFillColor(0);
	label->SetFillStyle(0);
	label->SetBorderSize(0);     
	label->SetTextAlign(11);          
	label->SetTextSize(Float_t(0.04));
	TText * labeltext = 0;
	labeltext = label->AddText("LHC#font[12]{b} 2011 Data");
	labeltext = label->AddText("#sqrt{s} = 7TeV");
	labeltext = label->AddText(footer);
	return label;
}

int main (int argc, char* argv[] )
{
	cout << endl << "Really, REALLY stupid plotter for producing comparative plots between Edinburgh fit outputs." << endl;
	cout << endl << "USAGE: merge_plot Output_1.root Output_2.root \"Title 1\" \"Title 2\" " << endl;
	gStyle->SetCanvasColor(0);
	gStyle->SetFillColor(0);
	gROOT->SetStyle("Plain");
	if( argc != 5) exit(-1);
	TCanvas* c3 = new TCanvas("throw2","throw");
	TFile* input_1 = new TFile( argv[1], "READ" );
	gDirectory->ls();
	TH2D* hist_1 = (TH2D*)gDirectory->Get("pllhist");//"Graph2D");//pllhist");//Graph2D_from_nllhist");
	TGraph2D* graph_1 = (TGraph2D*)gDirectory->Get("pllhist");//"Graph2D");//pllhist");//Graph2D_from_nllhist");
	graph_1->Draw();
	c3->Update();
	TFile* input_2 = new TFile( argv[2], "READ" );
	gDirectory->ls();
	TH2D* hist_2 = (TH2D*)gDirectory->Get("pllhist");//"fcnew");//pllhist");//lr_data");
	TGraph2D* graph_2 = (TGraph2D*)gDirectory->Get("pllhist");//"fcnew");//pllhist");//lr_data");
	graph_2->Draw();
	hist_2->Draw();
	c3->Update();

	TString Plot_Title_1( argv[3] );
	TString Plot_Title_2( argv[4] );

	TCanvas* c1 = new TCanvas("Output_Plot","Output_Plot",1680,1050);
	double pllconts[3] = {1.15,2.305,3.0};//,4.61};
	//double pll2[3] = {2.3,4.61,6};//,9.62};
	//double pll2[3] = {0.68,0.9,0.95};
	double confs[3] = {68.0,90.0,95.0};//,99.0};
	TList* contLevel = NULL;
	TGraph* curv     = NULL;
	TGraph* gc    = NULL;
	//gStyle->SetCanvasColor(0);
	//gStyle->SetPalette(1);
	//gROOT->SetStyle("Plain");
	//gROOT->ForceStyle();
//	gStyle->SetFrameBorderMode(0);
//	gStyle->SetCanvasBorderMode(0);
//	gStyle->SetPadBorderMode(0);
//	gStyle->SetPadColor(0);
//	gStyle->SetCanvasColor(0);
//	gStyle->SetStatColor(0);
//	gStyle->SetTitleFillColor(0);
//	gStyle->SetFillColor(0);
//	gStyle->SetFrameFillColor(0);
	//gStyle->SetFillStyle(0);
	//gROOT->ForceStyle();

	hist_1->SetContour(3,pllconts);
	hist_1->SetLineWidth(1);
	c1->cd();
	hist_1->Draw("cont2");
//	hist_1->GetXaxis()->SetRangeUser(-3.1,2);
	hist_1->SetContour(3,pllconts);
//	c1->cd();
	hist_1->Draw("cont2");
	//	StyleTH2D(hist_2);
//	c3->cd();
	hist_2->Draw();
//	hist_2->GetXaxis()->SetRangeUser(-3.1,2);
//	hist_2->Draw();
	hist_1->Draw("cont2");
//	hist_2->Draw("SAME");//cont2");
	hist_2->SetContour(3,pllconts);//pll2);//pllconts);//pll2);
	hist_2->SetLineColor(2);
	hist_2->SetLineWidth(1);
//	c3->cd();
	c1->Update();
//	hist_2->Draw("SAMEcont2");
//	hist_2->GetXaxis()->SetRangeUser(-3.1,1);
//	hist_2->SetContour(3,pll2);//pllconts);//pll2);
//	c1->cd();
//	hist_2->Draw("SAMEcont2");
//	c1->Update();

	TCanvas* c2 = new TCanvas( "Throw","Throw" );
	hist_1->Draw("CONT LIST");
	c2->Update();

	//addLHCbLabel("NLL Scan")->Draw();
	TObjArray *contObjArr = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	int TotalConts = contObjArr->GetSize();
	c1->cd();
	TLegend *leg = new TLegend(0.65,0.7,1.1,0.9);
	leg->SetHeader( Plot_Title_1 );
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	for(int i = 0; i < TotalConts; i++){
		TString confname = "";
		double cl = confs[i];
		confname +=cl;
		confname += "% C.L.";
		contLevel = (TList*)contObjArr->At(i);
		for(int j =0; j<contLevel->GetSize(); j++){
			curv = (TGraph*)contLevel->At(j);
			gc = (TGraph*)curv->Clone();
			if( i!=3 ) gc->SetLineColor(Color_t(i+2));
			else gc->SetLineColor(Color_t(i+3));
			gc->Draw("L");
		}
		leg->AddEntry(gc,confname, "L");
	}
	c1->cd();
	TLegend *leg2 = new TLegend(0.13,0.7,0.5,0.9);
	leg2->SetHeader( Plot_Title_2 );
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);

	c2->cd();
	hist_2->Draw("CONT LIST");
	c2->Update();
	contObjArr = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TotalConts = contObjArr->GetSize();
	for(int i = 0; i < TotalConts; i++){
		TString confname = "";
		double cl = confs[i];
		confname +=cl;
		confname += "% C.L.";
		contLevel = (TList*)contObjArr->At(i);
		for(int j =0; j<contLevel->GetSize(); j++){
			curv = (TGraph*)contLevel->At(j);
			gc = (TGraph*)curv->Clone();
			if( i!=3 ) gc->SetLineColor(Color_t(i+2));
			else gc->SetLineColor(Color_t(i+3));
			gc->SetLineStyle(Style_t(i+2));
			c1->cd();
			gc->Draw("L");
			c2->cd();
		}
		leg2->AddEntry(gc,confname, "L");

	}
//	c1->cd();

//	Double_t X_min_1 = strtod( argv[3], NULL );
//	Double_t Y_min_1 = strtod( argv[4], NULL );
//	Double_t X_min_2 = strtod( argv[5], NULL );
//	Double_t Y_min_2 = strtod( argv[6], NULL );
	
//	c1->cd();
//	TMarker* new_mark = new TMarker( X_min_1, Y_min_1, 20 );
//	new_mark->SetMarkerSize(1);
//	new_mark->SetMarkerColor(8);
	//new_mark->Draw("SAME");
//	TMarker* new_mark2 = new TMarker( X_min_2, Y_min_2, 20 );
//	new_mark2->SetMarkerSize(1);
//	new_mark2->SetMarkerColor(9);
	//new_mark2->Draw("SAME");
//	c1->Update();
	
	c1->cd();
	leg->Draw();
	leg2->Draw();
	addLHCbLabel("NLL Contours")->Draw();
	c1->Update();
	c1->Print("Output.png");
	c1->Print("Output.pdf");

	input_1->Close();
	input_2->Close();
	return 0;
}


