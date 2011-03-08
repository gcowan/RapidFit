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

typedef vector<Float_t> array1f;
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
	label->SetTextSize(0.04);
	TText * labeltext = 0;
	labeltext = label->AddText("LHC#font[12]{b} 2010 Data");
	labeltext = label->AddText("#sqrt{s} = 7TeV");
	labeltext = label->AddText(footer);
	return label;
}

TCanvas *makeTempCanvas(TH2* hist, TString labelname){
	TCanvas *tempCanvas = new TCanvas(labelname+" temperature",labelname+" temperature",1024,768);
	hist->Draw("colz");
	addLHCbLabel(labelname)->Draw();
	return tempCanvas;
}

TCanvas *makeContCanvas(TH2* hist, TString labelname){
	TCanvas *contourCanvas = new TCanvas(labelname+" contour",labelname+" contour",1024,768);
	hist->Draw("cont1z");
	addLHCbLabel(labelname)->Draw();
	return contourCanvas;
}
TCanvas *makeConfCanvas(TH2* hist, TString labelname,UInt_t nconts, double* conts, double* confs, bool pub){
	//TH2* hist = (TH2*)_hist->Clone("confhist");
	TString pname = "color";
	if(pub){pname = "pub";}
	TCanvas *confCanvas = new TCanvas(pname + " " + labelname+" CL contours",pname + " " + labelname+" CL contours",1024,768);
	hist->SetContour(nconts,conts);
	hist->Draw("cont LIST");
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
	hist->Draw("AXIS");
	for(int i = 0; i < TotalConts; i++){
		TString confname = "";
		double cl = confs[i];
		confname +=cl;
		confname += "\% C.L.";
		contLevel = (TList*)contObjArr->At(i);
		for(int j =0; j<contLevel->GetSize(); j++){
			curv = (TGraph*)contLevel->At(j);
			gc = (TGraph*)curv->Clone();
			if(pub){	
				gc->SetLineStyle(i+1);
			}else{
				gc->SetLineColor(i+2);
			}
			gc->Draw("L");
		}
		leg->AddEntry(gc,confname, "L");
	}
	addLHCbLabel(labelname)->Draw();
	leg->Draw();
	confCanvas->Update();
	//hist->Delete();
	hist->SetContour(20);
	return confCanvas;


}


int main(int argc, char *argv[]){

	if(argc !=5 ){
		cout << "Plots 2D FCScans from rapidfit flat ntuples. Usage:" << endl;
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
	cout << "INPUT NTUPLE CONTAINS: " << entries << "entries" << endl;


	TString param1valstr = param1string + "_value";
	TString param1genstr = param1string + "_gen";
	TString param1errstr = param1string + "_error";

	TString param2valstr = param2string + "_value";
	TString param2genstr = param2string + "_gen";
	TString param2errstr = param2string + "_error";

	TString FRstr = "Fit_Status";
	TString NLLstr = "NLL";


	//The first entry in every file will be the data result floated globally
	//Gen values will be -9999, 
	//Par values will be floated
	//Par errors !=0
	//The second entry in every file is the data result fixed at the current gridpoint. 
	//Gen values will be -9999, 
	//Par values will be set to the current gridpoint
	//Par errors will be 0
	//Entries 3-> Ntoyspersj+3 will be the fixed toy values
	//Gen values will be set to the current gridpoint
	//Par values will be set to the current gridpoint
	//Par errors will be 0
	//Entries Ntoyspersj+3->2*Ntoyspersj+3 will be the floated toy values
	//Gen values will be set to the current gridpoint
	//Par values will be floated
	//Par errors !=0

	TString notgen = "-9999.";
	Float_t param1val, param2val;
	Float_t nll, nlldatabest;
	array1f param1gridpoints, param2gridpoints, dataRatiogridpoints, clgridpoints;
	allresults->SetBranchAddress(NLLstr,&nll);
	// Find the global minimum from the first entry:

	allresults->GetEntry(0);
	nlldatabest = nll;
	cout << "GLOBAL DATA MINIMUM NLL: " << nlldatabest << endl;

	// Find the positions of the gridpoints we scanned, and find the profileLL at each point for data:
	TString datafixedstr = param1genstr + "==" + notgen + "&&" + param2genstr + "==" + notgen + "&&"+ param1errstr + "==0.0&&" +param2errstr +"==0.0&&" + FRstr + "==3.0";
	TTree* datafixed = allresults->CopyTree(datafixedstr);
	datafixed->SetBranchAddress(param1valstr,&param1val);
	datafixed->SetBranchAddress(param2valstr,&param2val);
	datafixed->SetBranchAddress(NLLstr,&nll);
	for(Long64_t i = 1; i < datafixed->GetEntries(); i++){
		datafixed->GetEntry(i);
		// We've found a new gridpoint
		// But is it unique?
		bool unique = true;
		for(UInt_t j = 0; j< param1gridpoints.size(); j++){
			if(param1gridpoints[j] == param1val && param2gridpoints[j] == param2val){unique = false;}
		}
		if(unique){
			param1gridpoints.push_back(param1val);
			param2gridpoints.push_back(param2val);
			dataRatiogridpoints.push_back(nll-nlldatabest);
			//cout << "GRIDPOINT FOUND: " << param1val << " " << param2val << " " << nll << endl;
		}
	}
	Float_t NLLtoyfixed, NLLtoyfloat;
	delete datafixed;
	output->cd();

	bool hastoys = false;
	TString isatoy = param1genstr +"!="+ notgen + "&&" + param2genstr +"!=" +notgen;
	TTree* toys = allresults->CopyTree(isatoy);
	cout << "FOUND " << toys->GetEntries() << " TOYS" << endl;
	//Is this just an LLscan, or were toys for FC generated? 
	if(toys->GetEntries() != 0){hastoys = true;}
	if(hastoys){
		//Find the toys generated at the previously discovered gridpoints
		for(UInt_t i = 0; i<param1gridpoints.size(); i++){
			//cout << "step " << i << " of " <<  param1gridpoints.size() << endl;
			// Build up cutstrings for the toys
			TString toyscutstr = FRstr + "==3.0&&" +param1genstr + "==";
			toyscutstr += param1gridpoints[i];
			toyscutstr += "&&" + param2genstr + "==";
			toyscutstr += param2gridpoints[i];
			toyscutstr += "&&";
			//Reduce the ntuple to only toys generated at this gridpoint, and only toys that had floated params
			TString floatedtoyscutstr = toyscutstr;
			floatedtoyscutstr += param1errstr + "!=0.0&&" + param2errstr + "!=0.0";
			TTree *floatedtoys = toys->CopyTree(floatedtoyscutstr);
			UInt_t floatedtoystot = floatedtoys->GetEntries();
			//Do the same for toys that had fixed params
			TString fixedtoyscutstr = toyscutstr;
			fixedtoyscutstr += param1errstr + "==0.0&&" + param2errstr + "==0.0";
			TTree* fixedtoys = toys->CopyTree(fixedtoyscutstr);
			UInt_t fixedtoystot = fixedtoys->GetEntries();
			if(fixedtoystot != floatedtoystot){
				//We've got serious problems! 
				cout << "Different number of toy fits found for gridpoint " << param1gridpoints[i] << "," << param2gridpoints[i] << " (" << fixedtoystot << "," << floatedtoystot << ") Can't continue!" << endl;

				exit(1);
			}
			if(fixedtoystot != 0){
			//Loop over the toys, pulling out the NLL ratio
			floatedtoys->SetBranchAddress(NLLstr,&NLLtoyfloat);	
			fixedtoys->SetBranchAddress(NLLstr,&NLLtoyfixed);
			UInt_t toyNLLsmaller = 0;
			for(UInt_t j = 0; j<fixedtoystot; j++){
				floatedtoys->GetEntry(j);	
				fixedtoys->GetEntry(j);	
				//FIXME THE LINE BELOW IS THE FELDMAN-COUSINS ORDERING METHOD USED BY CDF/HEIDELBERG: FIXME 
				//if the toyratio is smaller than the data ratio at this point, increment:
				if((NLLtoyfixed-NLLtoyfloat)<dataRatiogridpoints[i]){toyNLLsmaller++;}
			}
			//The C.L. is the percentage of toys that were smaller
			double cl = ((Double_t)toyNLLsmaller/(Double_t)fixedtoystot);
			clgridpoints.push_back(cl);
			}else{ 
				cout << param1gridpoints[i] << " " << param2gridpoints[i] << " WARNING: Found no toys here. " << endl;
				clgridpoints.push_back(-9999.0);
			}

			delete floatedtoys;
			delete fixedtoys;
		}
	}
	delete toys;
	//We now have 4 vectors: The gridpoints in x,y and the conf. limits in z or the profile LL in z. We need to make TGraphs from these bad boys. 

	//Copy vectors to arrays: 
	int npoints = param1gridpoints.size();
	Double_t* p2points = new Double_t [npoints];
	copy( param2gridpoints.begin(), param2gridpoints.end(),p2points);
	Double_t* p1points = new Double_t [npoints];
	copy( param1gridpoints.begin(), param1gridpoints.end(),p1points);
	Double_t* clpoints = new Double_t [npoints];
	copy( clgridpoints.begin(), clgridpoints.end(),clpoints);
	Double_t* pllpoints = new Double_t [npoints];
	copy( dataRatiogridpoints.begin(), dataRatiogridpoints.end(),pllpoints);

	output->cd();

	int np = (int)sqrt((double)npoints);

	//We get the data profile likelihood free, so let's plot it: 
	TGraph2D *pllgraph = new TGraph2D(npoints, p2points, p1points, pllpoints);
	pllgraph->SetName("pllgraph");
	pllgraph->SetNpx(np);
	pllgraph->SetNpy(np);

	TH2D* pllhist = pllgraph->GetHistogram();
	pllhist->SetName("pllhist");
	pllhist->SetTitle("");
	pllhist->GetXaxis()->SetTitle(param2string);
	pllhist->GetYaxis()->SetTitle(param1string);
	pllhist->Write();


	double pllconts[4] = {1.15,2.36,3.0,4.61};
	double fcconts[4] = {0.68,0.9,0.95,0.99};
	double confs[4] = {68.0,90.0,95.0,99.0};

	TCanvas *plltemp = makeTempCanvas(pllhist,"Profile Likelihood");
	plltemp->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_temp.pdf");
	plltemp->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_temp.png");
	TCanvas *pllcont = makeContCanvas(pllhist,"Profile Likelihood");
	pllcont->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_cont.pdf");
	pllcont->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_cont.png");
	TCanvas *pllconf = makeConfCanvas(pllhist,"Profile Likelihood",4,pllconts,confs,false);
	pllconf->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_conf.pdf");
	pllconf->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_conf.png");
	TCanvas *pllpub = makeConfCanvas(pllhist,"Profile Likelihood",4,pllconts,confs,true);
	pllpub->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_pub.pdf");
	pllpub->Print(outputdir+"/"+param1string+"_"+param2string+"_pll_pub.png");

	// Now let's plot the feldman cousins in the same way: 
	if(hastoys){
		TGraph2D *fcgraph = new TGraph2D(npoints, p2points, p1points, clpoints);
		fcgraph->SetName("fcgraph");
		fcgraph->SetNpx(np);
		fcgraph->SetNpy(np);
		TH2D* fchist = fcgraph->GetHistogram();
		fchist->SetName("fchist");
		fchist->SetTitle("");
		fchist->GetXaxis()->SetTitle(param2string);
		fchist->GetYaxis()->SetTitle(param1string);
		fchist->Write();

		TCanvas *fctemp = makeTempCanvas(fchist,"Feldman Cousins");
		fctemp->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_temp.pdf");
		fctemp->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_temp.png");
		TCanvas *fccont = makeContCanvas(fchist,"Feldman Cousins");
		fccont->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_cont.pdf");
		fccont->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_cont.png");
		TCanvas *fcconf = makeConfCanvas(fchist,"Feldman Cousins",4,fcconts,confs,false);
		fcconf->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_conf.pdf");
		fcconf->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_conf.png");
		TCanvas *fcpub = makeConfCanvas(fchist,"Feldman Cousins",4,fcconts,confs,true);
		fcpub->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_pub.pdf");
		fcpub->Print(outputdir+"/"+param1string+"_"+param2string+"_fc_pub.png");
		output->Write();
		output->Close();
	}

}

