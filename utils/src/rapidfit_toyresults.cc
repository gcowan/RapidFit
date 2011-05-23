/* stackergen: Part of the simpletools package
 *  * (c) Conor Fitzpatrick, 2008
 *   *
 *    * If you find this program useful in whole or in part 
 *     * please cite this paper: 
 *      *
 *       * Feel free to send bugreports, feature requests, patches etc to:
 *        * conor.fitzpatrick@cern.ch
 *         *
 *          */

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
#include<TF1.h>

using std::cout;
using std::endl;

typedef vector<Double_t> array1d;
typedef vector<TString> array1s;

inline TString prettyPrint(Double_t value){
	char pretty[20];
	TString prettyString;
	sprintf (pretty, "%1.1g",value);
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

int main(int argc, char *argv[]){

	if(argc !=3){
		cout << "Plots, fits and tabulates the mean, width, pull of RapidFit toyresults automagically. Usage:" << endl;
		cout <<  argv[0] << " <flatntuple.root> <outputdir>"    << endl;
		exit(1);
	}
	gROOT->SetStyle("Plain");
	gStyle->SetTitleStyle(1);
	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(1111);

	array1s param;
	array1d mean;
	array1d emean;
	array1d sigma;
	array1d esigma;
	array1d pull;
	array1d epull;
	array1d pullw;
	array1d epullw;

	TFile * input = TFile::Open( argv[1] );
	TString outputdir = argv[2];
	gSystem->mkdir( outputdir );
	TFile * output = new TFile(outputdir+"/results.root","RECREATE");
	input->cd();
	TNtuple* results;
	input->GetObject("RapidFitResult", results);
	if(!results){
		input->GetObject("RapidFitResult/RapidFitResult",results);
		if(!results){
			cout << "Couldn't find ntuple RapidFitResult in TFile" << endl;
			exit(1);
		}
	}
	TObjArray *vars = results->GetListOfLeaves();
	for(int i = 0; i < vars->GetEntries(); i++){
		TString varname = (vars->At(i)->GetName());
		if(varname.EndsWith("_value")){
			TString var_value = varname;
			varname.Resize(varname.Sizeof()-7);
			cout << varname << endl;
			TString var_pull = varname;
			var_pull += "_pull";
			TString var_error = varname;
			var_error += "_error";

			output->cd();
			results->Draw(var_pull+">>pull","Fit_Status==3");
			TH1F *pullhist = (TH1F*)gDirectory->Get("pull");
			pullhist->Sumw2();

			if(pullhist->GetRMS() > 0.0){
				param.push_back(varname);
				TCanvas *pullCanvas = new TCanvas(varname + " Pull",varname + " Pull",2048,1536);
				pullhist->SetTitle(varname + " Pull");
				pullhist->Draw();

				TF1 *pullgauss = new TF1(varname + "Pull Gaussian","gaus",pullhist->GetXaxis()->GetXmin(),pullhist->GetXaxis()->GetXmax());
				pullhist->Fit(pullgauss);
				pull.push_back(pullgauss->GetParameter(1));
				epull.push_back(pullgauss->GetParError(1));
				pullw.push_back(pullgauss->GetParameter(2));
				epullw.push_back(pullgauss->GetParError(2));
				pullgauss->Draw("SAME");
				pullCanvas->Draw();
				pullCanvas->Write();
				pullCanvas->Print(outputdir+"/"+varname+"_pull.pdf");
				pullCanvas->Print(outputdir+"/"+varname+"_pull.eps");
				pullgauss->Delete();	

				results->Draw(var_error+">>error","Fit_Status==3");
				TCanvas *errCanvas = new TCanvas(varname + " Error",varname + " Error",2048,1536);
				TH1F *errhist = (TH1F*)gDirectory->Get("error");
				errhist->Sumw2();
				errhist->SetTitle(varname + " Error");
				errhist->Draw();
				errCanvas->Draw();
				errCanvas->Write();
				errCanvas->Print(outputdir+"/"+varname+"_error.pdf");
				errCanvas->Print(outputdir+"/"+varname+"_error.eps");
				errhist->Delete();


				results->Draw(var_value+">>value","Fit_Status==3");
				TH1F *valhist = (TH1F*)gDirectory->Get("value");
				valhist->Sumw2();
				valhist->SetTitle(varname + " Value");
				TCanvas *valCanvas = new TCanvas(varname + " Value",varname + " Value",2048,1536);
				valhist->Draw();
				TF1 *valgauss = new TF1(varname + "Value Gaussian","gaus",valhist->GetXaxis()->GetXmin(),valhist->GetXaxis()->GetXmax());
				valhist->Fit(valgauss);
				mean.push_back(valgauss->GetParameter(1));
				emean.push_back(valgauss->GetParError(1));
				sigma.push_back(valgauss->GetParameter(2));
				esigma.push_back(valgauss->GetParError(2));
				valgauss->Draw("SAME");
				valCanvas->Draw();
				valCanvas->Write();
				valCanvas->Print(outputdir+"/"+varname+"_value.pdf");
				valCanvas->Print(outputdir+"/"+varname+"_value.eps");
				valhist->Delete();
				valgauss->Delete();


			}
			pullhist->Delete();
		}
	}


	output->Write();
	output->Close();

	cout << "\\begin{table}" << endl;
	cout << "\\begin{tabular}{lcccc}" << endl;
	cout << "Parameter	&	Mean	&	Std. Dev 	&	Pull mean	&	Pull Std. Dev	\\\\" << endl;
	for(UInt_t i =0; i< param.size(); i++){
		TString parname = param[i];
		parname.ReplaceAll("_","\\_");

		TString alertstr = " { ";
		TString strmean = prettyPrint(mean[i], emean[i]);
		TString strsigma = prettyPrint(sigma[i], esigma[i]);


		bool warnpull = (fabs(pull[i])/epull[i] > 3.0);
		TString strpull;
		if(warnpull){
			UInt_t wpullcol	= 0;
			strpull = "\\textcolor{";
			if (pull[i]<0.0){ strpull+= "blue!";}else{strpull+="red!";}
			wpullcol = (UInt_t)(fabs(200*pull[i]));
			if(wpullcol>100){wpullcol=100;}
			strpull+= wpullcol;
			if(wpullcol!=100){
				strpull+= "!black!";
				strpull+= 100;
			}
				strpull+= "}{";

		}
		strpull += prettyPrint(pull[i],epull[i]);
		if(warnpull) strpull += " } ";

		bool warnpullw = (fabs(pullw[i] -1)/epullw[i] > 3.0);
		TString strpullw;

		if(warnpullw){
		UInt_t wpullcolw = 0;
		strpullw = "\\textcolor{";
		if (pullw[i]<1.0){ 
			strpullw+= "blue!";
			wpullcolw = (UInt_t)((1.0 - pullw[i])*200);
		}else{
			strpullw+="red!";
			wpullcolw = (UInt_t)(fabs(1.0 - pullw[i])*200);
		}
		if(wpullcolw>100){wpullcolw=100;}
			strpullw+= wpullcolw;
			if(wpullcolw!=100){
				strpullw+= "!black!";
				strpullw+= 100;
			}
				strpullw+= "}{";
		}
		strpullw += prettyPrint(pullw[i],epullw[i]);
		if(warnpullw) strpullw += " } ";

		cout << parname << " 	& 	" << strmean << " 	& 	" << strsigma << "	 & 	" << strpull << "	 & 	" << strpullw << "	\\\\" << endl;

		//		cout << param[i] << " & 	$ " << prettyPrint(mean[i]) << " \\pm " << prettyPrint(emean[i]) << " $	&	$ " << prettyPrint(sigma[i]) << " \\pm " << prettyPrint(esigma[i]) << " $   &       $ " << prettyPrint(pull[i]) << " \\pm " << prettyPrint(epull[i]) <<  " $   &       $ " << prettyPrint(pullw[i]) << " \\pm " << prettyPrint(epullw[i]) << " $	\\\\" << endl;
	}
	cout << "\\end{tabular}" << endl;
	cout << "\\end{table}" << endl;
}
