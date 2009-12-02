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

int main(int argc, char *argv[]){

	if(argc !=3){
		cout << "Binary, slighlty easier version of Grieg's python script. Usage:" << endl;
		cout << "rapidresults PullPlots.root resultdir"	<< endl;
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
	TIter nextkey( gDirectory->GetListOfKeys() );
	TKey *key;
	TString outputdir = argv[2];
	gSystem->mkdir( outputdir );
	TFile * output = new TFile(outputdir+"/results.root","RECREATE");
	input->cd();
	while ( (key = (TKey*)nextkey())) {
		TObject *obj = key->ReadObj();
		if ( obj->IsA()->InheritsFrom( "TTree" ) ) {
			TTree *tree = (TTree*)obj;
			TString varname = tree->GetName();
			if(varname != "fitInfo"){
				param.push_back(varname);
				TObjArray *leaves = tree->GetListOfLeaves();
				Int_t nLeaves = leaves->GetEntries();
				for(int i = 0; i < nLeaves; i++){
					TString subVarname = leaves->At(i)->GetName();
					tree->Draw(subVarname+">>tmp");
					TH1F *hist = (TH1F*)gDirectory->Get("tmp");
					hist->SetName(varname+"_"+subVarname);
					hist->SetTitle(varname+" "+subVarname);
					hist->Sumw2();
					TCanvas *c = new TCanvas(varname+"_"+subVarname,varname+"_"+subVarname,2048,1536);
					hist->Draw();
					output->cd();
					hist->Write();
					c->Update();
					if(subVarname != "error"){
						TF1 *gauss = new TF1("Gaussian_"+subVarname,"gaus",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
						hist->Fit(gauss);
						gauss->Draw("SAME");
						if(subVarname == "value"){
							mean.push_back(gauss->GetParameter(1));
							emean.push_back(gauss->GetParError(1));
							sigma.push_back(gauss->GetParameter(2));
							esigma.push_back(gauss->GetParError(2));
						}
						if(subVarname == "pull"){
							pull.push_back(gauss->GetParameter(1));
							epull.push_back(gauss->GetParError(1));
							pullw.push_back(gauss->GetParameter(2));
							epullw.push_back(gauss->GetParError(2));
						}
						c->Update();
						gauss->Delete();
					}
					c->Print(outputdir+"/"+varname+"_"+subVarname+".pdf");
					c->Write();
					hist->Delete();
					input->cd();
				}
			}
		}
	}
	output->Write();
	output->Close();

	cout << "\\begin{center}" << endl;
	cout << "\\begin{tabular}{|l|c|c|c|c|}" << endl;

	cout << "\\hline" << endl;
	cout << "Parameter	&	Mean	&	Std. Dev 	&	Pull mean	&	Pull Std. Dev	\\\\" << endl;

	cout << "\\hline" << endl;
	for(UInt_t i =0; i< param.size(); i++){
		cout << param[i] << " & 	$ " << mean[i] << " \\pm " << emean[i] << " $	&	$ " << sigma[i] << " \\pm " << esigma[i] << " $   &       $ " << pull[i] << " \\pm " << epull[i] <<  " $   &       $ " << pullw[i] << " \\pm " << epullw[i] << " $	\\\\" << endl;
	}
	
	cout << "\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "\\end{center}" << endl;
}
