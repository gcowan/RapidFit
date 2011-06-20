#include "RooStats/SPlot.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooTreeDataStore.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooMinuit.h"
using std::cout;
using std::endl;
using namespace RooFit;
using namespace RooStats;

Width_t lwidth = 2;
double msize = 1;
int ncpu = 3;


// Macro to add pretty LHCb label to rooplots
void addLHCbLabel(RooPlot * frame1, TString legend){
	TPaveText * label = new TPaveText(0.14, 0.88, 0.14, 0.73,"BRNDC");
	label->SetFillColor(0);
	label->SetBorderSize(0);
	label->SetTextAlign(11);
	label->SetTextSize((Float_t)0.04);
	TText * labeltext = 0;
	labeltext = label->AddText("LHC#font[12]{b} 2011 Data");
	labeltext = label->AddText("#sqrt{s} = 7TeV");
	labeltext = label->AddText(legend);
	frame1->addObject(label);
}

// Macro to add parameters to rooplots (less ugly than the RooFit one)
void paramOn(RooArgSet * params, RooPlot * frame1, char* label, int numEntries, double chidof)
{
	double xmax = 0.14;
	double xmin = 0.14;
	double ymax = 0.73;
	double dy = 0.040;
	double ymin = ymax -dy;
	RooRealVar * var = 0;
	TIterator* pIter = params->createIterator() ;
	while((var=(RooRealVar*)pIter->Next())) {
		if(!var->isConstant()) ymin-= dy; 
	}
	if(string(label) != "") ymin-= dy; 
	TPaveText * box = new TPaveText(xmin, ymax, xmax, ymin,"BRNDC");
	box->SetFillColor(0);
	box->SetBorderSize(0);
	box->SetTextAlign(11);
	box->SetTextSize((Float_t)0.025);
	TText * text = 0;
	if(string(label) != "") text = box->AddText(label);
	char formatter [50];
	int sigfigs = 3;                
	RooLinkedList cmdList;    
	const RooCmdArg* formatCmd = static_cast<RooCmdArg*>(cmdList.FindObject("FormatArgs"));    
	pIter->Reset();    
	while((var=(RooRealVar*)pIter->Next())) {    
		if ( string(var->GetName()) == "#mu" ) sigfigs = 5;    
		TString *formatted = "NELU" ? var->format(sigfigs, "NELU") : var->format(*formatCmd) ;    
		text = box->AddText(formatted->Data());    
		delete formatted;    
	}                      
	sprintf(formatter, "Entries = %d", numEntries);
	if(numEntries > 0) text = box->AddText(formatter);
	if(chidof > 0.0){ 
		sprintf(formatter, "#chi^{2}/NDOF =  %4.3f", chidof);    
		text = box->AddText(formatter);                                                                                                    }    
		frame1->addObject(box);                                                                            
}

// Perform a fixed and profile likelihood in 1D for a given observable: 
RooPlot *LikelihoodScan(RooAbsReal * nll, RooRealVar *arg){
	TString legend = arg->getTitle() + " Likelihood Scan";
	RooUniformBinning binning = (RooUniformBinning&)arg->getBinning();
	RooPlot *scanplot = arg->frame(Bins(20));
	scanplot->SetTitle("");
	nll->plotOn(scanplot,ShiftToZero(),LineWidth(lwidth));
	RooCurve *curve = scanplot->getCurve();
	scanplot->SetMinimum(-10.0);
	scanplot->SetMaximum(curve->getYAxisMax());
	RooAbsReal* pll = nll->createProfile(*arg) ;
	pll->plotOn(scanplot,LineColor(kRed), LineWidth(lwidth)) ;
	addLHCbLabel(scanplot, legend);
	return scanplot;
}

// Make a nice plot 
RooPlot *MakePlot(RooAbsPdf * model, RooAbsData * data, RooAbsPdf * sig, RooAbsPdf *bkg, RooRealVar *arg, RooArgSet *params, TString name){
	RooPlot *varplot = arg->frame();
	varplot->SetTitle("");
	data->plotOn(varplot,LineWidth(lwidth),LineColor(kBlack),LineWidth(lwidth));
	model->plotOn(varplot,LineWidth(lwidth),LineColor(kBlack),LineWidth(lwidth));
	Double_t chi = varplot->chiSquare();
	model->plotOn(varplot,Components(RooArgSet(*sig)),LineColor(kBlue),LineStyle(kDashed),LineWidth(lwidth));
	model->plotOn(varplot,Components(RooArgSet(*bkg)),LineColor(kRed),LineStyle(kDashed),LineWidth(lwidth));
	addLHCbLabel(varplot, name);
	paramOn(params, varplot, (char*)"", data->numEntries(), chi);
	return varplot;
}


int main(int argc, char *argv[]){
	
	// Cosmetics: Make RooFit look good:
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
	gStyle->SetStatX(999.);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLineStyleString(2,"[12 12]");
	
	// Definitions of the mass parameter in the ntuple, ranges and binning to use:
	TString bsmassString = "mass";
	TString bsmassUnit = "MeV/c^{2}";
	TString bsmassTitle = "M(J/#psi #phi)";
	Double_t bsmassMin = 5200.0;
	Double_t bsmassMax = 5550.0;
	Double_t bsmassMid = 5367.0;
	UInt_t bsmassBins = 50;

	// Try to parse the input arguments
	if(argc !=4){
		cout << "Sytnax: " << argv[0] <<" <signal.root> <signal path> <output.root>" << endl;
		return EXIT_FAILURE;
	}

	// Open input file and get ntuple
	TFile *inFile = new TFile(argv[1]);
	TTree * inTuple = (TTree*)inFile->Get(argv[2]);

	// Declare mass variable to be extracted from ntuple to fit to
	RooRealVar * obs_bs_mass = new RooRealVar(bsmassString,bsmassTitle,bsmassMid,bsmassMin,bsmassMax,bsmassUnit);
	obs_bs_mass->setBins(bsmassBins);
	// Load the ntuple into a RooDataset, taking only the mass column of the ntuple (as this is much faster to fit to)
	RooArgList *fitArgs = new RooArgList(*obs_bs_mass);
	RooDataSet * fitData = new RooDataSet("rawfitData","raw input fitDataset",inTuple,*fitArgs);
	Double_t entries = fitData->numEntries();
	cout << " SAMPLE CONTAINS: " << entries << " EVENTS" << endl;

	// Also load the ntuple fully so that when we sWeight it we get all the ntuple columns, not just the mass (don't fit to this though)
	// This makes lots of warnings as a RooRealVar casts everything as a double even if it's a float/int/whatever. This is unfortunate but not problematic
	RooArgList *allArgs = new RooArgList(*obs_bs_mass);
	TObjArray *members = inTuple->GetListOfLeaves();
	RooRealVar *aVar;
	for(int i = 0; i < members->GetEntries(); i++){
		TString var = members->At(i)->GetName();
		if(!var.Contains("COV") && !var.Contains("ERR") && (!var.Contains(bsmassString))){
			aVar = new RooRealVar(var,var,-1.0);
			allArgs->add(*aVar,false);
		}

	}
	RooDataSet * allData = new RooDataSet("fitfitData","fit input fitDataset",inTuple,*allArgs);



	//create an ouput file
	TString outputdir = argv[3];
	gSystem->mkdir( outputdir );
	TFile * outputFile = new TFile(outputdir+"/result.root","RECREATE");

	//Bs FIT PARAMS

	// Signal
	RooArgSet *bs_sig_vars = new RooArgSet();
	RooRealVar *bs_sig_sigma1 = new RooRealVar("bs_sig_sigma","bs_sig_sigma",6.4016,1,20, "MeV/c^{2}");
	RooRealVar *bs_sig_mean1 = new RooRealVar("bs_sig_mean","bs_sig_mean",5.3671e+03,5350,5400, "MeV/c^{2}");
	bs_sig_vars->add(RooArgList(*bs_sig_sigma1,*bs_sig_mean1));

	// Background
	RooArgSet *bs_bkg_vars = new RooArgSet();
	RooRealVar *bs_bkg_coeff = new RooRealVar("bs_bkg_coeff","bs_bkg_coeff",-5.9060e-04,-0.01,0.01);
	bs_bkg_vars->add(RooArgList(*bs_bkg_coeff));

	// Both
	RooArgSet *bs_vars = new RooArgSet(*bs_sig_vars,*bs_bkg_vars);

	// Yields	
	RooRealVar *Nsig_bs = new RooRealVar("Nsig","Nsig",10.0,0.0,entries);
	RooRealVar *Nbkg_bs = new RooRealVar("Nbkg","Nbkg",10.0,0.0,entries);
	RooArgSet *bs_yields = new RooArgSet(*Nsig_bs,*Nbkg_bs);
	bs_vars->add(*bs_yields);

	// Bs FIT PDF
	RooAbsPdf *bs_sig = new RooGaussian("bs_sig","bs_sig",*obs_bs_mass,*bs_sig_mean1,*bs_sig_sigma1);
	RooAbsPdf *bs_bkg = new RooExponential("bs_bkg","bs_bkg",*obs_bs_mass,*bs_bkg_coeff);
	RooAbsPdf *bs = new RooAddPdf("bs","bs",RooArgList(*bs_sig,*bs_bkg),RooArgList(*Nsig_bs,*Nbkg_bs));

	//Here we perform the fit:
	//Create the NLL var
	RooAbsReal * bsnll = bs->createNLL(*fitData,NumCPU(ncpu),Extended(kTRUE));
	//Create the minuit instance
	RooMinuit *bsm = new RooMinuit(*bsnll);
	//Call migrad to minimise it
	bsm->migrad();
	//Save the fit result
	RooFitResult *bsresult = bsm->save();
	bsresult->Print();

	//Plot the fit
	TCanvas* masscanv = new TCanvas("masscanv","masscanv",1024,768);
	MakePlot(bs, fitData, bs_sig, bs_bkg, obs_bs_mass, bs_vars, "B_{s} Mass")->Draw();
	masscanv->Write();
	masscanv->Print(outputdir+"/mass.pdf");
	masscanv->Close();

	//Perform likelihood scans of all floated params in the fit:
	RooArgList scanparams = bsresult->floatParsFinal();
	scanparams.printLatex(Format("NEAU",AutoPrecision(1)));
	RooRealVar *scanpar;
	for(int i = 0; i< scanparams.getSize(); i++){
		scanpar = (RooRealVar*)scanparams.at(i);
		TCanvas * scanvas = new TCanvas(scanpar->getTitle()+"_scan",scanpar->getTitle()+"_scan",1024,768);
		LikelihoodScan(bsnll, scanpar)->Draw();
		scanvas->Write();
		scanvas->Print(outputdir+"/"+scanpar->getTitle()+"_scan.pdf");
		scanvas->Close();
	}

	//Now the important bit: Make the sWeights and write them to the output file:
	RooStats::SPlot *splot = new RooStats::SPlot("sData","An SPlot", *allData, bs, *bs_yields);
	((RooTreeDataStore*)allData->store()->tree())->Write();
	delete splot;
	// Close the output file
	outputFile->Write();
	outputFile->Close();
}
