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
#include "RooCBShape.h"
#include "RooProdPdf.h"
using std::cout;
using std::endl;
using namespace RooFit;
using namespace RooStats;

Width_t lwidth = 2;
double msize = 1;
int ncpu = 3;

void MakeFixed(RooArgList * params, bool constant){
	RooRealVar * arg;
	for(int i = 0; i< params->getSize(); i++){
		arg = (RooRealVar*)params->at(i);
		arg->setConstant(constant);
	}
}

void MakeFixed(RooArgSet * params, bool constant){
	MakeFixed(new RooArgList(*params), constant);
}

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

RooPlot *MakePlot(RooAbsPdf * model, RooAbsData * data, RooArgList* pdfs, RooRealVar *arg, RooArgSet *params, TString name){
	RooPlot *varplot = arg->frame();
	varplot->SetTitle("");
	data->plotOn(varplot,LineWidth(lwidth),LineColor(kBlack),LineWidth(lwidth));
	model->plotOn(varplot,LineWidth(lwidth),LineColor(kBlack),LineWidth(lwidth));
	Double_t chi = varplot->chiSquare();
	for(int i = 0; i<pdfs->getSize();i++){
		model->plotOn(varplot,Components(RooArgList(*pdfs->at(i))),LineColor(i+2),LineStyle(kDashed),LineWidth(lwidth));
	}
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

	TString jpsimassString = "mdau1";
	TString jpsimassUnit = "MeV/c^{2}";
	TString jpsimassTitle = "M(#mu #mu)";
	Double_t jpsimassMin = 3030.0;
	Double_t jpsimassMax = 3160.0;
	Double_t jpsimassMid = 3075.0;
	UInt_t jpsimassBins = 50;


	// Try to parse the input arguments
	if(argc !=5 && argc!=4){
		cout << "Sytnax: " << argv[0] <<" <signal.root> <signal path> <output.root> <n Bs Gaus>" << endl;
		return EXIT_FAILURE;
	}

	bool twogaus = false;
	bool jpsi = false;
	if(argc==5){
		if(atoi(argv[4])>1){
			twogaus=true;
		}
		if(atoi(argv[4])>2){
			jpsi=true;
		}
	}
	// Open input file and get ntuple
	TFile *inFile = new TFile(argv[1]);
	TTree * inTuple = (TTree*)inFile->Get(argv[2]);

	// Declare mass variable to be extracted from ntuple to fit to

	RooRealVar * obs_bs_mass = new RooRealVar(bsmassString,bsmassTitle,bsmassMid,bsmassMin,bsmassMax,bsmassUnit);
	RooRealVar * obs_jpsi_mass = new RooRealVar(jpsimassString,jpsimassTitle,jpsimassMid,jpsimassMin,jpsimassMax,jpsimassUnit);
	obs_bs_mass->setBins(bsmassBins);
	obs_jpsi_mass->setBins(jpsimassBins);

	// Load the ntuple into a RooDataset, taking only the mass column of the ntuple (as this is much faster to fit to)
	RooArgList *fitArgs = new RooArgList(*obs_bs_mass,*obs_jpsi_mass);
	RooDataSet * fitData = new RooDataSet("rawfitData","raw input fitDataset",inTuple,*fitArgs);
	Double_t entries = fitData->numEntries();
	cout << " SAMPLE CONTAINS: " << entries << " EVENTS" << endl;

	// Also load the ntuple fully so that when we sWeight it we get all the ntuple columns, not just the mass (don't fit to this though)
	// This makes lots of warnings as a RooRealVar casts everything as a double even if it's a float/int/whatever. This is unfortunate but not problematic
	RooArgList *allArgs = new RooArgList(*obs_bs_mass,*obs_jpsi_mass);
	TObjArray *members = inTuple->GetListOfLeaves();
	RooRealVar *aVar;
	for(int i = 0; i < members->GetEntries(); i++){
		TString var = members->At(i)->GetName();
		if(!var.Contains("COV") && !var.Contains("ERR") && (!var.Contains(bsmassString)) && (!var.Contains(jpsimassString))){
			aVar = new RooRealVar(var,var,-1.0);
			allArgs->add(*aVar,false);
		}

	}
	RooDataSet * allData = new RooDataSet("fitfitData","fit input fitDataset",inTuple,*allArgs);

	//create an ouput file
	TString outputdir = argv[3];
	gSystem->mkdir( outputdir );
	TFile * outputFile = new TFile(outputdir+"/result.root","RECREATE");

	// Yields	
	RooRealVar *Nsig_jpsi = new RooRealVar("Nsig_jpsi","Nsig_jpsi",10.0,0.0,entries);
	RooRealVar *Nbkg_jpsi = new RooRealVar("Nbkg_jpsi","Nbkg_jpsi",10.0,0.0,entries);
	RooRealVar *Nsig_bs = new RooRealVar("Nsig_bs","Nsig_ns",10.0,0.0,entries);
	RooRealVar *Nbkg_bs = new RooRealVar("Nbkg_bs","Nbkg_ns",10.0,0.0,entries);

	RooRealVar *Nsig = new RooRealVar("Nsig","Nsig",10.0,0.0,entries);
	RooRealVar *Nprompt = new RooRealVar("Nprompt","Nprompt",10.0,0.0,entries);
	RooRealVar *Nnojpsi = new RooRealVar("Nnojpsi","Nnojpsi",10.0,0.0,entries);
	RooRealVar *Nbkg = new RooRealVar("Nbkg","Nbkg",10.0,0.0,entries);

	RooArgList *yields = new RooArgList();

	if(jpsi){

		yields->add(RooArgList(*Nsig,*Nprompt,*Nnojpsi,*Nbkg));
	}else{
		yields->add(RooArgList(*Nsig,*Nbkg));

	}

	//Bs FIT PARAMS

	// Signal
	RooArgSet *bs_sig_vars = new RooArgSet();
	RooRealVar *bs_sig_sigma1 = new RooRealVar("bs_sig_sigma1","bs_sig_sigma1",6.0,0.0,16.0, "MeV/c^{2}");
	RooRealVar *bs_sig_mean1 = new RooRealVar("bs_sig_mean","bs_sig_mean",5.3671e+03,5360,5370, "MeV/c^{2}");
	RooRealVar *bs_sig_sigma_r1 = new RooRealVar("bs_sig_sigma_r1","bs_sig_sigma_r1",2.14);
	RooFormulaVar *bs_sig_sigma2 = new RooFormulaVar("bs_sig_sigma2","bs_sig_sigma1*bs_sig_sigma_r1",RooArgList(*bs_sig_sigma1,*bs_sig_sigma_r1));
	//RooRealVar *bs_sig_sigma2 = new RooRealVar("bs_sig_sigma2","bs_sig_sigma2",5.0,0.1,8.0, "MeV/c^{2}");
	RooRealVar *bs_sig_f1 = new RooRealVar("bs_sig_f1","bs_sig_f1",0.83);
	bs_sig_vars->add(RooArgList(*bs_sig_sigma1,*bs_sig_mean1));
	if(twogaus){
		bs_sig_vars->add(RooArgList(*bs_sig_sigma_r1,*bs_sig_f1));
	}

	// Background
	RooArgSet *bs_bkg_vars = new RooArgSet();
	RooRealVar *bs_bkg_coeff = new RooRealVar("bs_bkg_coeff","bs_bkg_coeff",-5.9060e-04,-0.005,0.005);
	bs_bkg_vars->add(RooArgList(*bs_bkg_coeff));

	// Both
	RooArgSet *bs_vars = new RooArgSet(*bs_sig_vars,*bs_bkg_vars);


	// Jpsi FIT PARAMS

	// Signal
	RooArgSet *jpsi_sig_vars = new RooArgSet();
	RooRealVar *jpsi_sig_mean = new RooRealVar("jpsi_sig_mean","jpsi_sig_mean",3095,3000,3200, "MeV/c^{2}");
	RooRealVar *jpsi_sig_sigma1 = new RooRealVar("jpsi_sig_sigma1","jpsi_sig_sigma1",0.011, -100,100, "MeV/c^{2}");
	//RooRealVar *jpsi_sig_sigma2 = new RooRealVar("jpsi_sig_sigma2","jpsi_sig_sigma2",0.017, 0.0,100, "MeV/c^{2}");
	RooRealVar *jpsi_sig_aa = new RooRealVar("jpsi_sig_aa", "jpsi_sig_aa", 1.975);
	RooRealVar *jpsi_sig_ab = new RooRealVar("jpsi_sig_ab", "jpsi_sig_ab", -0.0011);
	RooRealVar *jpsi_sig_ac = new RooRealVar("jpsi_sig_ac", "jpsi_sig_ac", -0.00018);
	//	RooRealVar *jpsi_sig_na = new RooRealVar("jpsi_sig_na", "jpsi_sig_na", 1.039);
	//	RooRealVar *jpsi_sig_nb = new RooRealVar("jpsi_sig_nb", "jpsi_sig_nb", -0.041);
	//	RooRealVar *jpsi_sig_nc = new RooRealVar("jpsi_sig_nc", "jpsi_sig_nc", 0.00227);
	RooRealVar *jpsi_sig_sa = new RooRealVar("jpsi_sig_sa","jpsi_sig_sa",-0.0546);
	RooRealVar *jpsi_sig_sb = new RooRealVar("jpsi_sig_sb","jpsi_sig_sb",0.26);
	RooRealVar *jpsi_sig_n = new RooRealVar("jpsi_sig_n","jpsi_sig_n",1.0);

	//RooRealVar *jpsi_sig_f1 = new RooRealVar("jpsi_sig_f1","jpsi_sig_f1",0.6,0.0,1.0);

	//jpsi_sig_vars->add(RooArgList(*jpsi_sig_mean,*jpsi_sig_sigma1,*jpsi_sig_sigma2,*jpsi_sig_f1));
	jpsi_sig_vars->add(RooArgList(*jpsi_sig_mean,*jpsi_sig_sigma1));
	jpsi_sig_vars->add(RooArgList(*jpsi_sig_aa,*jpsi_sig_ab,*jpsi_sig_ac));
	jpsi_sig_vars->add(RooArgList(*jpsi_sig_sa,*jpsi_sig_sb,*jpsi_sig_n));
	//jpsi_sig_vars->add(RooArgList( *jpsi_sig_na,*jpsi_sig_nb,*jpsi_sig_nc,*jpsi_sig_n));
	// Background
	RooArgSet *jpsi_bkg_vars = new RooArgSet();
	RooRealVar *jpsi_bkg_coeff = new RooRealVar("jpsi_bkg_coeff","jpsi_bkg_coeff",-1.3931e-03,-0.01,0.01);
	jpsi_bkg_vars->add(RooArgList(*jpsi_bkg_coeff));

	//Both
	RooArgSet *jpsi_vars = new RooArgSet(*jpsi_sig_vars,*jpsi_bkg_vars);


	// Bs FIT PDFs

	RooAbsPdf *bs_sig1 = new RooGaussian("bs_sig1","bs_sig1",*obs_bs_mass,*bs_sig_mean1,*bs_sig_sigma1);
	RooAbsPdf *bs_sig2 = new RooGaussian("bs_sig2","bs_sig2",*obs_bs_mass,*bs_sig_mean1,*bs_sig_sigma2);
	RooAbsPdf *bs_sig = bs_sig1;
	if(twogaus){
		bs_sig = new RooAddPdf("bs_sig","bs_sig",*bs_sig1,*bs_sig2,*bs_sig_f1);
	}

	RooAbsPdf *bs_bkg = new RooExponential("bs_bkg","bs_bkg",*obs_bs_mass,*bs_bkg_coeff);


	//Bs MASS FIT
	RooAbsPdf *bspdf = new RooAddPdf("bspdf","bspdf",RooArgSet(*bs_sig,*bs_bkg),RooArgSet(*Nsig_bs,*Nbkg_bs));
	//Here we perform the fit:
	//Create the NLL var
	RooAbsReal * bsnll = bspdf->createNLL(*fitData,NumCPU(ncpu),Extended(kTRUE));
	//Create the minuit instance
	RooMinuit *bsm = new RooMinuit(*bsnll);
	//Call migrad to minimise it
	bsm->migrad();
	//Save the fit result
	RooFitResult *bsresult = bsm->save();
	bsresult->Print();

	//Plot the fit
	TCanvas* bsmasscanv = new TCanvas("bsmasscanv","bsmasscanv",1024,768);
	if(twogaus){
	MakePlot(bspdf, fitData, new RooArgList(*bs_sig1,*bs_sig2,*bs_bkg), obs_bs_mass, bs_vars, "B_{s} Mass")->Draw();
	}else{
	MakePlot(bspdf, fitData, new RooArgList(*bs_sig,*bs_bkg), obs_bs_mass, bs_vars, "B_{s} Mass")->Draw();
	}
	bsmasscanv->Write();
	bsmasscanv->Print(outputdir+"/bsmass.pdf");
	bsmasscanv->Close();

	//Perform likelihood scans of all floated params in the fit:
	RooArgList bsscanparams = bsresult->floatParsFinal();
	bsscanparams.printLatex(Format("NEAU",AutoPrecision(1)));
	RooRealVar *bsscanpar;
	for(int i = 0; i< bsscanparams.getSize(); i++){
		bsscanpar = (RooRealVar*)bsscanparams.at(i);
		TCanvas * bsscanvas = new TCanvas(bsscanpar->getTitle()+"_scan",bsscanpar->getTitle()+"_scan",1024,768);
		LikelihoodScan(bsnll, bsscanpar)->Draw();
		bsscanvas->Write();
		bsscanvas->Print(outputdir+"/"+bsscanpar->getTitle()+"_scan.pdf");
		bsscanvas->Close();
	}
	
	//Fix the non-yield components of the fit
	MakeFixed(bs_vars,true);


	//Jpsi FIT PDFs
	RooFormulaVar *jpsi_sig_a1 = new RooFormulaVar("jpsi_sig_a1", "jpsi_sig_a1",  "jpsi_sig_aa + jpsi_sig_ab*jpsi_sig_sigma1 + jpsi_sig_ac*jpsi_sig_sigma1*jpsi_sig_sigma1",  RooArgSet(*jpsi_sig_aa, *jpsi_sig_ab, *jpsi_sig_ac,*jpsi_sig_sigma1));
	//RooFormulaVar *jpsi_sig_a2 = new RooFormulaVar("jpsi_sig_a2", "jpsi_sig_a2",  "jpsi_sig_aa + jpsi_sig_ab*jpsi_sig_sigma2 + jpsi_sig_ac*jpsi_sig_sigma2*jpsi_sig_sigma2",  RooArgSet(*jpsi_sig_aa, *jpsi_sig_ab, *jpsi_sig_ac, *jpsi_sig_sigma2)); 
	RooFormulaVar *jpsi_sig_mean1 = new RooFormulaVar("jpsi_sig_mean1","jpsi_sig_mean1","jpsi_sig_mean+jpsi_sig_sigma1*jpsi_sig_sa+jpsi_sig_sigma1*jpsi_sig_sigma1*jpsi_sig_sb",RooArgSet(*jpsi_sig_mean,*jpsi_sig_sigma1,*jpsi_sig_sa,*jpsi_sig_sb));
	//RooFormulaVar *jpsi_sig_mean2 = new RooFormulaVar("jpsi_sig_mean2","jpsi_sig_mean2","jpsi_sig_mean+jpsi_sig_sigma2*jpsi_sig_sa+jpsi_sig_sigma2*jpsi_sig_sigma2*jpsi_sig_sb",RooArgSet(*jpsi_sig_mean,*jpsi_sig_sigma2,*jpsi_sig_sa,*jpsi_sig_sb));


	RooAbsPdf *jpsi_sig = new RooCBShape("jpsi_sig_cb1","jpsi_sig_cb1",*obs_jpsi_mass, *jpsi_sig_mean1,*jpsi_sig_sigma1,*jpsi_sig_a1,*jpsi_sig_n);
	//RooAbsPdf *jpsi_sig_cb1 = new RooCBShape("jpsi_sig_cb1","jpsi_sig_cb1",*obs_jpsi_mass, *jpsi_sig_mean1,*jpsi_sig_sigma1,*jpsi_sig_a1,*jpsi_sig_n);
	//RooAbsPdf *jpsi_sig_cb2 = new RooCBShape("jpsi_sig_cb2","jpsi_sig_cb2",*obs_jpsi_mass, *jpsi_sig_mean2,*jpsi_sig_sigma2,*jpsi_sig_a2,*jpsi_sig_n);

	//RooAbsPdf *jpsi_sig = new RooAddPdf("jpsi_sig","jpsi_sig",*jpsi_sig_cb1,*jpsi_sig_cb2,*jpsi_sig_f1);

	RooAbsPdf *jpsi_bkg = new RooExponential("jpsi_bkg","jpsi_bkg",*obs_jpsi_mass,*jpsi_bkg_coeff);

	if(jpsi){
	//Jpsi MASS FIT
	RooAbsPdf *jpsipdf = new RooAddPdf("jpsipdf","jpsipdf",RooArgSet(*jpsi_sig,*jpsi_bkg),RooArgSet(*Nsig_jpsi,*Nbkg_jpsi));
	//Here we perform the fit:
	//Create the NLL var
	RooAbsReal * jpsinll = jpsipdf->createNLL(*fitData,NumCPU(ncpu),Extended(kTRUE));
	//Create the minuit instance
	RooMinuit *jpsim = new RooMinuit(*jpsinll);
	//Call migrad to minimise it
	jpsim->migrad();
	//Save the fit result
	RooFitResult *jpsiresult = jpsim->save();
	jpsiresult->Print();

	//Plot the fit
	TCanvas* jpsimasscanv = new TCanvas("jpsimasscanv","jpsimasscanv",1024,768);
	MakePlot(jpsipdf, fitData, new RooArgList(*jpsi_sig,*jpsi_bkg), obs_jpsi_mass, jpsi_vars, "J/#psi Mass")->Draw();
	//MakePlot(jpsipdf, fitData, new RooArgList(*jpsi_sig_cb1,*jpsi_sig_cb2,*jpsi_bkg), obs_jpsi_mass, jpsi_vars, "J/#psi Mass")->Draw();
	jpsimasscanv->Write();
	jpsimasscanv->Print(outputdir+"/jpsimass.pdf");
	jpsimasscanv->Close();
	

	//Perform likelihood scans of all floated params in the fit:
	RooArgList jpsiscanparams = jpsiresult->floatParsFinal();
	jpsiscanparams.printLatex(Format("NEAU",AutoPrecision(1)));
	RooRealVar *jpsiscanpar;
	for(int i = 0; i< jpsiscanparams.getSize(); i++){
		jpsiscanpar = (RooRealVar*)jpsiscanparams.at(i);
		TCanvas * jpsiscanvas = new TCanvas(jpsiscanpar->getTitle()+"_scan",jpsiscanpar->getTitle()+"_scan",1024,768);
		LikelihoodScan(jpsinll, jpsiscanpar)->Draw();
		jpsiscanvas->Write();
		jpsiscanvas->Print(outputdir+"/"+jpsiscanpar->getTitle()+"_scan.pdf");
		jpsiscanvas->Close();
	}

	//Fix the non-yield components of the fit
	MakeFixed(jpsi_vars,true);
	}
	//PRODPDFs
	RooArgList* pdfs = new RooArgList();
	RooAbsPdf *bs_sig_jpsi_sig = new RooProdPdf("bs_sig_jpsi_sig","bs_sig_jpsi_sig",RooArgSet(*bs_sig,*jpsi_sig));
	RooAbsPdf *bs_bkg_jpsi_sig = new RooProdPdf("bs_bkg_jpsi_sig","bs_bkg_jpsi_sig",RooArgSet(*bs_bkg,*jpsi_sig));
	RooAbsPdf *bs_sig_jpsi_bkg = new RooProdPdf("bs_sig_jpsi_bkg","bs_sig_jpsi_bkg",RooArgSet(*bs_sig,*jpsi_bkg));
	RooAbsPdf *bs_bkg_jpsi_bkg = new RooProdPdf("bs_bkg_jpsi_bkg","bs_bkg_jpsi_bkg",RooArgSet(*bs_bkg,*jpsi_bkg));
	if(jpsi){
		pdfs->add(RooArgList(*bs_sig_jpsi_sig,*bs_bkg_jpsi_sig,*bs_sig_jpsi_bkg,*bs_bkg_jpsi_bkg));
	}else{
		pdfs->add(RooArgList(*bs_sig,*bs_bkg));
	}

	RooAbsPdf *pdf = new RooAddPdf("pdf","pdf",*pdfs,*yields);
	//Here we perform the fit:
	//Create the NLL var
	RooAbsReal * nll = pdf->createNLL(*fitData,NumCPU(ncpu),Extended(kTRUE));
	//Create the minuit instance
	RooMinuit *m = new RooMinuit(*nll);
	//Call migrad to minimise it
	m->migrad();
	//Save the fit result
	RooFitResult *result = m->save();
	result->Print();


	//Plot the fit
	TCanvas* bsprodcanv = new TCanvas("bsprodcanv","bsprodcanv",1024,768);
	MakePlot(pdf, fitData, pdfs, obs_bs_mass, new RooArgSet(*yields), "B_{s} Mass")->Draw();
	bsprodcanv->Write();
	bsprodcanv->Print(outputdir+"/bsprod.pdf");
	bsprodcanv->Close();

	if(jpsi){
		TCanvas* jpsiprodcanv = new TCanvas("jpsiprodcanv","jpsiprodcanv",1024,768);
		MakePlot(pdf, fitData, pdfs, obs_jpsi_mass, new RooArgSet(*yields), "J/#psi Mass")->Draw();
		jpsiprodcanv->Write();
		jpsiprodcanv->Print(outputdir+"/jpsiprod.pdf");
		jpsiprodcanv->Close();
	}

	//Perform likelihood scans of all floated params in the fit:
	RooArgList scanparams = result->floatParsFinal();
	scanparams.printLatex(Format("NEAU",AutoPrecision(1)));
	RooRealVar *scanpar;
	for(int i = 0; i< scanparams.getSize(); i++){
		scanpar = (RooRealVar*)scanparams.at(i);
		TCanvas * scanvas = new TCanvas(scanpar->getTitle()+"_scan",scanpar->getTitle()+"_scan",1024,768);
		LikelihoodScan(nll, scanpar)->Draw();
		scanvas->Write();
		scanvas->Print(outputdir+"/"+scanpar->getTitle()+"_scan.pdf");
		scanvas->Close();
	}


	//Now the important bit: Make the sWeights and write them to the output file:
	RooStats::SPlot *splot = new RooStats::SPlot("sData","An SPlot", *allData, pdf, *yields);
	((RooTreeDataStore*)allData->store()->tree())->Write();
	delete splot;
	// Close the output file
	outputFile->Write();
	outputFile->Close();

}
