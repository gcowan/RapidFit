//	ROOT Headers
#include "TRandom.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TMatrixD.h"
#include "TDirectory.h"
#include "TSystem.h"
//	RapidFit Headers
#include "PerEventAngularAcceptance.h"
//	System Headers
#include <typeinfo>
#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

//Default constructor
PerEventAngularAcceptance::PerEventAngularAcceptance() :
	//	Did the author intend for ALL of these to be persistant in the heap?
	tupleFile(), decaytree(), evtList(), nev(), particles(), resonances(), partsAndRes(),
	p(), pInB(), pInJpsi(), cosThetaMuHistos(), cosThetaMuHistosPredicted(), cosThetaMuTmpHistos(),
	cosThetaMuSumHistos(), cosThetaMuHistosEffRatio(), cosThetaKHistos(), cosThetaKHistosPredicted(),
	cosThetaKTmpHistos(), cosThetaKSumHistos(), cosThetaKHistosEffRatio(), pHistos(), pHistosPredicted(),
	pTmpHistos(), pSumHistos(), pHistosEffRatio(), pHistosEffRatioPrev(), ptHistos(), ptHistosPredicted(),
	ptTmpHistos(), ptSumHistos(), ptHistosEffRatio(), pptHistos(), pptHistosPredicted(), pptTmpHistos(),
	pptSumHistos(), pptHistosEffRatio(), boosttoJpsi(), boostBackToLabFromJpsi(), boosttoB(),
	boostBackToLabFromB(), GeV(), pmin(), pmax(), pbinWidth(), ptmin(), ptmax(), ptbinWidth(), amin(), amax(),
	aymax(), pmin_jpsi(), ptmin_jpsi(), cosbins(), pbins(), ptbins(), w_sum_gen(), muonplus_PX(), muonplus_PY(),
	muonplus_PZ(), muonplus_PE(), muonminus_PX(), muonminus_PY(), muonminus_PZ(), muonminus_PE(), kaon_PX(),
	kaon_PY(), kaon_PZ(), kaon_PE(), Bu_MM(), Bu_TAU(), Bu_BKGCAT(), b_Bu_BKGCAT(), b_Bu_MM(), b_Bu_TAU(),
	b_muonplus_PX(), b_muonplus_PY(), b_muonplus_PZ(), b_muonplus_PE(),b_muonminus_PX(), b_muonminus_PY(),
	b_muonminus_PZ(), b_muonminus_PE(), b_kaon_PX(), b_kaon_PY(), b_kaon_PZ(), b_kaon_PE()
{
}

//Destructor
PerEventAngularAcceptance::~PerEventAngularAcceptance()
{
}

//Constructor with correct argument
PerEventAngularAcceptance::PerEventAngularAcceptance( string fileName, string tupleName, string outFileName ) : 
	//	Did the author intend for ALL of these to be persistant in the heap?
	tupleFile(), decaytree(), evtList(), nev(), particles(), resonances(), partsAndRes(),
	p(), pInB(), pInJpsi(), cosThetaMuHistos(), cosThetaMuHistosPredicted(), cosThetaMuTmpHistos(),
	cosThetaMuSumHistos(), cosThetaMuHistosEffRatio(), cosThetaKHistos(), cosThetaKHistosPredicted(),
	cosThetaKTmpHistos(), cosThetaKSumHistos(), cosThetaKHistosEffRatio(), pHistos(), pHistosPredicted(),
	pTmpHistos(), pSumHistos(), pHistosEffRatio(), pHistosEffRatioPrev(), ptHistos(), ptHistosPredicted(),
	ptTmpHistos(), ptSumHistos(), ptHistosEffRatio(), pptHistos(), pptHistosPredicted(), pptTmpHistos(),
	pptSumHistos(), pptHistosEffRatio(), boosttoJpsi(), boostBackToLabFromJpsi(), boosttoB(),
	boostBackToLabFromB(), GeV(), pmin(), pmax(), pbinWidth(), ptmin(), ptmax(), ptbinWidth(), amin(), amax(),
	aymax(), pmin_jpsi(), ptmin_jpsi(), cosbins(), pbins(), ptbins(), w_sum_gen(), muonplus_PX(), muonplus_PY(),
	muonplus_PZ(), muonplus_PE(), muonminus_PX(), muonminus_PY(), muonminus_PZ(), muonminus_PE(), kaon_PX(),
	kaon_PY(), kaon_PZ(), kaon_PE(), Bu_MM(), Bu_TAU(), Bu_BKGCAT(), b_Bu_BKGCAT(), b_Bu_MM(), b_Bu_TAU(),
	b_muonplus_PX(), b_muonplus_PY(), b_muonplus_PZ(), b_muonplus_PE(),b_muonminus_PX(), b_muonminus_PY(),
	b_muonminus_PZ(), b_muonminus_PE(), b_kaon_PX(), b_kaon_PY(), b_kaon_PZ(), b_kaon_PE()
{
	string null_s = outFileName; null_s = "";
	setupNtuple( fileName, tupleName );

	gRandom->SetSeed(1);

	GeV = 1000.;
	pbins = 40;
	pmin = 3.*GeV;
	pmax = 100.*GeV;
	pbinWidth = (pmax - pmin)/pbins;
	ptbins = 10;
	ptmin = 0.5*GeV;
	ptmax = 5.*GeV;
	ptbinWidth = (ptmax - ptmin)/ptbins;
	cosbins = 50;

	pmin_jpsi = 0.*GeV;
	ptmin_jpsi = 1.*GeV;

	amin = 0.02;
	amax = 0.28;
	aymax = 0.25;

	particles.push_back("kaon");
	particles.push_back("muonplus");
	particles.push_back("muonminus");

	resonances.push_back("Jpsi");
	resonances.push_back("B");

	partsAndRes.push_back("kaon");
	partsAndRes.push_back("muonplus");
	partsAndRes.push_back("muonminus");
	partsAndRes.push_back("Jpsi");
	partsAndRes.push_back("B");

	TVector3 boosttoJpsi(0,0,0), boostBackToLabFromJpsi(0,0,0);
	TVector3 boosttoB(0,0,0), boostBackToLabFromB(0,0,0);

	vector<string>::iterator particleIterator;
	for ( particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator )
	{
		string name = *particleIterator;
		cosThetaMuHistos[name]         = new TH1F((name+"_cosMu").c_str(),      name.c_str(), cosbins, -1., 1);
		cosThetaMuTmpHistos[name]      = new TH1F((name+"_cosMu_tmp").c_str(),  name.c_str(), cosbins, -1., 1);
		cosThetaMuSumHistos[name]      = new TH1F((name+"_cosMu_sum").c_str(),  name.c_str(), cosbins,  -1., 1);
		cosThetaMuHistosPredicted[name]= new TH1F((name+"_cosMu_pre").c_str(),  name.c_str(), cosbins,  -1., 1);
		cosThetaMuHistosEffRatio[name] = new TH1F((name+"_cosMu_eff").c_str(),  name.c_str(), cosbins,  -1., 1);
		cosThetaKHistos[name]          = new TH1F((name+"_cosK").c_str(),     name.c_str(), cosbins,  -1., 1);
		cosThetaKTmpHistos[name]       = new TH1F((name+"_cosK_tmp").c_str(), name.c_str(), cosbins,  -1., 1);
		cosThetaKSumHistos[name]       = new TH1F((name+"_cosK_sum").c_str(), name.c_str(), cosbins,  -1., 1);
		cosThetaKHistosPredicted[name] = new TH1F((name+"_cosK_pre").c_str(), name.c_str(), cosbins,  -1., 1);
		cosThetaKHistosEffRatio[name]  = new TH1F((name+"_cosK_eff").c_str(), name.c_str(), cosbins,  -1., 1);
		pHistos[name]              = new TH1F((name+"_p").c_str(),      name.c_str(), pbins, pmin, pmax);
		pTmpHistos[name]           = new TH1F((name+"_p_tmp").c_str(),  name.c_str(), pbins, pmin, pmax);
		pSumHistos[name]           = new TH1F((name+"_p_sum").c_str(),  name.c_str(), pbins, pmin, pmax);
		pHistosPredicted[name]     = new TH1F((name+"_p_pre").c_str(),  name.c_str(), pbins, pmin, pmax);
		pHistosEffRatio[name]      = new TH1F((name+"_p_eff").c_str(),  name.c_str(), pbins, pmin, pmax);
		pHistosEffRatioPrev[name]  = new TH1F((name+"_p_eff_prev").c_str(),  name.c_str(), pbins, pmin, pmax);
		ptHistos[name]             = new TH1F((name+"_pt").c_str(),     name.c_str(), ptbins, ptmin, ptmax);
		ptTmpHistos[name]          = new TH1F((name+"_pt_tmp").c_str(), name.c_str(), ptbins, ptmin, ptmax);
		ptSumHistos[name]          = new TH1F((name+"_pt_sum").c_str(), name.c_str(), ptbins, ptmin, ptmax);
		ptHistosPredicted[name]    = new TH1F((name+"_pt_pre").c_str(), name.c_str(), ptbins, ptmin, ptmax);
		ptHistosEffRatio[name]     = new TH1F((name+"_pt_eff").c_str(), name.c_str(), ptbins, ptmin, ptmax);
		pptHistos[name]            = new TH2F((name+"_p_pt").c_str(),     name.c_str(), pbins, pmin, pmax, ptbins, ptmin, ptmax);
		pptTmpHistos[name]         = new TH2F((name+"_p_pt_tmp").c_str(), name.c_str(), pbins, pmin, pmax, ptbins, ptmin, ptmax);
		pptSumHistos[name]         = new TH2F((name+"_p_pt_sum").c_str(), name.c_str(), pbins, pmin, pmax, ptbins, ptmin, ptmax);
		pptHistosPredicted[name]   = new TH2F((name+"_p_pt_pre").c_str(), name.c_str(), pbins, pmin, pmax, ptbins, ptmin, ptmax);
		pptHistosEffRatio[name]    = new TH2F((name+"_p_pt_eff").c_str(), name.c_str(), pbins, pmin, pmax, ptbins, ptmin, ptmax);
		cosThetaMuHistos[name]->Sumw2();
		cosThetaMuHistosPredicted[name]->Sumw2();
		cosThetaKHistos[name]->Sumw2();
		cosThetaKHistosPredicted[name]->Sumw2();
		pHistos[name]->Sumw2();
		pTmpHistos[name]->Sumw2();
		pSumHistos[name]->Sumw2();
		pHistosPredicted[name]->Sumw2();
		pHistosEffRatio[name]->Sumw2();
		pHistosEffRatioPrev[name]->Sumw2();
		pHistosEffRatioPrev[name]->SetBit(TH1::kIsAverage);
		ptHistos[name]->Sumw2();
		ptTmpHistos[name]->Sumw2();
		ptSumHistos[name]->Sumw2();
		ptHistosPredicted[name]->Sumw2();
		ptHistosEffRatio[name]->Sumw2();
		pptHistos[name]->Sumw2();
		pptTmpHistos[name]->Sumw2();
		pptSumHistos[name]->Sumw2();
		pptHistosPredicted[name]->Sumw2();
		pptHistosEffRatio[name]->Sumw2();
	}
}

void PerEventAngularAcceptance::setupNtuple( string tupleFileName, string tupleName )
{
	tupleFile = new TFile( tupleFileName.c_str(), "READ" );
	decaytree = (TNtuple*)tupleFile->Get( tupleName.c_str() );
	if ( decaytree == NULL )
	{
		cerr << "Ntuple not found. Are you specifying the correct file and ntuple path?" << endl;
		exit(1);
	}

	string cut = " Bu_M > 5200 && Bu_M < 5360";
	cut += " && Bu_TAU > -0.2 && Bu_TAU < 20";
	cut += " && muonplus_P > 3000. && muonplus_PT > 500.";
	cut += " && muonminus_P > 3000. && muonminus_PT > 500.";
	cut += " && kaon_P > 3000. && kaon_PT > 500.";
	cut += " && Jpsi_P > 0. && Jpsi_PT > 1000.";
	cut += " && abs(Jpsi_M-3096.) < 50.";
	cut += " && Bu_BKGCAT==0";

	string cut1 = "1==1";

	decaytree->Draw(">>evtList", cut1.c_str());
	evtList = (TEventList*)gDirectory->Get("evtList");
	nev = evtList->GetN();
	nev = 10000;
	cout << "Number of events passing cut: " << nev << endl;

	decaytree->SetMakeClass(1);
	decaytree->SetBranchStatus("*",1);
	//decaytree->SetBranchAddress("Bu_BKGCAT",   &B_BKGCAT, &b_B_BKGCAT);
	//decaytree->SetBranchAddress("Bu_MM",   &B_MM, &b_B_MM);
	//decaytree->SetBranchAddress("Bu_TAU",   &B_TAU, &b_B_TAU);

	decaytree->SetBranchAddress("muonplus_PX",   &muonplus_PX, &b_muonplus_PX);
	decaytree->SetBranchAddress("muonplus_PY",   &muonplus_PY, &b_muonplus_PY);
	decaytree->SetBranchAddress("muonplus_PZ",   &muonplus_PZ, &b_muonplus_PZ);
	decaytree->SetBranchAddress("muonplus_PE",   &muonplus_PE, &b_muonplus_PE);

	decaytree->SetBranchAddress("muonminus_PX",   &muonminus_PX, &b_muonminus_PX);
	decaytree->SetBranchAddress("muonminus_PY",   &muonminus_PY, &b_muonminus_PY);
	decaytree->SetBranchAddress("muonminus_PZ",   &muonminus_PZ, &b_muonminus_PZ);
	decaytree->SetBranchAddress("muonminus_PE",   &muonminus_PE, &b_muonminus_PE);

	decaytree->SetBranchAddress("kaon_PX",   &kaon_PX, &b_kaon_PX);
	decaytree->SetBranchAddress("kaon_PY",   &kaon_PY, &b_kaon_PY);
	decaytree->SetBranchAddress("kaon_PZ",   &kaon_PZ, &b_kaon_PZ);
	decaytree->SetBranchAddress("kaon_PE",   &kaon_PE, &b_kaon_PE);

	decaytree->SetBranchAddress("Bu_BKGCAT",   &Bu_BKGCAT, &b_Bu_BKGCAT);
	decaytree->SetBranchAddress("Bu_MM",   &Bu_MM, &b_Bu_MM);
	decaytree->SetBranchAddress("Bu_TAU",   &Bu_TAU, &b_Bu_TAU);
}

void PerEventAngularAcceptance::fillEffHistos( int iteration )
{
	//Fill the histograms. Start with a flat efficiency.
	cout << "Filling the efficiency histograms for iteration " << iteration << endl;
	//Need these since we will overwrite the existing histograms
	map<string, TH1F *> pTmp;
	map<string, TH1F *> ptTmp;
	map<string, TH2F *> pptTmp;

	vector<string>::iterator particleIterator;
	for ( particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator )
	{
		string particle = *particleIterator;

		pHistosEffRatioPrev[particle] = (TH1F*)pHistosEffRatio[particle]->Clone("pEffPrev");
		pTmp[particle] = (TH1F*)pHistosEffRatio[particle]->Clone("ptmp");
		ptTmp[particle] = (TH1F*)ptHistosEffRatio[particle]->Clone("pttmp");
		pptTmp[particle] = (TH2F*)pptHistosEffRatio[particle]->Clone("ppttmp");

		pHistosEffRatio[particle]->Reset();
		ptHistosEffRatio[particle]->Reset();
		pptHistosEffRatio[particle]->Reset();

		int pBinNumber = 1;
		int ptBinNumber = 1;
		double weight = 0.0;
		//First do the p histo
		for ( double pbin = pmin; pbin < pmax; pbin += pbinWidth )
		{
			weight = iteration == 1 ? 0.85 : pTmp[particle]->GetBinContent(pBinNumber);
			pHistosEffRatio[particle]->Fill(pbin, weight);
			// Now also do the 2D histogram
			ptBinNumber = 1;
			for ( double ptbin = ptmin; ptbin < ptmax; ptbin += ptbinWidth )
			{
				//cout << pbin << " " << ptbin << " " << pBinNumber << " " << ptBinNumber << endl;
				weight = iteration == 1 ? 0.85 : pptTmp[particle]->GetBinContent(pBinNumber, ptBinNumber);
				pptHistosEffRatio[particle]->Fill(pbin, ptbin, weight);
				ptBinNumber += 1;
			}
			pBinNumber += 1;
		}
		// Now do the pt one on it's own
		ptBinNumber = 1;
		for ( double ptbin = ptmin; ptbin < ptmax; ptbin += ptbinWidth )
		{
			weight = iteration == 1 ? 0.85 : ptTmp[particle]->GetBinContent(ptBinNumber);
			ptHistosEffRatio[particle]->Fill(ptbin, weight);
			ptBinNumber += 1;
		}
		if (iteration == 1) {
			pHistosEffRatioPrev[particle] = (TH1F*)pHistosEffRatio[particle]->Clone("pEffPrev");
		}
	}
}

void PerEventAngularAcceptance::loopOnReconstructedBs()
{
	vector<string>::iterator particleIterator;
	for ( particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator )
	{
		string particle = *particleIterator;
		cosThetaMuHistos[particle]->Reset();
		cosThetaKHistos[particle]->Reset();
		pHistos[particle]->Reset();
		ptHistos[particle]->Reset();
		pptHistos[particle]->Reset();
		pHistosPredicted[particle]->Reset();
		ptHistosPredicted[particle]->Reset();
		pptHistosPredicted[particle]->Reset();
	}

	for (int i = 0; i < nev; ++i)
	{
		if (i % 10000 == 0) cout << "Processing event " << i << endl;
		decaytree->GetEntry(evtList->GetEntry(i));

		p["kaon"]      = new TLorentzVector(kaon_PX, kaon_PY, kaon_PZ, kaon_PE);
		p["muonminus"] = new TLorentzVector(muonminus_PX, muonminus_PY, muonminus_PZ, muonminus_PE);
		p["muonplus"]  = new TLorentzVector(muonplus_PX, muonplus_PY, muonplus_PZ, muonplus_PE);

		// Jpsi and B momentum
		p["Jpsi"] = new TLorentzVector(*p["muonminus"] + *p["muonplus"]);
		p["B"] = new TLorentzVector(*p["Jpsi"] + *p["kaon"]);

		if(fabs(Bu_BKGCAT-0)<DOUBLE_TOLERANCE) continue; 
		if(Bu_MM < 5200. || Bu_MM > 5360.) continue;
		if(Bu_TAU < -0.2 || Bu_TAU > 20.) continue;

		if(fabs(p["Jpsi"]->M()-3096)>50) continue;

		if(p["muonplus"]->P() < pmin) continue;
		if(p["muonplus"]->Pt() < ptmin) continue;
		if(p["muonminus"]->P() < pmin) continue;
		if(p["muonminus"]->Pt() < ptmin) continue;
		if(p["kaon"]->P() < pmin) continue;
		if(p["kaon"]->Pt() < ptmin) continue;

		if(p["Jpsi"]->P() < pmin_jpsi) continue;
		if(p["Jpsi"]->Pt() < ptmin_jpsi) continue;

		// want the reconstructed B's to be inside the fiducial volume of LHCb
		vector<int> inacc = fiducialCuts( *p["muonplus"], *p["muonminus"], *p["kaon"], *p["Jpsi"] );

		if (inacc[0]+inacc[1]+inacc[2]+inacc[3] < 4) continue;
		/*cout << i << " " << Bu_BKGCAT << " " << Bu_MM << " "
		  << Bu_TAU << " "
		  << p["Jpsi"]->M()-3096 << Bu_TAU << " "
		  << p["muonplus"]->P() << " "
		  << p["muonplus"]->Pt() << " "
		  << p["muonminus"]->P() << " "
		  << p["muonminus"]->Pt() << " "
		  << p["kaon"]->P() << " "
		  << p["kaon"]->Pt() << " "
		  << p["Jpsi"]->P() << " "
		  << p["Jpsi"]->Pt() << " "

		  << endl;
		  */

		// Fill some histograms of the reconstructed quantities
		for ( particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator )
		{
			string particle = *particleIterator;
			pHistos[particle]->Fill(p[particle]->P());
			ptHistos[particle]->Fill(p[particle]->Pt());
			pptHistos[particle]->Fill(p[particle]->P(), p[particle]->Pt());
		}

		// calculate the reconstructed angles
		double csth_mu = coshel(*p["muonplus"], *p["Jpsi"], *p["B"]);
		double csth_k = cosk(*p["kaon"], *p["B"]);

		// Fill some histograms of these. Do it for each particle - this is redundant, but can't be bothered changing now
		for ( particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator )
		{
			string particle = *particleIterator;
			cosThetaMuHistos[particle]->Fill(csth_mu);
			cosThetaKHistos[particle]->Fill(csth_k);
		}

		// Now boost
		boost();

		// Now calculate f given p_B
		generateEventsAndCalculateFgivenB();

		//cout << w_sum_gen << endl;
		// Now do the sum of the generated distributions
		vector<string>::iterator particleIterator;
		for ( particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator )
		{
			string particle = *particleIterator;
			cosThetaMuHistosPredicted[particle]->Add(cosThetaMuTmpHistos[particle], 1/w_sum_gen);
			cosThetaKHistosPredicted[particle]->Add(cosThetaKTmpHistos[particle], 1/w_sum_gen);
			pHistosPredicted[particle]->Add(pTmpHistos[particle], 1/w_sum_gen);
			ptHistosPredicted[particle]->Add(ptTmpHistos[particle], 1/w_sum_gen);
			pptHistosPredicted[particle]->Add(pptTmpHistos[particle], 1/w_sum_gen);
			cosThetaMuTmpHistos[particle]->Reset();
			cosThetaKTmpHistos[particle]->Reset();
			pTmpHistos[particle]->Reset();
			ptTmpHistos[particle]->Reset();
			pptTmpHistos[particle]->Reset();
		}
	}
	// Now divide the predicted distributions by the actual reconstructed distribution to get the efficiency
	for ( particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator )
	{
		string particle = *particleIterator;
		//pHistosPredicted[particle]->Divide(pHistosEffRatioPrev[particle]);//, 1., 1., "b");
		pHistosEffRatio[particle]->Divide(pHistos[particle], pHistosPredicted[particle], 1., 1., "b");
		pHistosEffRatio[particle]->Multiply(pHistosEffRatioPrev[particle]);
		ptHistosEffRatio[particle]->Divide(ptHistos[particle], ptHistosPredicted[particle]);
		pptHistosEffRatio[particle]->Divide(pptHistos[particle], pptHistosPredicted[particle]);
		delete p[particle];
		delete pInJpsi[particle];
		delete pInB[particle];
	}
}

void PerEventAngularAcceptance::generateEventsAndCalculateFgivenB()
{
	//Generate events for each B momentum and then form the distribution
	//p, weighted by the product of the individual particle efficiencies,
	//giving f(p|p_B). What about the normalisation of f?
	vector<string>::iterator particleIterator;
	for ( particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator )
	{
		string particle = *particleIterator;
		cosThetaMuTmpHistos[particle]->Reset();
		cosThetaKTmpHistos[particle]->Reset();
		pTmpHistos[particle]->Reset();
		ptTmpHistos[particle]->Reset();
		pptTmpHistos[particle]->Reset();
	}

	map<string, TLorentzVector *> newP;

	w_sum_gen = 0.;

	while (w_sum_gen < 50)
	{
		double cosThetaK = gRandom->Uniform(-1., 1.);
		double phiK = gRandom->Uniform(0., 2.*TMath::Pi());
		// Yuehong says the K has an isotropic distribution in the rest frame of the B. So why this functional form?
		double x = pInB["kaon"]->P() * sqrt(1 - cosThetaK*cosThetaK)*cos(phiK);
		double y = pInB["kaon"]->P() * sqrt(1 - cosThetaK*cosThetaK)*sin(phiK);
		double z = pInB["kaon"]->P() * cosThetaK;
		double e = pInB["kaon"]->E();

		newP["kaon"] = new TLorentzVector(x, y, z, e );
		newP["kaon"]->Boost(boostBackToLabFromB);

		// Jpsi must have opposite momentum to the K in the B rest frame
		// If I use the form below, when I do the boose of p['B'] to the NewJpsi frame, I get NAN.
		newP["Jpsi"] = new TLorentzVector(-x, -y, -z, pInB["Jpsi"]->E());
		newP["Jpsi"]->Boost(boostBackToLabFromB);

		// Now boost this Jpsi (which is in the lab frame) to it's rest frame and work out the muon angular distributions
		TVector3 boosttoNewJpsi = -(newP["Jpsi"]->BoostVector());
		TVector3 boostBackToLabFromNewJpsi = newP["Jpsi"]->BoostVector();

		double cosThetaMu = gRandom->Uniform(-1., 1.);
		double phiMu = gRandom->Uniform(0., 2*TMath::Pi());
		double xprim = pInJpsi["muonplus"]->P() * sqrt(1 - cosThetaMu*cosThetaMu)*cos(phiMu);
		double yprim = pInJpsi["muonplus"]->P() * sqrt(1 - cosThetaMu*cosThetaMu)*sin(phiMu);
		double zprim = pInJpsi["muonplus"]->P() * cosThetaMu;

		//cout << cosThetaMu << " " << phiMu << " " << xprim <<  " " << yprim << " " << zprim << endl;

		// Now do the Jpsi. Which particles do we need here? Check!!!!!!!!
		TLorentzVector p4BInNewJpsi = *p["B"];
		p4BInNewJpsi.Boost(boosttoNewJpsi);

		//cout << p4BInNewJpsi.P() << endl;

		TVector3 uzInNewJpsi = -((p4BInNewJpsi.Vect()).Unit());
		TVector3 uzInLab(0,0,1);
		TVector3 uyInNewJpsi = (uzInNewJpsi.Cross(uzInLab) ).Unit();
		TVector3 uxInNewJpsi = (uyInNewJpsi.Cross(uzInNewJpsi) ).Unit();
		/*
		   cout << uzInNewJpsi.X() << endl;
		   cout << uzInLab.X() << endl;
		   cout << uyInNewJpsi.X() << endl;
		   cout << uxInNewJpsi.X() << endl;
		   */
		TVector3 v3newmup = xprim*uxInNewJpsi + yprim*uyInNewJpsi + zprim*uzInNewJpsi;
		TVector3 v3newmum = -xprim*uxInNewJpsi - yprim*uyInNewJpsi - zprim*uzInNewJpsi;

		//cout << v3newmup.X() << " " << v3newmup.Y() << " " << v3newmup.Z() << endl;

		newP["muonplus"]  = new TLorentzVector(v3newmup.Px(),v3newmup.Py(),v3newmup.Pz(), pInJpsi["muonplus"]->E());
		newP["muonminus"] = new TLorentzVector(v3newmum.Px(),v3newmum.Py(),v3newmum.Pz(), pInJpsi["muonminus"]->E());

		newP["muonplus"]->Boost(boostBackToLabFromNewJpsi);
		newP["muonminus"]->Boost(boostBackToLabFromNewJpsi);

		//cout << newP["muonplus"].P() << " " << newP["muonplus"].Pt() << " " << xprim <<  " " << yprim << " " << zprim << endl;

		// apply the cuts to the generated events, which now have momentum in lab frame (sec 3.1, point 3)
		vector<int> inacc = fiducialCuts( *newP["muonplus"], *newP["muonminus"], *newP["kaon"], *newP["Jpsi"] );

		// Now form the distribution for p, pt of each of the particles
		// Weight must be given by the product of the efficiencies of each particle
		// In the first iteration, the efficiencies are taken to be flat
		double w_angular = 1. - cosThetaMu*cosThetaMu; // why is this angular term needed?
		double w_mup = effMuonP(*newP["muonplus"]) * inacc[0];
		double w_mum = effMuonM(*newP["muonminus"]) * inacc[1];
		double w_k = effKaon(*newP["kaon"]) * inacc[2];

		double w_ev = w_angular * inacc[3] * w_mup * w_mum * w_k;
		w_sum_gen += w_ev;

		// I really don't know why this is necessary
		map<string, double> w_no_eff;
		w_no_eff["muonplus"]  = w_angular * inacc[3] * w_mum * w_k   * inacc[0];
		w_no_eff["muonminus"] = w_angular * inacc[3] * w_mup * w_k   * inacc[1];
		w_no_eff["kaon"]      = w_angular * inacc[3] * w_mup * w_mum * inacc[2];

		/* 
		   cout << w_angular << " "<< effMuonP(*newP["muonplus"]) << " " << effMuonM(*newP["muonminus"]) << " " << effKaon(*newP["kaon"]) << endl;
		   cout << w_angular << " "<< inacc[0] << " " << inacc[1] << " " << inacc[2] << " " << inacc[3] << endl;
		   cout << w_angular << " " << w_sum_gen << endl;
		   */

		vector<string>::iterator particleIterator;
		for ( particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator )
		{
			string particle = *particleIterator;
			cosThetaMuTmpHistos[particle]->Fill(cosThetaMu, w_ev);
			cosThetaKTmpHistos[particle]->Fill(cosThetaK, w_ev);
			pTmpHistos[particle]->Fill(newP[particle]->P(), w_ev);//w_no_eff[particle]);//w_ev);
			ptTmpHistos[particle]->Fill(newP[particle]->Pt(), w_no_eff[particle]);//w_ev);
			pptTmpHistos[particle]->Fill(newP[particle]->P(), newP[particle]->Pt(), w_no_eff[particle]);//w_ev);
			delete newP[particle];
		}
		delete newP["Jpsi"];
	}
}

void PerEventAngularAcceptance::boost()
{
	// In what follows, we will need to boost to the B rest frame
	// Inside the generation loop, we will also have to boost to the
	// Jpsi rest frame to work out the muon angular distribution so
	// do some Jpsi stuff here.
	boosttoJpsi = -(p["Jpsi"]->BoostVector());
	boostBackToLabFromJpsi = p["Jpsi"]->BoostVector();

	// Now do the boost
	vector<string>::iterator particleIterator;
	for ( particleIterator = partsAndRes.begin(); particleIterator != partsAndRes.end(); ++particleIterator )
	{
		string particle = *particleIterator;
		// Make copy of the vector, otherwise the p vector will be modified by the boost
		pInJpsi[particle] = new TLorentzVector(*p[particle]);
		pInJpsi[particle]->Boost(boosttoJpsi);
	}
	// We need to boost to the B rest frame, generate some events
	// and then boost back to the lab frame.
	boosttoB = -(p["B"]->BoostVector());
	boostBackToLabFromB = p["B"]->BoostVector();

	// Now do the boost
	for ( particleIterator = partsAndRes.begin(); particleIterator != partsAndRes.end(); ++particleIterator )
	{
		string particle = *particleIterator;
		// Make copy of the vector, otherwise the p vector will be modified by the boost
		pInB[particle] = new TLorentzVector(*p[particle]);
		pInB[particle]->Boost(boosttoB);
	}
}


void PerEventAngularAcceptance::writeHistos()
{
	string outFileName = "out.root";
	TFile * outFile = new TFile( outFileName.c_str(), "RECREATE");
	outFile->SetCompressionLevel(9);
	vector<string>::iterator particleIterator;
	for ( particleIterator = particles.begin(); particleIterator != particles.end(); ++particleIterator )
	{
		string particle = *particleIterator;  
		cosThetaMuHistos[particle]->Write();
		cosThetaKHistos[particle]->Write();
		cosThetaMuHistosPredicted[particle]->Write();
		cosThetaKHistosPredicted[particle]->Write();
		pHistos[particle]->Write();
		ptHistos[particle]->Write();
		pptHistos[particle]->Write();
		pHistosPredicted[particle]->Write();
		ptHistosPredicted[particle]->Write();
		pptHistosPredicted[particle]->Write();
		pHistosEffRatio[particle]->Write();
		pHistosEffRatioPrev[particle]->Write();
		ptHistosEffRatio[particle]->Write();
		pptHistosEffRatio[particle]->Write();
	}
	outFile->Close();
	delete outFile;
}  

double PerEventAngularAcceptance::effMuonP(TLorentzVector P)
{
	double p = P.P();
	double pt = P.Pt();
	double py = P.Py();
	double pz = P.Pz();

	//cout << pz <<  " " << p << " " << pt << " " << pt/p << " " << fabs(py/p) << endl;

	if (pz<0) return 0;
	if (p<pmin) return 0;
	if (pt<ptmin) return 0;
	if (pt/p<amin || pt/p>amax || fabs(py/p)>aymax) return 0;

	int j1 = int ((p-pmin)/pbinWidth) + 1;
	if (j1 > pbins) j1 = pbins;
	int j2 = int ((pt-ptmin)/ptbinWidth) + 1;
	if (j2 > ptbins) j2 = ptbins;

	return pptHistosEffRatio["muonplus"]->GetBinContent(j1, j2);
}

double PerEventAngularAcceptance::effMuonM(TLorentzVector P)
{
	double p = P.P();
	double pt = P.Pt();
	double py = P.Py();
	double pz = P.Pz();

	if (pz < 0) return 0;
	if (p < pmin) return 0;
	if (pt < ptmin) return 0;
	if (pt/p < amin || pt/p > amax || fabs(py/p) > aymax) return 0;

	int j1 = int ((p-pmin)/pbinWidth) + 1;
	if(j1 > pbins) j1 = pbins;
	int j2 = int ((pt-ptmin)/ptbinWidth) + 1;
	if(j2 > ptbins) j2 = ptbins;

	return pptHistosEffRatio["muonminus"]->GetBinContent(j1, j2);
}

double PerEventAngularAcceptance::effKaon(TLorentzVector P)
{

	double p = P.P();
	double pt = P.Pt();
	double py = P.Py();
	double pz = P.Pz();

	if (pz < 0) return 0;
	if (p < pmin) return 0;
	if (pt < ptmin) return 0;
	if (pt/p < amin || pt/p > amax || fabs(py/p) > aymax) return 0;

	int j1 = int ((p-pmin)/pbinWidth) + 1;
	if(j1 > pbins) j1 = pbins;
	int j2 = int ((pt-ptmin)/ptbinWidth) + 1;
	if(j2 > ptbins) j2 = ptbins;

	return pptHistosEffRatio["kaon"]->GetBinContent(j1, j2);
}

double PerEventAngularAcceptance::cosk(TLorentzVector particle, TLorentzVector parent) {

	TVector3 boosttoparent = -(parent.BoostVector());

	particle.Boost(boosttoparent);

	TVector3 particle3 = particle.Vect();
	TVector3 zaxis(0,0,1);

	double numerator = particle3.Dot(zaxis);
	double denominator = (particle3.Mag())*(zaxis.Mag());
	double temp = 2;

	if(fabs(denominator-0)<DOUBLE_TOLERANCE) temp= numerator/denominator;

	return temp;
}

double PerEventAngularAcceptance::coshel(TLorentzVector particle, TLorentzVector parent,
		TLorentzVector grandparent) {

	TVector3 boosttoparent = -(parent.BoostVector());

	particle.Boost(boosttoparent);
	grandparent.Boost(boosttoparent);

	TVector3 particle3 = particle.Vect();
	TVector3 grandparent3 = -grandparent.Vect();

	Float_t numerator = Float_t(particle3.Dot(grandparent3));
	Float_t denominator = Float_t((particle3.Mag())*(grandparent3.Mag()));
	Float_t temp = 2;
	if(fabs(denominator-0)<DOUBLE_TOLERANCE) temp= numerator/denominator;

	return temp;
}

vector<int> PerEventAngularAcceptance::fiducialCuts( TLorentzVector newv4mup, TLorentzVector newv4mum, TLorentzVector newv4k, TLorentzVector newv4Jpsi )
{
	// Calculate if the particles are in the acceptance or not
	int inacc_mup = 1;
	int inacc_mum = 1;
	int inacc_k = 1;
	int inacc_jpsi = 1;

	// Now apply the cuts
	if (newv4mup.P() < pmin) inacc_mup = 0;
	if (newv4mup.Pt()< ptmin) inacc_mup = 0;
	if (newv4mup.Pz() < 0) inacc_mup = 0;
	if (newv4mup.Pt()/newv4mup.P() < amin || newv4mup.Pt()/newv4mup.P() > amax || fabs(newv4mup.Py()/newv4mup.P()) > aymax) inacc_mup = 0;

	if (newv4mum.P() < pmin) inacc_mum = 0;
	if (newv4mum.Pt() < ptmin) inacc_mum = 0;
	if (newv4mum.Pz() < 0) inacc_mum = 0;
	if (newv4mum.Pt()/newv4mum.P() < amin || newv4mum.Pt()/newv4mum.P() > amax || fabs(newv4mum.Py()/newv4mum.P()) > aymax) inacc_mum = 0;

	if (newv4k.P() < pmin) inacc_k = 0;
	if (newv4k.Pt() < ptmin) inacc_k = 0;
	if (newv4k.Pz() < 0) inacc_k = 0;
	if (newv4k.Pt()/newv4k.P() < amin || newv4k.Pt()/newv4k.P() > amax || fabs(newv4k.Py()/newv4k.P()) > aymax) inacc_k = 0;

	if (newv4Jpsi.P() < pmin_jpsi) inacc_jpsi = 0;
	if (newv4Jpsi.Pt() < ptmin_jpsi) inacc_jpsi = 0;

	vector<int> res;
	res.push_back(inacc_mup);
	res.push_back(inacc_mum);
	res.push_back(inacc_k);
	res.push_back(inacc_jpsi);

	return res;
}
