/**
        @class PerEventAngularAcceptance 

        @author Greig Cowan greig.cowan@cern.ch
	      @date 2010-05-03
*/

#pragma once
#ifndef PEREVENTANGULARACCEPTANCE_H
#define PEREVENTANGULARACCEPTANCE_H

//	ROOT Headers
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TBranch.h"
#include "TEventList.h"
#include "TH1F.h"
#include "TH2F.h"
//	System Headers
#include <vector>
#include <map>

using namespace::std;

class PerEventAngularAcceptance
{
	public:
		PerEventAngularAcceptance( string fileName, string outFileName, string tupleName );
		PerEventAngularAcceptance();
		~PerEventAngularAcceptance();

		void setupNtuple( string, string);
		void fillEffHistos( int iteration );
		void loopOnReconstructedBs();
		void writeHistos();
		void boost();
		void generateEventsAndCalculateFgivenB();
		double effMuonP( TLorentzVector P );
		double effMuonM( TLorentzVector P );
		double effKaon( TLorentzVector P );
		vector<int> fiducialCuts( TLorentzVector newv4mup, TLorentzVector newv4mum, TLorentzVector newv4k, TLorentzVector newv4Jpsi );
		double cosk( TLorentzVector particle, TLorentzVector parent);
		double coshel( TLorentzVector particle, TLorentzVector parent, TLorentzVector grandparent);

		TFile * tupleFile;
		TNtuple * decaytree;
		TEventList * evtList;
		int nev;
		vector<string> particles;
		vector<string> resonances;
		vector<string> partsAndRes;

		// Will hold the momenta for each particle in different rest frames
		map<string, TLorentzVector *> p;       // Lab frame
		map<string, TLorentzVector *> pInB;    // B rest frame
		map<string, TLorentzVector *> pInJpsi; // Jpsi rest frame

		// Histograms that we will need to calculate the efficiencies
		map<string, TH1F *> cosThetaMuHistos;
		map<string, TH1F *> cosThetaMuHistosPredicted;
		map<string, TH1F *> cosThetaMuTmpHistos; // used to store quantities during event generation
		map<string, TH1F *> cosThetaMuSumHistos; // Weighted sum of all the 
		map<string, TH1F *> cosThetaMuHistosEffRatio; // The efficiency we are after 
		map<string, TH1F *> cosThetaKHistos;
		map<string, TH1F *> cosThetaKHistosPredicted;
		map<string, TH1F *> cosThetaKTmpHistos;
		map<string, TH1F *> cosThetaKSumHistos;
		map<string, TH1F *> cosThetaKHistosEffRatio;
		map<string, TH1F *> pHistos;
		map<string, TH1F *> pHistosPredicted;
		map<string, TH1F *> pTmpHistos; // used to store quantities during event generation
		map<string, TH1F *> pSumHistos; // Weighted sum of all the 
		map<string, TH1F *> pHistosEffRatio; // The efficiency we are after 
		map<string, TH1F *> pHistosEffRatioPrev; // The efficiency we are after 
		map<string, TH1F *> ptHistos;
		map<string, TH1F *> ptHistosPredicted;
		map<string, TH1F *> ptTmpHistos;
		map<string, TH1F *> ptSumHistos;
		map<string, TH1F *> ptHistosEffRatio;
		map<string, TH2F *> pptHistos;
		map<string, TH2F *> pptHistosPredicted;
		map<string, TH2F *> pptTmpHistos;
		map<string, TH2F *> pptSumHistos;
		map<string, TH2F *> pptHistosEffRatio;

		TVector3 boosttoJpsi, boostBackToLabFromJpsi;
		TVector3 boosttoB, boostBackToLabFromB;

		double GeV;
		double pmin, pmax, pbinWidth;
		double ptmin, ptmax, ptbinWidth;
		double amin, amax, aymax;
		double pmin_jpsi, ptmin_jpsi;

		int cosbins, pbins, ptbins;
		
		double w_sum_gen;

		// Required to read stuff from the ntuple
		float muonplus_PX, muonplus_PY, muonplus_PZ, muonplus_PE;
		float muonminus_PX, muonminus_PY, muonminus_PZ, muonminus_PE;
		float kaon_PX, kaon_PY, kaon_PZ, kaon_PE;
		float Bu_MM, Bu_TAU, Bu_BKGCAT;

		TBranch   *b_Bu_BKGCAT;
		TBranch   *b_Bu_MM;
		TBranch   *b_Bu_TAU;

		TBranch   *b_muonplus_PX;
		TBranch   *b_muonplus_PY;
		TBranch   *b_muonplus_PZ;
		TBranch   *b_muonplus_PE;
		TBranch   *b_muonminus_PX;
		TBranch   *b_muonminus_PY;
		TBranch   *b_muonminus_PZ;
		TBranch   *b_muonminus_PE;
		TBranch   *b_kaon_PX;
		TBranch   *b_kaon_PY;
		TBranch   *b_kaon_PZ;
		TBranch   *b_kaon_PE;

	protected:
	  
	private:
		//	Uncopyable!
		PerEventAngularAcceptance ( const PerEventAngularAcceptance& );
		PerEventAngularAcceptance& operator = ( const PerEventAngularAcceptance& );
};
#endif
