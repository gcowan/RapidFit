
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TNtuple.h"
#include "TString.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TGraph2D.h"

#include "EdStyle.h"
#include "TTree_Processing.h"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace::std;

int main( int argc, char* argv[] )
{
	if( argc < 4 ) return -1;

	TString file_full_analysis = argv[1];
	TString file_CPEven_analysis = argv[2];
	TString file_CPOdd_analysis = argv[3];

	TString TTree_Name = "RapidFitResult";

	TFile* All_Data = new TFile( file_full_analysis, "READ" );
	TTree* full_fit = (TTree*) gDirectory->Get( TTree_Name );

	TFile* Even_Data = new TFile( file_CPEven_analysis, "READ" );
	TTree* even_fit = (TTree*) gDirectory->Get( TTree_Name );

	TFile* Odd_Data = new TFile( file_CPOdd_analysis, "READ" );
	TTree* odd_fit = (TTree*) gDirectory->Get( TTree_Name );

	unsigned int total_entries = (unsigned)full_fit->GetEntries();

	vector<double>* full_dG;
	vector<double>* full_dG_error;
	vector<double>* full_fit_status;

	vector<double>* tau_l;
	vector<double>* tau_l_err;
	vector<double>* tau_l_fit_status;
	vector<double>* tau_h;
	vector<double>* tau_h_err;
	vector<double>* tau_h_fit_status;

	vector<double>* diff_val = new vector<double>();

	vector<double>* weighted_dG = new vector<double>();
	vector<double>* weighted_dG_error = new vector<double>();

	full_dG = TTree_Processing::Buffer_Branch( full_fit, "deltaGamma_value", "" );
	full_dG_error = TTree_Processing::Buffer_Branch( full_fit, "deltaGamma_error", "" );

	full_fit_status = TTree_Processing::Buffer_Branch( full_fit, "Fit_Status", "" );

	tau_l = TTree_Processing::Buffer_Branch( even_fit, "tau_value", "" );
	tau_l_err = TTree_Processing::Buffer_Branch( even_fit, "tau_error", "" );
	tau_l_fit_status = TTree_Processing::Buffer_Branch( even_fit, "Fit_Status", "" );

	tau_h = TTree_Processing::Buffer_Branch( odd_fit, "tau_value", "" );
	tau_h_err = TTree_Processing::Buffer_Branch( odd_fit, "tau_error", "" );
	tau_h_fit_status = TTree_Processing::Buffer_Branch( even_fit, "Fit_Status", "" );

	unsigned int max_events = full_dG->size();
	if( tau_l->size() < max_events ) max_events = tau_l->size();
	if( tau_h->size() < max_events ) max_events = tau_h->size();

	for( unsigned int i=0; i< max_events; ++i )
	{
		if( fabs( (*full_fit_status)[i] - 3. ) > 1E-3 ) continue;
		if( fabs( (*tau_l_fit_status)[i] - 3. ) > 1E-3 ) continue;
		if( fabs( (*tau_h_fit_status)[i] - 3. ) > 1E-3 ) continue;

		double this_dG = (1./((*tau_l)[i]) ) - (1./((*tau_h)[i]) );
		weighted_dG->push_back( this_dG );

		double this_err = this_dG * sqrt( ((*tau_l_err)[i]/(*tau_l)[i])*((*tau_l_err)[i]/(*tau_l)[i]) + ((*tau_h_err)[i]/(*tau_h)[i])*((*tau_h_err)[i]/(*tau_h)[i]) );
		weighted_dG_error->push_back( fabs(this_err) );
	}

	for( unsigned int i=0; i< weighted_dG->size(); ++i )
	{
		if( fabs( (*full_fit_status)[i] - 3. ) > 1E-3 ) continue;
		if( fabs( (*tau_l_fit_status)[i] - 3. ) > 1E-3 ) continue;
		if( fabs( (*tau_h_fit_status)[i] - 3. ) > 1E-3 ) continue;
		cout << "true: " << (*full_dG)[i] << " \\pm " << (*full_dG_error)[i] << "\t\t";
		cout << "weighted: " << (*weighted_dG)[i] << " \\pm  " << (*weighted_dG_error)[i] << endl;
		diff_val->push_back( (*full_dG)[i] - (*weighted_dG)[i] );
	}

	//exit(0);

	TH1D* diff_th = new TH1D( "diff_th", "diff_th", 60, -0.15, 0.15 );

	for( unsigned int i=0; i< diff_val->size(); ++i )
	{
		diff_th->Fill( (*diff_val)[i] );
	}

	TCanvas* c1 = EdStyle::RapidFitCanvas( "canv", "canv" );

	diff_th->Draw();

	c1->Update();

	c1->Print("test.pdf");

	return 0;
}


