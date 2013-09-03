
#include "EdStyle.h"
#include "TTree_Processing.h"
#include "ROOT_File_Processing.h"
#include "Histo_Processing.h"
#include "Template_Functions.h"

#include "TROOT.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TPaveStats.h"

#include <iomanip>
#include <vector>
#include <cstdlib>

using namespace::std;

int main( int argc, char* argv[] )
{
	if( argc < 3 || argc > 4 )
	{
		cout << "Usage: " << argv[0] << " gammaL.root gammaH.root" << endl << endl;
		exit(0);
	}

	//      Setup the Canvas and such
	EdStyle* RapidFit_Style = new EdStyle();
	RapidFit_Style->SetStyle();

	RapidFit_Style->SetColzPlotStyle();

	gStyle->SetOptFit(111);
	gStyle->SetPalette(1);

	//	Open First TTree in the file
	TTree* gammaLTree = ROOT_File_Processing::GetFirstTree( argv[1] );

	vector<double>* gammaL_tau_value = TTree_Processing::Buffer_Branch( gammaLTree, "1./tau_value" );
	vector<double>* gammaL_tau_error = TTree_Processing::Buffer_Branch( gammaLTree, "(1./tau_value)*(tau_error/tau_value)" );
	vector<double>* gammaL_FitStatus = TTree_Processing::Buffer_Branch( gammaLTree, "Fit_Status" );

	//	Open First TTree in the file
	TTree* gammaHTree = ROOT_File_Processing::GetFirstTree( argv[2] );

	vector<double>* gammaH_tau_value = TTree_Processing::Buffer_Branch( gammaHTree, "1./tau_value" );
	vector<double>* gammaH_tau_error = TTree_Processing::Buffer_Branch( gammaHTree, "(1./tau_value)*(tau_error/tau_value)" );
	vector<double>* gammaH_FitStatus = TTree_Processing::Buffer_Branch( gammaHTree, "Fit_Status" );

	cout << "Storing GLGH Data" << endl;
	TFile* newFile = new TFile( "GLGHData.root", "RECREATE" );
	TTree* GLGHTree = new TTree( "GLGHTree", "GLGHTree" );

	TTree_Processing::AddBranch( GLGHTree, "gammaH_value", *gammaH_tau_value );
	TTree_Processing::AddBranch( GLGHTree, "gammaH_error", *gammaH_tau_error );
	TTree_Processing::AddBranch( GLGHTree, "gammaL_value", *gammaL_tau_value );
	TTree_Processing::AddBranch( GLGHTree, "gammaL_error", *gammaL_tau_error );

	GLGHTree->Write("",TObject::kOverwrite);
	newFile->Write("",TObject::kOverwrite);

	cout << "Stored!" << endl;

	TTree* gammaTree = NULL;

	vector<double>* gamma_value = NULL;
	vector<double>* deltaGamma_value = NULL;
	vector<double>* FitStatus_value = NULL;

	if( argc == 4 )
	{
		gammaTree = ROOT_File_Processing::GetFirstTree( argv[3] );
		gamma_value = TTree_Processing::Buffer_Branch( gammaTree, "gamma_value" );
		deltaGamma_value = TTree_Processing::Buffer_Branch( gammaTree, "deltaGamma_value" ); 
		FitStatus_value = TTree_Processing::Buffer_Branch( gammaTree, "Fit_Status" );
	}

	vector<double> gammaS_value, gammaS_error, gammaS_pull, deltaGammaS_value, deltaGammaS_error, deltaGammaS_pull;

	vector<double> gamma_diff, deltaGamma_diff;

	vector<double> angGammaTrudiff, untagGammaTrudiff, angdGammaTrudiff, untagdGammaTrudiff;

	TGraphErrors* corr_graph = new TGraphErrors( gammaH_tau_value->size(), &((*gammaL_tau_value)[0]), &((*gammaH_tau_value)[0]),
			&((*gammaL_tau_error)[0]), &((*gammaH_tau_error)[0]) );

	TGraphErrors* corr_graph2 = new TGraphErrors( gammaH_tau_value->size(), &((*gammaH_tau_value)[0]), &((*gammaL_tau_value)[0]),
			&((*gammaH_tau_error)[0]), &((*gammaL_tau_error)[0]) );

	double GL_sum=0., GH_sum=0.;
	double GL_sum_err=0., GH_sum_err=0.;
	double GL_sum_sq=0., GH_sum_sq=0.;
	double GL_sum_err_sq=0., GH_sum_err_sq=0.;
	double GLGH_prod_sum=0.;
	double GLGH_prod_err_sum=0.;
	double GL_avr=0., GH_avr=0.;
	double GL_avr_err=0., GH_avr_err=0.;
	double N = (double)(gammaL_tau_value->size());

	for( unsigned int i=0; i< gammaL_tau_value->size(); ++i )
	{
		GL_sum += (*gammaL_tau_value)[i];
		GL_sum_err += (*gammaL_tau_error)[i]*(*gammaL_tau_error)[i]/((*gammaL_tau_value)[i]*(*gammaL_tau_value)[i]);
		GH_sum += (*gammaH_tau_value)[i];
		GH_sum_err += (*gammaH_tau_error)[i]*(*gammaH_tau_error)[i]/((*gammaH_tau_value)[i]*(*gammaH_tau_value)[i]);
		GL_sum_sq += (*gammaL_tau_value)[i]*(*gammaL_tau_value)[i];
		GL_sum_err_sq += (*gammaL_tau_value)[i]*(*gammaL_tau_value)[i]*sqrt( ((*gammaL_tau_error)[i]*(*gammaL_tau_error)[i])/((*gammaL_tau_value)[i]*(*gammaL_tau_value)[i]) *2. );
		GH_sum_sq += (*gammaH_tau_value)[i]*(*gammaH_tau_value)[i];
		GH_sum_err_sq += (*gammaH_tau_value)[i]*(*gammaH_tau_value)[i]*sqrt( ((*gammaH_tau_error)[i]*(*gammaH_tau_error)[i])/((*gammaH_tau_value)[i]*(*gammaH_tau_value)[i]) *2. );
		GLGH_prod_sum += (*gammaL_tau_value)[i]*(*gammaH_tau_value)[i];
		GLGH_prod_err_sum += (*gammaL_tau_value)[i]*(*gammaH_tau_value)[i]*sqrt( ((*gammaL_tau_error)[i]*(*gammaL_tau_error)[i])/((*gammaL_tau_value)[i]*(*gammaL_tau_value)[i])
					 + ((*gammaH_tau_error)[i]*(*gammaH_tau_error)[i])/((*gammaH_tau_value)[i]*(*gammaH_tau_value)[i]) );
	}

	
	cout << "a:" << endl;
	GL_sum_err = sqrt( GL_sum_err ); cout << GL_sum << " / " << GL_sum_err << endl;
	GH_sum_err = sqrt( GH_sum_err ); cout << GH_sum << " / " << GH_sum_err << endl;
	GL_sum_err_sq = sqrt( GL_sum_err_sq ); cout << GL_sum_sq << " / " << GL_sum_err_sq << endl;
	GH_sum_err_sq = sqrt( GH_sum_err_sq ); cout << GH_sum_sq << " / " << GH_sum_err_sq << endl;
	GLGH_prod_err_sum = sqrt( GLGH_prod_err_sum ); cout << GLGH_prod_sum << " / " << GLGH_prod_err_sum << endl;

	cout << "b:" << endl;

	GL_avr = GL_sum / (double)(gammaL_tau_value->size());
	GH_avr = GH_sum / (double)(gammaH_tau_value->size());

	GL_avr_err = GL_sum_err / (double)(gammaL_tau_value->size());
	GH_avr_err = GH_sum_err / (double)(gammaH_tau_value->size());

	cout << GL_avr << " / " << GL_avr_err << endl;
	cout << GH_avr << " / " << GH_avr_err << endl;

	double corr = N*GLGH_prod_sum - GL_sum*GH_sum;
	corr /= sqrt( (N*GL_sum_sq - GL_sum*GL_sum) * (N*GH_sum_sq - GH_sum*GH_sum) );

	double corr_err_numer = sqrt( GLGH_prod_err_sum*GLGH_prod_err_sum + GL_sum*GH_sum*sqrt((GL_sum_err*GL_sum_err)/(GL_sum_sq*GL_sum_sq) + (GH_sum_err*GH_sum_err)/(GH_sum_sq*GH_sum_sq) ));

	double corr_err_denom = sqrt( ( N*GL_sum_err_sq*GL_sum_err_sq + ( (GL_sum_err*GL_sum_err) *2. ) ) 
				    + ( N*GH_sum_err_sq*GH_sum_err_sq + ( (GH_sum_err*GH_sum_err) *2. ) ) );

	cout << corr_err_numer << " // " << corr_err_denom << endl;

	cout << N*GLGH_prod_sum - GL_sum*GH_sum << " // " << sqrt( (N*GL_sum_sq - GL_sum*GL_sum) * (N*GH_sum_sq - GH_sum*GH_sum) ) << endl;

	double corr_err = fabs(corr)* sqrt( corr_err_numer*corr_err_numer/((N*GLGH_prod_sum - GL_sum*GH_sum)*(N*GLGH_prod_sum - GL_sum*GH_sum)) + corr_err_denom*corr_err_denom/((N*GL_sum_sq - GL_sum*GL_sum) * (N*GH_sum_sq - GH_sum*GH_sum)) );
	
	GL_avr = GL_sum / (double)(gammaL_tau_value->size());
	GH_avr = GH_sum / (double)(gammaH_tau_value->size());

//	double corr = N*GLGH_prod_sum - GL_sum*GH_sum;
//	corr /= sqrt( (N*GL_sum_sq - GL_sum*GL_sum) * (N*GH_sum_sq - GH_sum*GH_sum) );


	cout << endl << setprecision(10) << "GL / GH Correlation: " << corr << " \\pm " << corr_err << endl << endl;

	double rho=-0.788;//85;//-0.7850;

	for( unsigned int i=0; i< gammaH_tau_error->size(); ++i )
	{
		(*gammaH_tau_error)[i] *=0.9;
	}

	for( unsigned int i=0; i< gammaL_tau_value->size(); ++i )
	{
		double gammaS_val=0., gammaS_err=0., delGammaS_val=0., delGammaS_err=0.;

		gammaS_val = 0.5 * ( (*gammaL_tau_value)[i] + (*gammaH_tau_value)[i] );
		gammaS_err = //( (*gammaL_tau_value)[i] + (*gammaH_tau_value)[i] )

			sqrt( (*gammaH_tau_error)[i]*(*gammaH_tau_error)[i]*0.5*0.5 + 0.5*0.5*(*gammaL_tau_error)[i]*(*gammaL_tau_error)[i] + 0.5*0.5*2.*rho*(*gammaH_tau_error)[i]*(*gammaL_tau_error)[i] );

			//			gammaS_val * sqrt( ((*gammaL_tau_error)[i]*(*gammaL_tau_error)[i])/((*gammaL_tau_value)[i]*(*gammaL_tau_value)[i]) + ((*gammaH_tau_error)[i]*(*gammaH_tau_error)[i])/((*gammaH_tau_value)[i]*(*gammaH_tau_value)[i]) - 0.8*(*gammaH_tau_error)[i]*(*gammaL_tau_error)[i]/((*gammaL_tau_value)[i]*(*gammaH_tau_value)[i])*2. );
/*				gammaS_val* ( sqrt( ((*gammaL_tau_error)[i]*(*gammaL_tau_error)[i])/((*gammaL_tau_value)[i]*(*gammaL_tau_value)[i])
				 		 + (*gammaH_tau_error)[i]*(*gammaH_tau_error)[i] )/((*gammaH_tau_value)[i]*(*gammaH_tau_value)[i])
				 		 + rho*2.*(*gammaH_tau_error)[i]*(*gammaL_tau_error)[i]/((*gammaL_tau_value)[i]*(*gammaH_tau_value)[i] ));
*/
		delGammaS_val = (*gammaL_tau_value)[i] - (*gammaH_tau_value)[i];
		delGammaS_err = sqrt( (*gammaL_tau_error)[i]*(*gammaL_tau_error)[i] + (*gammaH_tau_error)[i]*(*gammaH_tau_error)[i] + 2.*(fabs(rho))*(*gammaL_tau_error)[i]*(*gammaH_tau_error)[i] );
//		delGammaS_err = sqrt( (*gammaL_tau_error)[i]*(*gammaL_tau_error)[i]/((*gammaL_tau_value)[i]*(*gammaL_tau_value)[i])+(*gammaH_tau_error)[i]*(*gammaH_tau_error)[i]/((*gammaH_tau_value)[i]*(*gammaH_tau_value)[i]) - rho*2.*(*gammaH_tau_error)[i]*(*gammaL_tau_error)[i]/((*gammaH_tau_value)[i]*(*gammaL_tau_value)[i]) );

		gammaS_value.push_back( gammaS_val );
		gammaS_error.push_back( gammaS_err );
		gammaS_pull.push_back( (gammaS_val-0.6715)/gammaS_err );

		deltaGammaS_value.push_back( delGammaS_val );
		deltaGammaS_error.push_back( delGammaS_err );
		deltaGammaS_pull.push_back( (delGammaS_val-0.1)/delGammaS_err );
	}

	int thisOptStat = gStyle->GetOptStat();

	gStyle->SetOptStat(0);

	TH1* gamma_histo = new TH1D( "gamma_dist", "gamma_dist", gammaS_value.size(), get_minimum(gammaS_value), get_maximum(gammaS_value) );
	TH1* gammaErr_histo = new TH1D( "gammaErr_dist", "gammaErr_dist", gammaS_error.size(), get_minimum(gammaS_error), get_maximum(gammaS_error) );
	TH1* gammaPull_histo = new TH1D( "gammaPull_dist", "gammaPull_dist", gammaS_pull.size(), get_minimum(gammaS_pull), get_maximum(gammaS_pull) );
	for( unsigned int i=0; i< gammaS_value.size(); ++i ) gamma_histo->Fill( gammaS_value[i] );
	for( unsigned int i=0; i< gammaS_error.size(); ++i ) gammaErr_histo->Fill( gammaS_error[i] );
	for( unsigned int i=0; i< gammaS_pull.size(); ++i ) gammaPull_histo->Fill( gammaS_pull[i] );

	TH1* dGamma_histo = new TH1D( "deltaGamma_dist", "deltaGamma_dist", deltaGammaS_value.size(), get_minimum(deltaGammaS_value), get_maximum(deltaGammaS_value) );
	TH1* dGammaErr_histo = new TH1D( "deltaGammaErr_dist", "deltaGammaErr_dist", deltaGammaS_error.size(), get_minimum(deltaGammaS_error), get_maximum(deltaGammaS_error) );
	TH1* dGammaPull_histo = new TH1D( "deltaGammaPull_dist", "deltaGammaPull_dist", deltaGammaS_pull.size(), get_minimum(deltaGammaS_pull), get_maximum(deltaGammaS_pull) );
	for( unsigned int i=0; i< gammaS_value.size(); ++i ) dGamma_histo->Fill( deltaGammaS_value[i] );
	for( unsigned int i=0; i< gammaS_error.size(); ++i ) dGammaErr_histo->Fill( deltaGammaS_error[i] );
	for( unsigned int i=0; i< gammaS_pull.size(); ++i ) dGammaPull_histo->Fill( deltaGammaS_pull[i] );

	cout << "gamma_output" << endl;

	TCanvas* c1 = EdStyle::RapidFitCanvas( "gamma_C", "gamma_C" );

	gamma_histo->Draw();
	c1->Update();
	Histogram_Processing::OptimallyRebin( gamma_histo );
	TString fit_type = Histogram_Processing::Best_Fit_Function( gamma_histo );
	Histogram_Processing::Silent_Fit( gamma_histo, fit_type );
	gamma_histo->Draw();
	c1->Update();
	TPaveStats* thisStats = (TPaveStats*)gamma_histo->GetListOfFunctions()->FindObject("stats");
	if( thisStats != NULL )
	{
		thisStats->SetFillStyle(3023);
		thisStats->SetTextColor(1);
		thisStats->Draw();
	}

	gamma_histo->GetXaxis()->SetTitle("#Gamma value [ps^{-1}]");
	gamma_histo->GetYaxis()->SetTitle("Candidates");
	c1->Update();
	c1->Print("gammaOutput.pdf");

	c1->cd();
	gammaErr_histo->Draw();
	c1->Update();
	Histogram_Processing::OptimallyRebin( gammaErr_histo );
	fit_type = Histogram_Processing::Best_Fit_Function( gammaErr_histo );
	Histogram_Processing::Silent_Fit( gammaErr_histo, fit_type );
	gammaErr_histo->Draw();
	c1->Update();
	thisStats = (TPaveStats*)gammaErr_histo->GetListOfFunctions()->FindObject("stats");
	if( thisStats != NULL )
	{
		thisStats->SetFillStyle(3023);
		thisStats->SetTextColor(1);
		thisStats->Draw();
	}

	gammaErr_histo->GetXaxis()->SetTitle("#Gamma error [ps^{-1}]");
	gammaErr_histo->GetYaxis()->SetTitle("Candidates");
	c1->Update();
	c1->Print("gammaErrOutput.pdf");

	c1->cd();
	gammaPull_histo->Draw();
	c1->Update();
	Histogram_Processing::OptimallyRebin( gammaPull_histo );
	fit_type = Histogram_Processing::Best_Fit_Function( gammaPull_histo );
	Histogram_Processing::Silent_Fit( gammaPull_histo, fit_type );
	gammaPull_histo->Draw();
	c1->Update();
	thisStats = (TPaveStats*)gammaPull_histo->GetListOfFunctions()->FindObject("stats");
	if( thisStats != NULL )
	{
		thisStats->SetFillStyle(3023);
		thisStats->SetTextColor(1);
		thisStats->Draw();
	}

	gammaPull_histo->GetXaxis()->SetTitle("#Gamma pull [ps^{-1}]");
	gammaPull_histo->GetYaxis()->SetTitle("Candidates");
	c1->Update();
	c1->Print("gammaPullOutput.pdf");

	cout << "dGamma_output" << endl;

	TCanvas* c2 = EdStyle::RapidFitCanvas( "dgamma_C", "dgamma_C" );

	dGamma_histo->Draw();
	c2->Update();
	Histogram_Processing::OptimallyRebin( dGamma_histo );
	fit_type = Histogram_Processing::Best_Fit_Function( dGamma_histo );
	Histogram_Processing::Silent_Fit( dGamma_histo, fit_type );
	dGamma_histo->Draw();
	c2->Update();
	thisStats = (TPaveStats*)dGamma_histo->GetListOfFunctions()->FindObject("stats");
	if( thisStats != NULL )
	{
		thisStats->SetFillStyle(3023);
		thisStats->SetTextColor(1);
		thisStats->Draw();
	}
	dGamma_histo->GetXaxis()->SetTitle("#Delta#Gamma value [ps^{-1}]");
	dGamma_histo->GetYaxis()->SetTitle("Candidates");
	c2->Update();
	c2->Print("deltaGammaOutput.pdf");

	c2->cd();
	dGammaErr_histo->Draw();
	c2->Update();
	Histogram_Processing::OptimallyRebin( dGammaErr_histo );
	fit_type = Histogram_Processing::Best_Fit_Function( dGammaErr_histo );
	Histogram_Processing::Silent_Fit( dGammaErr_histo, fit_type );
	dGammaErr_histo->Draw();
	c2->Update();
	thisStats = (TPaveStats*)dGammaErr_histo->GetListOfFunctions()->FindObject("stats");
	if( thisStats != NULL )
	{
		thisStats->SetFillStyle(3023);
		thisStats->SetTextColor(1);
		thisStats->Draw();
	}
	dGammaErr_histo->GetXaxis()->SetTitle("#Delta#Gamma error [ps^{-1}]");
	dGammaErr_histo->GetYaxis()->SetTitle("Candidates");
	c2->Update();
	c2->Print("deltaGammaErrOutput.pdf");

	c2->cd();
	dGammaPull_histo->Draw();
	c2->Update();
	Histogram_Processing::OptimallyRebin( dGammaPull_histo );
	fit_type = Histogram_Processing::Best_Fit_Function( dGammaPull_histo );
	Histogram_Processing::Silent_Fit( dGammaPull_histo, fit_type );
	dGammaPull_histo->Draw();
	c2->Update();
	thisStats = (TPaveStats*)dGammaPull_histo->GetListOfFunctions()->FindObject("stats");
	if( thisStats != NULL )
	{
		thisStats->SetFillStyle(3023);
		thisStats->SetTextColor(1);
		thisStats->Draw();
	}
	dGammaPull_histo->GetXaxis()->SetTitle("#Delta#Gamma pull [ps^{-1}]");
	dGammaPull_histo->GetYaxis()->SetTitle("Candidates");
	c2->Update();
	c2->Print("deltaGammaPullOutput.pdf");

	gStyle->SetOptStat( 0 );//thisOptStat );

	if( argc == 4 )
	{
		for( unsigned int i=0; i< gammaL_tau_value->size(); ++i )
		{
			if( (*gammaH_FitStatus)[i] == 3 )
			{
				if( (*gammaL_FitStatus)[i] == 3 )
				{
					if( (*FitStatus_value)[i] == 3 )
					{
						gamma_diff.push_back( gammaS_value[i] - (*gamma_value)[i] );

						angGammaTrudiff.push_back( gammaS_value[i] );
						untagGammaTrudiff.push_back( (*gamma_value)[i] );

						deltaGamma_diff.push_back( deltaGammaS_value[i] - (*deltaGamma_value)[i] );

						angdGammaTrudiff.push_back( deltaGammaS_value[i] );
						untagdGammaTrudiff.push_back( (*deltaGamma_value)[i] );
					}
				}
			}
		}

		int rebin_factor=1;

		TH1* delta_gamma_histo = new TH1D( "gammadiff_dist", "gammadiff_dist", gamma_diff.size(), get_minimum(gamma_diff), get_maximum(gamma_diff) );
		for( unsigned int i=0; i< gamma_diff.size(); ++i ) delta_gamma_histo->Fill( gamma_diff[i] );

		TH1* delta_dGamma_histo = new TH1D( "dGammadiff_dist", "dGammadiff_dist", deltaGamma_diff.size(), get_minimum(deltaGamma_diff), get_maximum(deltaGamma_diff) );
		for( unsigned int i=0; i< deltaGamma_diff.size(); ++i ) delta_dGamma_histo->Fill( deltaGamma_diff[i] );

		int gamma_bins=33;
		int dGamma_bins=33;

		//TH1* gamma_truth =new TH1D( "gamma_tru", "gamma_tru", angGammaTrudiff.size(), get_minimum(angGammaTrudiff), get_maximum(angGammaTrudiff) );
		TH1* gamma_truth =new TH1D( "gamma_tru", "gamma_tru", gamma_bins, get_minimum(angGammaTrudiff), get_maximum(angGammaTrudiff) );
		for( unsigned int i=0; i< angGammaTrudiff.size(); ++i ) gamma_truth->Fill( angGammaTrudiff[i] );
		//Histogram_Processing::OptimallyRebin( gamma_truth );
		gamma_truth->SetLineColor( 2 );
		//TH1* gamma_untag =new TH1D( "gamma_untag", "gamma_untag", untagGammaTrudiff.size(), get_minimum(angGammaTrudiff), get_maximum(angGammaTrudiff) );
		TH1* gamma_untag =new TH1D( "gamma_untag", "gamma_untag", gamma_bins, get_minimum(angGammaTrudiff), get_maximum(angGammaTrudiff) );
		for( unsigned int i=0; i< untagGammaTrudiff.size(); ++i ) gamma_untag->Fill( untagGammaTrudiff[i] );
		//Histogram_Processing::OptimallyRebin( gamma_untag );
		gamma_untag->SetLineColor( 3 );

		rebin_factor = int( ( (double) gamma_untag->GetNbinsX() ) / ( (double) gamma_truth->GetNbinsX() ) );
		//gamma_untag->Rebin( rebin_factor ); gamma_untag->SetLineColor( 3 );

		//TH1* dGamma_truth = new TH1D( "deltaGamma_tru", "deltaGamma_tru", angdGammaTrudiff.size(), get_minimum(angdGammaTrudiff), get_maximum(angdGammaTrudiff) );
		TH1* dGamma_truth = new TH1D( "deltaGamma_tru", "deltaGamma_tru", dGamma_bins, get_minimum(angdGammaTrudiff), get_maximum(angdGammaTrudiff) );
		for( unsigned int i=0; i< angdGammaTrudiff.size(); ++i ) dGamma_truth->Fill( angdGammaTrudiff[i] );
		//Histogram_Processing::OptimallyRebin( dGamma_truth );
		dGamma_truth->SetLineColor( 2 );
		//TH1* dGamma_untag =new TH1D( "deltaGamma_untag", "deltaGamma_untag", untagdGammaTrudiff.size(), get_minimum(angdGammaTrudiff), get_maximum(angdGammaTrudiff) );
		TH1* dGamma_untag =new TH1D( "deltaGamma_untag", "deltaGamma_untag", dGamma_bins, get_minimum(angdGammaTrudiff), get_maximum(angdGammaTrudiff) );
		for( unsigned int i=0; i< untagdGammaTrudiff.size(); ++i ) dGamma_untag->Fill( untagdGammaTrudiff[i] );
		//Histogram_Processing::OptimallyRebin( dGamma_untag );
		dGamma_untag->SetLineColor( 3 );

		rebin_factor = int( ( (double) dGamma_untag->GetNbinsX() ) / ( (double) dGamma_truth->GetNbinsX() ) );
		//dGamma_untag->Rebin( rebin_factor ); dGamma_untag->SetLineColor( 3 );

		cout << "diff_Gamma" << endl;

		TCanvas* c3 = EdStyle::RapidFitCanvas( "diffgamma_C", "diffgamma_C" );

		delta_gamma_histo->Draw();
		c3->Update();
		Histogram_Processing::OptimallyRebin( delta_gamma_histo );
		fit_type = Histogram_Processing::Best_Fit_Function( delta_gamma_histo );
		Histogram_Processing::Silent_Fit( delta_gamma_histo, fit_type );
		c3->Update();
		thisStats = (TPaveStats*)delta_gamma_histo->GetListOfFunctions()->FindObject("stats");
		if( thisStats != NULL )
		{
			thisStats->SetFillStyle(3023);
			thisStats->SetTextColor(1);
			thisStats->Draw();
		}
		delta_gamma_histo->GetXaxis()->SetTitle("#delta(#Gamma_{s}) [ps^{-1}]");
		delta_gamma_histo->GetYaxis()->SetTitle("Candidates");
		c3->Update();
		c3->Print("diff_Gamma.pdf");

		TCanvas* c4 = EdStyle::RapidFitCanvas( "diffdGamma_C", "diffdGamma_C" );

		delta_dGamma_histo->Draw();
		c4->Update();
		Histogram_Processing::OptimallyRebin( delta_dGamma_histo );
		fit_type = Histogram_Processing::Best_Fit_Function( delta_dGamma_histo );
		Histogram_Processing::Silent_Fit( delta_dGamma_histo, fit_type );
		c4->Update();
		thisStats = (TPaveStats*)delta_dGamma_histo->GetListOfFunctions()->FindObject("stats");
		if( thisStats != NULL )
		{
			thisStats->SetFillStyle(3023);
			thisStats->SetTextColor(1);
			thisStats->Draw();
		}
		delta_dGamma_histo->GetXaxis()->SetTitle("#delta(#Delta#Gamma_{s}) [ps^{-1}]");
		delta_dGamma_histo->GetYaxis()->SetTitle("Candidates");
		c4->Update();
		c4->Print("diff_deltaGamma.pdf");

		cout << "gamma_overlay" << endl;

		TCanvas* c5 = EdStyle::RapidFitCanvas( "gamma_overlay", "gamma_overlay" );

		gamma_untag->Draw();
		c5->Update();
		gamma_untag->GetXaxis()->SetTitle("#Gamma_{s} [ps^{-1}]");
		gamma_untag->GetYaxis()->SetTitle("Candidates");
		gamma_truth->Draw("SAME");
		c5->Update();
		TLegend* thisLegend = EdStyle::LHCbLegend();
		thisLegend->AddEntry( gamma_untag, "#Gamma_{s} Full Analysis", "l" );
		thisLegend->AddEntry( gamma_truth, "#Gamma_{s} Angular Moment Analysis", "l" );
		thisLegend->SetFillStyle(3023);
		thisLegend->SetTextColor(1);
		thisLegend->SetTextSize( thisLegend->GetTextSize()*0.5 );
		thisLegend->Draw();
		c5->Update();
		c5->Print("gamma_overlay.pdf");

		TCanvas* c6 = EdStyle::RapidFitCanvas( "dGamma_overlay", "dGamma_overlay" );

		dGamma_untag->Draw();
		c6->Update();
		dGamma_untag->GetXaxis()->SetTitle("#Delta#Gamma_{s} [ps^{-1}]");
		dGamma_untag->GetYaxis()->SetTitle("Candidates");
		dGamma_truth->Draw("SAME");
		c6->Update();
		TLegend* thisLegend2 = EdStyle::LHCbLegend();
		thisLegend2->AddEntry( dGamma_untag, "#Delta#Gamma_{s} Full Analysis", "l" );
		thisLegend2->AddEntry( dGamma_truth, "#Delta#Gamma_{s} Angular Moment Analysis", "l" );
		thisLegend2->SetFillStyle(3023);
		thisLegend2->SetTextColor(1);
		thisLegend2->SetTextSize( thisLegend2->GetTextSize()*0.5 );
		thisLegend2->Draw();
		c6->Update();
		c6->Print("deltaGamma_overlay.pdf");
	}
	//else
	//{
	gStyle->SetOptStat(0);

	cout << "Working on Correlation Graph" << endl;

	TCanvas* c7 = EdStyle::RapidFitCanvas( "corr_Graph", "corr_Graph" );

	corr_graph->Draw("AP");
	c7->Update();
	corr_graph->GetXaxis()->SetTitle("#Gamma_{L} [ps^{-1}]");
	corr_graph->GetYaxis()->SetTitle("#Gamma_{H} [ps^{-1}]");
	c7->Update();

	TString funcStr_Orig("("); funcStr_Orig+=1./corr;funcStr_Orig.Append("*x+[0])");

	TString funcStr("([1]*x+[0])");

	double xmin = corr_graph->GetXaxis()->GetXmin();
	double xmax = corr_graph->GetXaxis()->GetXmax();
	double ymin = corr_graph->GetYaxis()->GetXmin();
	double ymax = corr_graph->GetYaxis()->GetXmax();

	TF1* corr_funct = new TF1( "corr_func", funcStr, xmin, xmax );

	corr_funct->SetParameters( 1.18579, corr );

	corr_graph->SetLineWidth( 3 );

	corr_graph->Fit( corr_funct, "FEM" );

	corr_graph->SetLineWidth( 3 );

	double* x_func = new double[100];
	double* y_func = new double[100];

	double xstep = (xmax-xmin)/99.;
	for( unsigned int i=0; i< 100; ++i )
	{
		x_func[i] = xmin + xstep * i;
		y_func[i] = corr_funct->Eval( x_func[i] );
	}

	TGraph* thisCorrGraph = new TGraph( 100, x_func, y_func );
	thisCorrGraph->SetLineColor(3);

	thisCorrGraph->Draw("SAMEL");
	c7->Update();

	TF1* corr_funct2 = new TF1( "corr_func2", funcStr_Orig, xmin, xmax );

	corr_graph->Fit( corr_funct2, "FEM" );

	double* x_func2 = new double[100];
	double* y_func2 = new double[100];

	for( unsigned int i=0; i< 100; ++i )
	{
		x_func2[i] = xmin + xstep * i;
		y_func2[i] = corr_funct2->Eval( x_func2[i] );
	}

	TGraph* thisCorrGraph2 = new TGraph( 100, x_func2, y_func2 );

	thisCorrGraph2->SetLineColor( 2 );
	thisCorrGraph2->SetLineWidth( 3 );
	thisCorrGraph2->Draw("SAMEL");

	c7->Update();

	c7->Print("corr2DErr.pdf");

	cout << "build" << endl;
	TH2* thisHisto = new TH2D( "histo_colz2", "histo_colz2", 100, xmin, xmax, 100, ymin, ymax );
	cout << "fill" << endl;
	for( unsigned int i=0; i< gammaL_tau_value->size(); ++i )
	{
		thisHisto->Fill( gammaL_tau_value->at(i), gammaH_tau_value->at(i) );
	}

	int xbins =100;// /*100;*/Histogram_Processing::OptimumBinNumber( thisHisto, 1 );
	int ybins =100;// /*100;*/Histogram_Processing::OptimumBinNumber( thisHisto, 2 );

	cout << xbins << "\t\t" << ybins << endl;

	thisHisto = new TH2D( "histo_colz", "histo_colz", xbins, xmin, xmax, ybins, ymin, ymax );

	for( unsigned int i=0; i< gammaL_tau_value->size(); ++i )
	{
		thisHisto->Fill( gammaL_tau_value->at(i), gammaH_tau_value->at(i) );
	}

	TCanvas* c8 = EdStyle::RapidFitCanvas( "corr_Histo", "corr_Histo" );

	cout << "plot" << endl;
	thisHisto->Draw("colz");

	c8->Update();
	thisHisto->GetXaxis()->SetTitle("#Gamma_{L} [ps^{-1}]");
	thisHisto->GetYaxis()->SetTitle("#Gamma_{H} [ps^{-1}]");

	thisCorrGraph->Draw("SAMEL");
	thisCorrGraph2->Draw("SAMEL");

	c8->Update();
	c8->Print("corr2DHist.pdf");

	TFile* thisRoot = new TFile("corr2Doutput.root","RECREATE");

	thisHisto->Write("",TObject::kOverwrite);

	thisRoot->Write("",TObject::kOverwrite);

	//	See http://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall06/reading/gaussians.pdf 2.1 Bivariate Gaussian function
	//	http://www.maths.qmul.ac.uk/~ig/MTH5118/Notes11-09.pdf
	//
	TF2* corr_function = new TF2( "corr_matrix",
				"[5]/(2*TMath::Pi()*[0]*[1]*TMath::Sqrt(1 - [2]*[2])) * TMath::Exp( -(1./(2*(1-[2]*[2]))*(((x-[3])*(x-[3]))/([0]*[0])+((y-[4])*(y-[4]))/([1]*[1])-(2*[2]*(x-[3])*(y-[4]))/([0]*[1])) ))",
				xmin, xmax, ymin, ymax );

	corr_function->SetLineStyle( 1 );
	corr_function->SetLineWidth( 3 );

	corr_function->SetParameters( 0.01, 0.05, -0.8, 0.7215, 0.6215, 5e-04 );// xmin+((xmax-xmin)*0.5), ymin+((ymax-ymin)*0.5) );

	//corr_function->FixParameter( 0, 0.015 );
	//corr_function->FixParameter( 1, 0.03 );
	//corr_function->FixParameter( 2, -0.79 );
	//corr_function->FixParameter( 3, 0.7215 );//xmin+((xmax-xmin)*0.5)-0.005 );
	//corr_function->FixParameter( 4, 0.6215 );//ymin+((ymax-ymin)*0.5)-0.005 );

	//corr_function->FixParameter( 5, 5e-4 );

	int contour_num=3;
	//double* contours = new double[(unsigned)contour_num];
	//contours[0]=0.68;
	//contours[1]=0.9;
	//contours[2]=0.95;
	//corr_function->SetContour( contour_num, contours );

	gStyle->SetTitleOffset((Float_t)1.2,"X");
	gStyle->SetTitleOffset((Float_t)1.2,"Y");
	gROOT->UseCurrentStyle();
	gROOT->ForceStyle( true );


	TCanvas* c9 = EdStyle::RapidFitCanvas( "corr_func", "corr_func" );

	//corr_function->SetContour( 0, NULL );

	thisHisto->Fit( corr_function, "N" );//"EN" );

	double Integral = corr_function->Integral( xmin, xmax, ymin, ymax, 1E-12 );

	/*double* contours = new double[(unsigned)contour_num];
	  contours[0]=0.68*Integral;
	  contours[1]=0.9*Integral;
	  contours[2]=0.95*Integral;
	  corr_function->SetContour( contour_num, contours );*/

	corr_function->SetLineStyle( 1 );
	corr_function->SetLineWidth( 3 );

	thisHisto->GetXaxis()->SetNdivisions( 505, true );
	thisHisto->GetYaxis()->SetNdivisions( 505, true );

	thisHisto->Draw("");
	c9->Update();
	corr_function->SetLineStyle( 1 );
	int orig_contour_num = corr_function->GetContour();
	corr_function->SetContour( 10 );
	thisHisto->SetMarkerColor( 15 );
	c9->Update();

	corr_function->DrawCopy("SAME");

	c9->Update();

	c9->Print("corr_func_scatter.pdf");

	corr_function->SetContour( orig_contour_num );

	thisHisto->Draw("colz");
	c9->Update();

	corr_function->DrawCopy("SAME");

	c9->Update();

	c9->Print("corr_func.pdf");

	corr_function->SetLineWidth(3);

	thisHisto->GetXaxis()->SetTitleOffset((Float_t)1.2);
	thisHisto->GetYaxis()->SetTitleOffset((Float_t)1.2);

	thisHisto->Rebin2D();
	thisHisto->Rebin2D();
	thisHisto->Draw("lego2");
	
	corr_function->SetLineWidth( 3 );
	//corr_function->SetContour( 10 );
	corr_function->DrawCopy("CONT3 SAME LIST");
	c9->Update();
	corr_function->SetLineWidth( 3 );
	c9->Update();

	//TObjArray* contours= (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");

	TObjArray* contours2= (TObjArray*) gROOT->GetListOfSpecials();
	TList* contours = gPad->GetListOfPrimitives();

	cout << contours << endl;

	for( unsigned int i=0; i< (unsigned)contours->GetEntries(); ++i )
	{
		cout << i << "\t" << contours->At(i)->GetName() << endl;
		/*TList* thisContour = (TList*) contours->At(i);
		cout << thisContour << "\t" << i << endl;
		//	All closed contours in this program!
		TGraph* thisLine = (TGraph*) thisContour->First();
		cout << thisLine << endl;
		thisLine->SetLineWidth( 5 );*/
	}
	for( unsigned int i=0; i< (unsigned)contours2->GetEntries(); ++i )
	{
		cout << i << "\t" << contours2->At(i)->GetName() << endl;
	}
	//exit(0);
	c9->Update();

	c9->Print("corr_func3D.pdf");
	c9->Print("corr_func3D.C");

	TFile* corrfuncFile = new TFile("corrFunc.root","RECREATE");

	c9->Write("",TObject::kOverwrite);

	corrfuncFile->Write("",TObject::kOverwrite);
	//}

	return 0;
}

