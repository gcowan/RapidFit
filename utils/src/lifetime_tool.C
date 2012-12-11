
#include "EdStyle.h"
#include "TTree_Processing.h"
#include "ROOT_File_Processing.h"
#include "Histo_Processing.h"
#include "Template_Functions.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"

#include "TPaveStats.h"

#include <iomanip>
#include <vector>

using namespace::std;

int main( int argc, char* argv[] )
{
	if( argc != 3 )
	{
		cout << "Usage: " << argv[0] << " gammaL.root gammaH.root" << endl << endl;
	}

	//      Setup the Canvas and such
	EdStyle* RapidFit_Style = new EdStyle();
	RapidFit_Style->SetStyle();

	gStyle->SetOptFit(111);

	//	Open First TTree in the file
	TTree* gammaLTree = ROOT_File_Processing::GetFirstTree( argv[1] );

	vector<double>* gammaL_tau_value = TTree_Processing::Buffer_Branch( gammaLTree, "1./tau_value" );
	vector<double>* gammaL_tau_error = TTree_Processing::Buffer_Branch( gammaLTree, "tau_error/(tau_value*tau_value)" );
	vector<double>* gammaL_FitStatus = TTree_Processing::Buffer_Branch( gammaLTree, "Fit_Status" );

	//	Open First TTree in the file
	TTree* gammaHTree = ROOT_File_Processing::GetFirstTree( argv[2] );

	vector<double>* gammaH_tau_value = TTree_Processing::Buffer_Branch( gammaHTree, "1./tau_value" );
	vector<double>* gammaH_tau_error = TTree_Processing::Buffer_Branch( gammaHTree, "tau_error/(tau_value*tau_value)" );
	vector<double>* gammaH_FitStatus = TTree_Processing::Buffer_Branch( gammaHTree, "Fit_Status" );

	TTree* gammaTree = ROOT_File_Processing::GetFirstTree( argv[3] );

	vector<double>* gamma_value = TTree_Processing::Buffer_Branch( gammaTree, "gamma_value" );
	vector<double>* deltaGamma_value = TTree_Processing::Buffer_Branch( gammaTree, "deltaGamma_value" );
	vector<double>* FitStatus_value = TTree_Processing::Buffer_Branch( gammaTree, "Fit_Status" );

	vector<double> gammaS_value, gammaS_error, deltaGammaS_value, deltaGammaS_error;

	vector<double> gamma_diff, deltaGamma_diff;

	vector<double> angGammaTrudiff, untagGammaTrudiff, angdGammaTrudiff, untagdGammaTrudiff;

	TGraphErrors* corr_graph = new TGraphErrors( gammaH_tau_value->size(), &((*gammaL_tau_value)[0]), &((*gammaH_tau_value)[0]),
			&((*gammaL_tau_error)[0]), &((*gammaH_tau_error)[0]) );

	double GL_sum=0., GH_sum=0.;
	double GL_sum_sq=0., GH_sum_sq=0.;
	double GLGH_prod_sum=0.;
	double GL_avr=0., GH_avr=0.;
	double N = (double)(gammaL_tau_value->size());

	for( unsigned int i=0; i< gammaL_tau_value->size(); ++i )
	{
		GL_sum += (*gammaL_tau_value)[i];
		GH_sum += (*gammaH_tau_value)[i];
		GL_sum_sq += (*gammaL_tau_value)[i]*(*gammaL_tau_value)[i];
		GH_sum_sq += (*gammaH_tau_value)[i]*(*gammaH_tau_value)[i];
		GLGH_prod_sum += (*gammaL_tau_value)[i]*(*gammaH_tau_value)[i];
	}
	GL_avr = GL_sum / (double)(gammaL_tau_value->size());
	GH_avr = GH_sum / (double)(gammaH_tau_value->size());

	double corr = N*GLGH_prod_sum - GL_sum*GH_sum;
	corr /= sqrt( (N*GL_sum_sq - GL_sum*GL_sum) * (N*GH_sum_sq - GH_sum*GH_sum) );

	cout << endl << setprecision(10) << "GL / GH Correlation: " << corr << endl << endl;

	for( unsigned int i=0; i< gammaL_tau_value->size(); ++i )
	{
		double gammaS_val=0., gammaS_err=0., delGammaS_val=0., delGammaS_err=0.;

		gammaS_val = 0.5 * ( (*gammaL_tau_value)[i] + (*gammaH_tau_value)[i] );
		gammaS_err = 0.5 * ( sqrt( (*gammaL_tau_error)[i]*(*gammaL_tau_error)[i]+(*gammaH_tau_error)[i]*(*gammaH_tau_error)[i] ) );

		delGammaS_val = (*gammaL_tau_value)[i] - (*gammaH_tau_value)[i];
		delGammaS_err = sqrt( (*gammaL_tau_error)[i]*(*gammaL_tau_error)[i]+(*gammaH_tau_error)[i]*(*gammaH_tau_error)[i] );

		gammaS_value.push_back( gammaS_val );
		gammaS_error.push_back( gammaS_err );

		deltaGammaS_value.push_back( delGammaS_val );
		deltaGammaS_error.push_back( delGammaS_err );

		if( (*gammaH_FitStatus)[i] == 3 )
		{
			if( (*gammaL_FitStatus)[i] == 3 )
			{
				if( (*FitStatus_value)[i] == 3 )
				{

					gamma_diff.push_back( gammaS_val - (*gamma_value)[i] );

					angGammaTrudiff.push_back( gammaS_val );
					untagGammaTrudiff.push_back( (*gamma_value)[i] );

					deltaGamma_diff.push_back( delGammaS_val - (*deltaGamma_value)[i] );

					angdGammaTrudiff.push_back( delGammaS_val );
					untagdGammaTrudiff.push_back( (*deltaGamma_value)[i] );
				}
			}
		}
	}

	TH1* gamma_histo = new TH1D( "gamma_dist", "gamma_dist", gammaS_value.size(), get_minimum(gammaS_value), get_maximum(gammaS_value) );
	for( unsigned int i=0; i< gammaS_value.size(); ++i ) gamma_histo->Fill( gammaS_value[i] );

	TH1* dGamma_histo = new TH1D( "deltaGamma_dist", "deltaGamma_dist", deltaGammaS_value.size(), get_minimum(deltaGammaS_value), get_maximum(deltaGammaS_value) );
	for( unsigned int i=0; i< gammaS_value.size(); ++i ) dGamma_histo->Fill( deltaGammaS_value[i] );

	TH1* delta_gamma_histo = new TH1D( "gammadiff_dist", "gammadiff_dist", gamma_diff.size(), get_minimum(gamma_diff), get_maximum(gamma_diff) );
	for( unsigned int i=0; i< gamma_diff.size(); ++i ) delta_gamma_histo->Fill( gamma_diff[i] );

	TH1* delta_dGamma_histo = new TH1D( "dGammadiff_dist", "dGammadiff_dist", deltaGamma_diff.size(), get_minimum(deltaGamma_diff), get_maximum(deltaGamma_diff) );
	for( unsigned int i=0; i< gamma_diff.size(); ++i ) delta_dGamma_histo->Fill( deltaGamma_diff[i] );

	TH1* gamma_truth =new TH1D( "gamma_tru", "gamma_tru", angGammaTrudiff.size(), get_minimum(angGammaTrudiff), get_maximum(angGammaTrudiff) );
	for( unsigned int i=0; i< angGammaTrudiff.size(); ++i ) gamma_truth->Fill( angGammaTrudiff[i] );
	Histogram_Processing::OptimallyRebin( gamma_truth ); gamma_truth->SetLineColor( 2 );
	TH1* gamma_untag =new TH1D( "gamma_untag", "gamma_untag", untagGammaTrudiff.size(), get_minimum(untagGammaTrudiff), get_maximum(untagGammaTrudiff) );
	for( unsigned int i=0; i< untagGammaTrudiff.size(); ++i ) gamma_untag->Fill( untagGammaTrudiff[i] );
	Histogram_Processing::OptimallyRebin( gamma_untag ); gamma_untag->SetLineColor( 3 );

	TH1* dGamma_truth = new TH1D( "deltaGamma_tru", "deltaGamma_tru", angdGammaTrudiff.size(), get_minimum(angdGammaTrudiff), get_maximum(angdGammaTrudiff) );
	for( unsigned int i=0; i< angdGammaTrudiff.size(); ++i ) dGamma_truth->Fill( angdGammaTrudiff[i] );
	Histogram_Processing::OptimallyRebin( dGamma_truth ); dGamma_truth->SetLineColor( 2 );
	TH1* dGamma_untag =new TH1D( "deltaGamma_untag", "deltaGamma_untag", untagdGammaTrudiff.size(), get_minimum(untagdGammaTrudiff), get_maximum(untagdGammaTrudiff) );
	for( unsigned int i=0; i< untagdGammaTrudiff.size(); ++i ) dGamma_untag->Fill( untagdGammaTrudiff[i] );
	Histogram_Processing::OptimallyRebin( dGamma_untag ); dGamma_untag->SetLineColor( 3 );

	TCanvas* c1 = new TCanvas( "gamma_C", "gamma_C", 1680, 1050 );

	gamma_histo->Draw();
	c1->Update();
	Histogram_Processing::OptimallyRebin( gamma_histo );
	TString fit_type = Histogram_Processing::Best_Fit_Function( gamma_histo );
	Histogram_Processing::Silent_Fit( gamma_histo, fit_type );
	TPaveStats* thisStats = (TPaveStats*)gamma_histo->GetListOfFunctions()->FindObject("stats");
	thisStats->SetFillStyle(3023);
	thisStats->SetTextColor(1);
	thisStats->Draw();
	gamma_histo->GetXaxis()->SetTitle("#Gamma [ps^{-1}]");
	gamma_histo->GetYaxis()->SetTitle("Candidates");
	c1->Update();
	c1->Print("gammaOutput.pdf");

	TCanvas* c2 = new TCanvas( "dgamma_C", "dgamma_C", 1680, 1050 );

	dGamma_histo->Draw();
	c2->Update();
	Histogram_Processing::OptimallyRebin( dGamma_histo );
	fit_type = Histogram_Processing::Best_Fit_Function( dGamma_histo );
	Histogram_Processing::Silent_Fit( dGamma_histo, fit_type );
	thisStats = (TPaveStats*)dGamma_histo->GetListOfFunctions()->FindObject("stats");
	thisStats->SetFillStyle(3023);
	thisStats->SetTextColor(1);
	dGamma_histo->GetXaxis()->SetTitle("#Delta#Gamma [ps^{-1}]");
	dGamma_histo->GetYaxis()->SetTitle("Candidates");
	thisStats->Draw();
	c2->Update();
	c2->Print("deltaGammaOutput.pdf");

	TCanvas* c3 = new TCanvas( "diffgamma_C", "diffgamma_C", 1680, 1050 );

	delta_gamma_histo->Draw();
	c3->Update();
	Histogram_Processing::OptimallyRebin( delta_gamma_histo );
	fit_type = Histogram_Processing::Best_Fit_Function( delta_gamma_histo );
	Histogram_Processing::Silent_Fit( delta_gamma_histo, fit_type );
	thisStats = (TPaveStats*)delta_gamma_histo->GetListOfFunctions()->FindObject("stats");
	thisStats->SetFillStyle(3023);
	thisStats->SetTextColor(1);
	thisStats->Draw();
	delta_gamma_histo->GetXaxis()->SetTitle("#delta(#Gamma) [ps^{-1}]");
	delta_gamma_histo->GetYaxis()->SetTitle("Candidates");
	c3->Update();
	c3->Print("diff_Gamma.pdf");

	TCanvas* c4 = new TCanvas( "diffdGamma_C", "diffdGamma_C", 1680, 1050 );

	delta_dGamma_histo->Draw();
	c4->Update();
	Histogram_Processing::OptimallyRebin( delta_dGamma_histo );
	fit_type = Histogram_Processing::Best_Fit_Function( delta_dGamma_histo );
	Histogram_Processing::Silent_Fit( delta_dGamma_histo, fit_type );
	thisStats = (TPaveStats*)delta_dGamma_histo->GetListOfFunctions()->FindObject("stats");
	thisStats->SetFillStyle(3023);
	thisStats->SetTextColor(1);
	delta_dGamma_histo->GetXaxis()->SetTitle("#delta(#Delta#Gamma) [ps^{-1}]");
	delta_dGamma_histo->GetYaxis()->SetTitle("Candidates");
	thisStats->Draw();
	c4->Update();
	c4->Print("diff_deltaGamma.pdf");

	TCanvas* c5 = new TCanvas( "gamma_overlay", "gamma_overlay", 1680, 1050 );

	gamma_truth->Draw();
	gamma_truth->GetXaxis()->SetTitle("#Gamma [ps^{-1}]");
	gamma_truth->GetYaxis()->SetTitle("Candidates");
	gamma_untag->Draw("SAME");
	c5->Update();
	c5->Print("gamma_overlay.pdf");

	TCanvas* c6 = new TCanvas( "dGamma_overlay", "dGamma_overlay", 1680, 1050 );

	dGamma_truth->Draw();
	dGamma_truth->GetXaxis()->SetTitle("#Delta#Gamma [ps^{-1}]");
	dGamma_truth->GetYaxis()->SetTitle("Candidates");
	dGamma_untag->Draw("SAME");
	c6->Update();
	c6->Print("deltaGamma_overlay.pdf");

	TCanvas* c7 = new TCanvas( "corr_Graph", "corr_Graph", 1680, 1050 );

	corr_graph->Draw("AP");
	c7->Update();
	corr_graph->GetXaxis()->SetTitle("#Gamma_{L}");
	corr_graph->GetYaxis()->SetTitle("#Gamma_{H}");
	c7->Update();

	TString funcStr_Orig("("); funcStr_Orig+=corr;funcStr_Orig.Append("*x+[0])");

	TString funcStr("([1]*x+[0])");

	double xmin = corr_graph->GetXaxis()->GetXmin();
	double xmax = corr_graph->GetXaxis()->GetXmax();
	double ymin = corr_graph->GetYaxis()->GetXmin();
	double ymax = corr_graph->GetYaxis()->GetXmax();

	TF1* corr_funct = new TF1( "corr_func", funcStr, xmin, xmax );

	corr_funct->SetParameters( 1.18579, corr );

	corr_graph->Fit( corr_funct, "R" );

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

	corr_graph->Fit( corr_funct2, "R" );

	double* x_func2 = new double[100];
	double* y_func2 = new double[100];

	for( unsigned int i=0; i< 100; ++i )
	{
		x_func2[i] = xmin + xstep * i;
		y_func2[i] = corr_funct2->Eval( x_func2[i] );
	}

	TGraph* thisCorrGraph2 = new TGraph( 100, x_func2, y_func2 );

	thisCorrGraph2->SetLineColor(2);
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

	int xbins = Histogram_Processing::OptimumBinNumber( thisHisto );
	int ybins = Histogram_Processing::OptimumBinNumber( thisHisto, 2 );

	thisHisto = new TH2D( "histo_colz", "histo_colz", xbins, xmin, xmax, ybins, ymin, ymax );

	for( unsigned int i=0; i< gammaL_tau_value->size(); ++i )
	{
		thisHisto->Fill( gammaL_tau_value->at(i), gammaH_tau_value->at(i) );
	}

	TCanvas* c8 = new TCanvas( "corr_Histo", "corr_Histo", 1680, 1050 );

	cout << "plot" << endl;
	thisHisto->Draw("colz");

	c8->Update();

	thisCorrGraph->Draw("SAMEL");
	thisCorrGraph2->Draw("SAMEL");

	c8->Update();
	c8->Print("corr2DHist.pdf");

	return 0;
}

