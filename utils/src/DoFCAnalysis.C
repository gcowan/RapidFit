//	ROOT Headers
#include "TTree.h"
#include "TRandom3.h"
#include "TString.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMarker.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TObject.h"
#include "TLegend.h"
#include "THStack.h"
//	RapidFit Headers
#include "EdStyle.h"
//	utils Headers
#include "DoFCAnalysis.h"
#include "TTree_Processing.h"
#include "Histo_Processing.h"
#include "Template_Functions.h"
#include "RapidFit_Output_File.h"
#include "OutputPlots.h"
//	System Headers
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace::std;

//	Print some info about the toys at this coordinate
void FeldmanCousinsAnalysis::Print_Info_At_Grid_ij( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, vector<vector<Double_t> > this_Grid_Coordinate )
{
	cout << "AT: ";
	vector<string>::iterator param_i = controlled_parameter_name.begin();
	vector<Double_t>::iterator grid_ij = coordinate.begin();
	for( ; param_i != controlled_parameter_name.end(); ++param_i )
	{
		cout << *param_i << ": " << *grid_ij << " ";
	}
	cout << "\t" << "THERE ARE:\t" << this_Grid_Coordinate.size() << " TOYS!";// << endl;
}

//	Print some useful output about the number of toys at this coordinate
void FeldmanCousinsAnalysis::Print_Toy_Info( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, int Fixed, int Free )
{
	cout << "AT: ";
	vector<string>::iterator param_i = controlled_parameter_name.begin();
	vector<Double_t>::iterator grid_ij = coordinate.begin();
	for( ; param_i != controlled_parameter_name.end(); ++param_i )
	{
		cout << *param_i << ": " << *grid_ij << " ";
	}
	cout << "\t" << "THERE ARE:\t" << Fixed << " Fixed Toys and " << Free << " Free Toys.";// << endl;
}

//	Construct a string which describes this unique coordinate with no spaces
TString FeldmanCousinsAnalysis::Flat_Coord_String( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate )
{
	vector<string>::iterator param_i = controlled_parameter_name.begin();
	vector<Double_t>::iterator grid_ij = coordinate.begin();
	TString name;
	for( ; param_i != controlled_parameter_name.end(); ++param_i, ++grid_ij )
	{
		stringstream val_stream;
		val_stream << setprecision(4) << *grid_ij;
		TString val( val_stream.str() );
		name.Append("_"+TString(*param_i)+"_"+val);
	}
	return name;
}

//	Construct multiple plots for the nuicence parameter distributions at this point
void FeldmanCousinsAnalysis::Plot_Nuisance_Parameters( TTree* input_tree, TString& free_toys, TString& fixed_toys, vector<TString>& Free_Params,
							vector<string> controlled_parameter_name, vector<Double_t> coordinate, TRandom3* rand, double DLL_Data )
{
	TString grid_name = Flat_Coord_String( controlled_parameter_name, coordinate );
	vector<TString> param_var;
	vector<TString> Free_Params_copy = Free_Params;
	for( vector<string>::iterator cpar_i = controlled_parameter_name.begin(); cpar_i != controlled_parameter_name.end(); ++cpar_i )
	{
		Free_Params_copy.push_back( cpar_i->c_str() );
	}
	param_var.push_back( "_value" );
	param_var.push_back( "_error" );
	param_var.push_back( "_pull" );
	for( vector<TString>::iterator Param_i = Free_Params_copy.begin(); Param_i != Free_Params_copy.end(); ++Param_i )
	{
		vector<TString>::iterator var_i = param_var.begin();
		vector<TString>::iterator var_e = param_var.end();
		for( ; var_i != var_e; ++var_i )
		{
			TString rand_str; rand_str+=rand->Rndm();
			string rand_string( rand_str.Data() );
			replace( rand_string.begin(), rand_string.end(), '.', '_' );
			rand_str = rand_string.c_str();
			//vector<vector<Double_t> > Free_this_Grid_Coordinate = Plotter_Data( input_tree, *Param_i+*var_i, free_toys, rand );
			//vector<vector<Double_t> > Fixed_this_Grid_Coordinate = Plotter_Data( input_tree, *Param_i+*var_i, fixed_toys, rand );

			vector<Double_t>* param_data_free = TTree_Processing::Buffer_Branch( input_tree, *Param_i+*var_i, free_toys ); //Free_this_Grid_Coordinate[0];
			vector<Double_t>* nll_free = TTree_Processing::Buffer_Branch( input_tree, "NLL", free_toys );
			/*vector<Double_t>* nll_free_Asr = TTree_Processing::Buffer_Branch( input_tree, "As_sq_error", free_toys );
			  vector<Double_t>* nll_free_dpa = TTree_Processing::Buffer_Branch( input_tree, "delta_para_value", free_toys );
			  vector<Double_t>* nll_free_dpar = TTree_Processing::Buffer_Branch( input_tree, "delta_para_error", free_toys );
			  vector<Double_t>* nll_free_dpe = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_value*(delta_perp_value>=0.)*(delta_perp_value<=(2*3.14159))", free_toys );
			  vector<Double_t>* nll_free_dper = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_error", free_toys );
			  vector<Double_t>* nll_free_ds = TTree_Processing::Buffer_Branch( input_tree, "delta_s_value*(delta_s_value>=0.)*(delta_s_value<=(2*3.14159))", free_toys );
			  vector<Double_t>* nll_free_dsr = TTree_Processing::Buffer_Branch( input_tree, "delta_s_error", free_toys );
			  vector<Double_t>* nll_free_game = TTree_Processing::Buffer_Branch( input_tree, "gamma_error", free_toys );
			  vector<Double_t>* nll_free_dmr = TTree_Processing::Buffer_Branch( input_tree, "deltaM_error", free_toys ); */
			//vector<Double_t>* nll_free_phi = TTree_Processing::Buffer_Branch( input_tree, "Phi_s_value", free_toys );
			//vector<Double_t>* nll_free_phr = TTree_Processing::Buffer_Branch( input_tree, "Phi_s_error", free_toys );
			//vector<Double_t>* nll_free_p1e = TTree_Processing::Buffer_Branch( input_tree, "mistagP1_error", free_toys );
			//vector<Double_t>* nll_free_p0e = TTree_Processing::Buffer_Branch( input_tree, "mistagP0_error", free_toys );

			vector<Double_t>* param_data_fixed = TTree_Processing::Buffer_Branch( input_tree, *Param_i+*var_i, fixed_toys ); //Fixed_this_Grid_Coordinate[0];
			/*vector<Double_t>* nll_fixed_Asr = TTree_Processing::Buffer_Branch( input_tree, "As_sq_error", fixed_toys );
			  vector<Double_t>* nll_fixed_dpa = TTree_Processing::Buffer_Branch( input_tree, "delta_para_value", fixed_toys );
			  vector<Double_t>* nll_fixed_dpar = TTree_Processing::Buffer_Branch( input_tree, "delta_para_error", fixed_toys );
			  vector<Double_t>* nll_fixed_dpe = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_value*(delta_perp_value>=0.)*(delta_perp_value<=(2*3.14159))", fixed_toys );
			  vector<Double_t>* nll_fixed_dper = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_error", fixed_toys );*/
			vector<Double_t>* nll_fixed = TTree_Processing::Buffer_Branch( input_tree, "NLL", fixed_toys );
			/*vector<Double_t>* nll_fixed_ds = TTree_Processing::Buffer_Branch( input_tree, "delta_s_value*(delta_s_value>=0.)*(delta_s_value<=(2*3.14159))", fixed_toys );
			  vector<Double_t>* nll_fixed_dsr = TTree_Processing::Buffer_Branch( input_tree, "delta_s_error", fixed_toys );
			  vector<Double_t>* nll_fixed_game = TTree_Processing::Buffer_Branch( input_tree, "gamma_error", fixed_toys );
			  vector<Double_t>* nll_fixed_dmr = TTree_Processing::Buffer_Branch( input_tree, "deltaM_error", fixed_toys );
			//vector<Double_t>* nll_fixed_phi = TTree_Processing::Buffer_Branch( input_tree, "Phi_s_value", fixed_toys );
			vector<Double_t>* nll_fixed_phr = TTree_Processing::Buffer_Branch( input_tree, "Phi_s_error", fixed_toys );
			//vector<Double_t>* nll_fixed_p1e = TTree_Processing::Buffer_Branch( input_tree, "mistagP1_error", fixed_toys );
			//vector<Double_t>* nll_fixed_p0e = TTree_Processing::Buffer_Branch( input_tree, "mistagP0_error", fixed_toys );

			vector<Double_t>* Free_Toy_diff = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_value-delta_s_value", free_toys );
			vector<Double_t>* Fixed_Toy_diff = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_value-delta_s_value", fixed_toys );

			 */

			bool bad = false;

			vector<Double_t> free_2_plot, fixed_2_plot;

			for( unsigned int i=0; i != param_data_free->size(); ++i )
			{
				double strong_upper = 2.5, strong_lower = 1.;

				//if( (*nll_free_As)[i] > 0.012 && (*nll_fixed_As)[i] > 0.012 )
				//{
				//if( (*nll_free_dpa)[i] < strong_upper && (*nll_fixed_dpa)[i] < strong_upper )
				//{
				//if( fabs((*nll_free_phi)[i]) < 1 && fabs((*nll_fixed_phi)[i]) < 1 )
				//{
				//if( (*nll_free_dpa)[i] > 0 && (*nll_fixed_dpa)[i] > 0 )
				//{
				//if( fabs((*Free_Toy_diff)[i]) < 1. && fabs((*Fixed_Toy_diff)[i]) < 1. )
				//{
				/*if( (*nll_free_dpe)[i] > strong_lower && (*nll_fixed_dpe)[i] > strong_lower )
				  {
				  if( (*nll_free_dpe)[i] < strong_upper && (*nll_fixed_dpe)[i] < strong_upper )
				  {
				  if( (*nll_free_ds)[i] > strong_lower && (*nll_fixed_ds)[i] > strong_lower )
				  {
				  if( (*nll_free_ds)[i] < strong_upper && (*nll_fixed_ds)[i] < strong_upper )
				  {*/
				/*
				   if( (*nll_free_game)[i] < 0.00575 && (*nll_fixed_game)[i] < 0.00575 )
				   {
				   if( (*nll_free_dpar)[i] > 0.05 && (*nll_fixed_dpar)[i] > 0.05 )
				   {
				   if( (*nll_free_dpar)[i] < 0.35 && (*nll_fixed_dpar)[i] < 0.35 )
				   {
				   if( (*nll_free_dper)[i] > 0.15 && (*nll_fixed_dper)[i] > 0.15 )
				   {
				   if( (*nll_free_dper)[i] < 0.4 && (*nll_fixed_dper)[i] < 0.4 )
				   {
				   if( (*nll_free_dsr)[i] < 0.35 && (*nll_fixed_dsr)[i] < 0.35 )
				   {
				   if( (*nll_free_dsr)[i] > 0.15 && (*nll_fixed_dsr)[i] > 0.15 )
				   {
				   if( (*nll_free_Asr)[i] < 0.0125 && (*nll_fixed_Asr)[i] < 0.0125 )
				   {
				   if( (*nll_free_Asr)[i] > 0.0075 && (*nll_fixed_Asr)[i] > 0.0075 )
				   {
				   if( (*nll_free_dmr)[i] < 0.125 && (*nll_fixed_dmr)[i] < 0.125 )
				   {
				   if( (*nll_free_phr)[i] < 0.15 && (*nll_fixed_phr)[i] < 0.15 )
				   {*//*
				if( (*nll_free_p1e)[i] < 0.02425 && (*nll_fixed_p1e)[i] < 0.02425 )
				{
				if( (*nll_free_p1e)[i] > 0.02375 && (*nll_fixed_p1e)[i] > 0.02375 )
				{
				if( (*nll_free_p0e)[i] < 0.0088 && (*nll_fixed_p1e)[i] < 0.0088 )
				{
				if( (*nll_free_p0e)[i] > 0.0084 && (*nll_fixed_p1e)[i] > 0.0084 )
				{
				if( (*nll_free_dpe)[i] < strong_upper && (*nll_fixed_dpe)[i] < strong_upper )
				{
				if( (*nll_free_ds)[i] < strong_upper && (*nll_fixed_ds)[i] < strong_upper )
				{
				if( (*nll_free_dpa)[i] > 2.5 && (*nll_fixed_dpa)[i] > 2.5 )
				{
				if( (*nll_free_dpe)[i] > strong_lower && (*nll_fixed_dpe)[i] > strong_lower )
				{
				if( (*nll_free_ds)[i] > strong_lower && (*nll_fixed_ds)[i] > strong_lower )
				{*/
				if( ((*nll_fixed)[i]-(*nll_free)[i]) < 0 )
				{
					bad = true;
					cout << "BAD:\t" << *Param_i+*var_i << "\tFREE: " << (*param_data_free)[i] << "\tFIXED: " << (*param_data_fixed)[i] ;
					cout << "\t\tNLL: " << (*nll_free)[i] << "\t"  << (*nll_fixed)[i] << "\tDLL(data): " << DLL_Data;
					cout << "\tDLL(toy): " << ((*nll_fixed)[i]-(*nll_free)[i]) << endl;
				}
				else if( ((*nll_fixed)[i]-(*nll_free)[i]) > 2.*DLL_Data )//&& ((*nll_free_As)[i] > 0.012 && (*nll_fixed_As)[i] > 0.012 ) )
				{
					bad = true;
					cout << "BAD2:\t" << *Param_i+*var_i << "\tFREE: " << (*param_data_free)[i] << "\tFIXED: " << (*param_data_fixed)[i] ;
					cout << "\t\tNLL: " << (*nll_free)[i] << "\t"  << (*nll_fixed)[i] << "\tDLL(data): " << DLL_Data;
					cout << "\tDLL(toy): " << ((*nll_fixed)[i]-(*nll_free)[i]) << endl;
				}
				free_2_plot.push_back( (*param_data_free)[i] );
				fixed_2_plot.push_back( (*param_data_fixed)[i] );
				//}}}}}}}}}}}
				//}
			}

			if( bad ) cout << endl;

			TCanvas* free_c1 = EdStyle::RapidFitCanvas( "TCanv_free_"+rand_str, "TCanv_free_"+rand_str );
			TH1* free_data_th1 = Histogram_Processing::Get_TH1( free_2_plot, rand );//, get_optimal_histo_bins(param_data_free) );
			free_data_th1->Draw();
			free_c1->Update();
			Histogram_Processing::Silent_Print( free_c1, *Param_i+"_"+*var_i+"_"+grid_name+"_free.pdf" );
			Histogram_Processing::Silent_Print( free_c1, *Param_i+"_"+*var_i+"_"+grid_name+"_free.png" );
			Histogram_Processing::Silent_Print( free_c1, *Param_i+"_"+*var_i+"_"+grid_name+"_free.C" );

			TCanvas* fixed_c1 = EdStyle::RapidFitCanvas( "TCanv_fixed_"+rand_str, "TCanv_fixed_"+rand_str );
			TH1* fixed_data_th1 = Histogram_Processing::Get_TH1( fixed_2_plot, rand );//, get_optimal_histo_bins(param_data_fixed) );
			fixed_data_th1->Draw();
			fixed_c1->Update();
			Histogram_Processing::Silent_Print( fixed_c1, *Param_i+"_"+*var_i+"_"+grid_name+"_fixed.pdf" );
			Histogram_Processing::Silent_Print( fixed_c1, *Param_i+"_"+*var_i+"_"+grid_name+"_fixed.png" );
			Histogram_Processing::Silent_Print( fixed_c1, *Param_i+"_"+*var_i+"_"+grid_name+"_fixed.C" );

			TGraph* free_graph = new TGraph( free_data_th1 );
			TGraph* fixed_graph = new TGraph( fixed_data_th1 );

			TMultiGraph* multi_graph = new TMultiGraph( "TMulti_"+rand_str, "TMulti_"+rand_str );

			multi_graph->Add( free_graph );
			multi_graph->Add( fixed_graph );

			TCanvas* c1 = EdStyle::RapidFitCanvas( "TCanv_"+rand_str, "TCanv_"+rand_str );

			multi_graph->Draw( "APL" );

			c1->Update();
			Histogram_Processing::Silent_Print( c1, *Param_i+"_"+*var_i+"_"+grid_name+"_multi.pdf" );
			Histogram_Processing::Silent_Print( c1, *Param_i+"_"+*var_i+"_"+grid_name+"_multi.png" );
			Histogram_Processing::Silent_Print( c1, *Param_i+"_"+*var_i+"_"+grid_name+"_multi.C" );
		}
	}
}

//	Construct and write out the nll distributions of the toys to root file and as plots
void FeldmanCousinsAnalysis::NLL_dists( vector<Double_t>& param_data_free, vector<Double_t>& param_data_fixed, vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, TRandom3* rand )
{
	TString grid_name = Flat_Coord_String( controlled_parameter_name, coordinate );
	TString rand_str; rand_str += rand->Rndm();
	string rand_string( rand_str.Data() );
	replace( rand_string.begin(), rand_string.end(), '.', '_' );
	rand_str = rand_string.c_str();
	TCanvas* free_c1 = EdStyle::RapidFitCanvas( "TCanv_free_"+rand_str, "TCanv_free_"+rand_str );
	TH1* free_data_th1 = Histogram_Processing::Get_TH1( param_data_free, rand );//, get_optimal_histo_bins(param_data_free) );
	free_data_th1->Draw();
	free_c1->Update();
	Histogram_Processing::Silent_Print( free_c1, "NLL_"+grid_name+"_free.pdf" );
	Histogram_Processing::Silent_Print( free_c1, "NLL_"+grid_name+"_free.png" );
	Histogram_Processing::Silent_Print( free_c1, "NLL_"+grid_name+"_free.C" );

	TCanvas* fixed_c1 = EdStyle::RapidFitCanvas( "TCanv_fixed_"+rand_str, "TCanv_fixed_"+rand_str );
	TH1* fixed_data_th1 = Histogram_Processing::Get_TH1( param_data_fixed, rand );//, get_optimal_histo_bins(param_data_fixed) );
	fixed_data_th1->Draw();
	fixed_c1->Update();
	Histogram_Processing::Silent_Print( fixed_c1, "NLL_"+grid_name+"_fixed.pdf" );
	Histogram_Processing::Silent_Print( fixed_c1, "NLL_"+grid_name+"_fixed.png" );
	Histogram_Processing::Silent_Print( fixed_c1, "NLL_"+grid_name+"_fixed.C" );

	TGraph* free_graph = new TGraph( free_data_th1 );
	TGraph* fixed_graph = new TGraph( fixed_data_th1 );

	TMultiGraph* multi_graph = new TMultiGraph( "TMulti_"+rand_str, "TMulti_"+rand_str );

	free_graph->SetLineColor(1);
	multi_graph->Add( free_graph );
	multi_graph->Add( fixed_graph );

	TCanvas* c1 = EdStyle::RapidFitCanvas( "TCanv_"+rand_str, "TCanv_"+rand_str );

	multi_graph->Draw("APL");

	c1->Update();
	Histogram_Processing::Silent_Print( c1, "NLL_"+grid_name+"_multi.pdf" );
	Histogram_Processing::Silent_Print( c1, "NLL_"+grid_name+"_multi.png" );
	Histogram_Processing::Silent_Print( c1, "NLL_"+grid_name+"_multi.C" );
}

//	Print some useful information to the screen&file about this grid point
void FeldmanCousinsAnalysis::Output_GridPoint( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, double LOCAL_DATA_DLL, vector<Double_t>& Toy_DLL_Dist )
{
	TString grid_name = Flat_Coord_String( controlled_parameter_name, coordinate );

	TFile* grid_point = new TFile( "File_"+grid_name+".root", "RECREATE" );

	double dll_min = 0.;
	double dll_max = 0.;

	for( vector<Double_t>::iterator toy_i = Toy_DLL_Dist.begin(); toy_i !=Toy_DLL_Dist.end(); ++toy_i )
	{
		if( *toy_i < dll_min ) dll_min = *toy_i;
		if( *toy_i > dll_max ) dll_max = *toy_i;
	}

	if( LOCAL_DATA_DLL > dll_max ) dll_max = LOCAL_DATA_DLL;

	TH1D* grid_th1 = new TH1D( "TH1_"+grid_name,"TH1_"+grid_name, 100, dll_min, dll_max );//&(Toy_DLL_Dist[0]) );

	for( unsigned int i=0; i< Toy_DLL_Dist.size(); ++i )
	{
		grid_th1->Fill( Toy_DLL_Dist[i] );
	}

	TCanvas* c1 = EdStyle::RapidFitCanvas( "TCanvas_"+grid_name, "TCanvas_"+grid_name );

	grid_th1->Draw("PE");
	c1->Update();           //      STUPID ROOT

	double y_max = grid_th1->GetBinContent( grid_th1->GetMaximumBin() );

	TLine* data_DLL = new TLine( (Double_t)LOCAL_DATA_DLL, 0., (Double_t)LOCAL_DATA_DLL, (Double_t)y_max );

	data_DLL->Draw( "SAME" );
	c1->Update();

	c1->SetLogy();
	c1->Update();

	Histogram_Processing::Silent_Print( c1, "DLL_Dist"+grid_name+".pdf");
	Histogram_Processing::Silent_Print( c1, "DLL_Dist"+grid_name+".png");
	Histogram_Processing::Silent_Print( c1, "DLL_Dist"+grid_name+".C" );

	grid_th1->Write("",TObject::kOverwrite);
	c1->Write("",TObject::kOverwrite);

	grid_point->Close();
}

//	Return a pair of vectors first=CV of each control parameter second=CV_err of each control parameter
pair<vector<double>*, vector<double>* > FeldmanCousinsAnalysis::Get_Best_CV( TTree* input_tree, vector<string>& controlled_param_name )
{
	//      By definition the first event in a file is the Global Best fit hence can safely extract the NLL
	input_tree->Draw("NLL","","goff",1,0);
	double GLOBAL_BEST_NLL = input_tree->GetV1()[0];

	TString NLL_Cut_String = "(abs(NLL-"; NLL_Cut_String+=GLOBAL_BEST_NLL; NLL_Cut_String.Append( ")<"+double_tolerance+")" );

	vector<double>* Global_CV=new vector<double>();
	vector<double>* Global_CV_err = new vector<double>();

	for( vector<string>::iterator param_i = controlled_param_name.begin(); param_i != controlled_param_name.end(); ++param_i )
	{
		TString param=param_i->c_str();
		input_tree->Draw( param+"_value", NLL_Cut_String, "goff", 1, 0 );
		double param_val = input_tree->GetV1()[0];
		input_tree->Draw( param+"_error", NLL_Cut_String, "goff", 1, 0 );
		double param_err = input_tree->GetV1()[0];
		Global_CV->push_back( param_val );
		Global_CV_err->push_back( param_err );
	}
	return make_pair( Global_CV, Global_CV_err );
}

//	Perform the full analysis of a FC output file from RapidFit
void FeldmanCousinsAnalysis::DoFCAnalysis( TTree* input_tree, vector<string>& controlled_parameter_name, TRandom3* rand, vector<string>& OtherOptions )
{
	cout << "Starting FC Analysis" << endl;
	print( controlled_parameter_name );

	(void) OtherOptions;

	//	Just to make sure we never lose any informatio
	input_tree->SetEstimate( input_tree->GetEntries() );

	//	DLL coordinates have the same _gen as the global CV, but all have unique fixed _value
	vector<TString> Free_Params = RapidFit_Output_File::get_free_non_scanned_parameters( input_tree, controlled_parameter_name );

	//	This Draw String gives the returned control parameters from the cut string
	TString Grid_Draw_String = RapidFit_Output_File::Construct_Draw_String( controlled_parameter_name );

	//	This Cut String returns only the fit results which are 
	TString Grid_Cut_String = RapidFit_Output_File::Construct_Cut_String( input_tree, controlled_parameter_name, true );

	//	Output for the user
	cout << "Draw:\t" << Grid_Draw_String << endl;
	cout << "Cut: \t" << Grid_Cut_String << endl;

	//	By definition the first event in a file is the Global Best fit hence can safely extract the NLL
	input_tree->Draw("NLL","","goff",1,0);
	double GLOBAL_BEST_NLL = input_tree->GetV1()[0];

	vector<double>* Global_CV=NULL;
	vector<double>* Global_CV_err=NULL;

	pair<vector<double>*, vector<double>* > Global_Params = Get_Best_CV( input_tree, controlled_parameter_name );

	Global_CV = Global_Params.first; Global_CV_err = Global_Params.second;

	//	The results from the Plotter_Data are BY DEFINITION mangled, but they are a unique sorted set of grid_coordinates

	//      return the *UNIQUE* corrdinates contained in the input_tree from the Draw_String after applying the Cut_String
	vector<vector<Double_t> > Grid_Coordinates = TTree_Processing::Plotter_Data( input_tree, Grid_Draw_String, Grid_Cut_String, rand );
	Grid_Coordinates = rotate( Grid_Coordinates );

	//	Output for the user
	cout << "DIMENTION: " << Grid_Coordinates[0].size() << "\t" << "COORDINATES: " << Grid_Coordinates.size() << endl;

	//	Vectors to store the output from the analysis
	vector<double> DATA_DLL;
	vector<pair<vector<double>, double> > ALL_CL_FROM_FC;

	//	Make an output Directory for the FC analysis tool
	gSystem->mkdir("RapidFit_FC_distributions");
	gSystem->cd("RapidFit_FC_distributions");

	TDirectory* here = gDirectory;

	cout << endl;

	unsigned int coord=1;

	//vector<vector<Double_t> >::iterator grid_i = Grid_Coordinates.begin();
	//for( ; grid_i != Grid_Coordinates.end(); ++grid_i, ++coord )
	vector<vector<Double_t> >::iterator grid_i = Grid_Coordinates.end();
	for( --grid_i; grid_i != Grid_Coordinates.begin()-1; --grid_i, ++coord )
	{
		stringstream num; num << coord << "/" << Grid_Coordinates.size();
		cout << num.str() << setw(12);

		//	Cut String to select only a SINGLE fit result for the LL at this point
		TString data_cut = RapidFit_Output_File::Data_At_Grid_ij( controlled_parameter_name, *grid_i ); 

		//	Slect only toy Results which were generated at this grid point
		//	By definition the _gen value is equal to the fit minima from the LL scan result at this point
		//	Could even expand on this to use the _gen of all nuisence too, but thata extreme belt&braces
		TString toys_cut = RapidFit_Output_File::Fits_At_Grid_ij( controlled_parameter_name, *grid_i, true );

		//	Cut string which will return only free toys
		TString free_condition = RapidFit_Output_File::Construct_Fixed_condition( controlled_parameter_name, false );

		//	Cut string which will return only fixed toys
		TString fixed_condition = RapidFit_Output_File::Construct_Fixed_condition( controlled_parameter_name, true );

		//	Only Fixed/Free toys at this coordinate
		TString free_toys = RapidFit_Output_File::Fits_At_Grid_ij( controlled_parameter_name, *grid_i, true ) + "&&" + free_condition;// + "&&" + "(Phi_s_gen!=Phi_s_value)" + "&&(NLL>0)" + "&&(Phi_s_scan!=1)";
		//TString fixed_toys = toys_cut + "&&" + fixed_condition + "&&(NLL>0)" + "&&(Phi_s_scan!=1)";
		TString fixed_toys = RapidFit_Output_File::Fits_At_Grid_ij( controlled_parameter_name, *grid_i, true ) + "&&" + fixed_condition;// + "&&" + "(Phi_s_gen==Phi_s_value)" + "&&(NLL>0)" + "&&(Phi_s_scan!=1)";

		//	Some output for the user
		cout << "Free CUT:\t" << free_toys << endl;
		cout << "Fixed CUT:\t" << fixed_toys << endl;
		cout << "data_cut:\t" << data_cut << endl;

		//	this data will be mangled, however we only want the single unique result so don't care
		vector<vector<Double_t> > Data_Coordinate = TTree_Processing::Plotter_Data( input_tree, "NLL", data_cut, rand );

		//	Get the NLL for the fixed/free toys in a completely UNMANGLED way
		vector<Double_t>* Fixed_Toy_NLL = TTree_Processing::Buffer_Branch( input_tree, "NLL", fixed_toys );
		/*vector<Double_t>* Fixed_Toy_As = TTree_Processing::Buffer_Branch( input_tree, "As_sq_value", fixed_toys );
		  vector<Double_t>* Fixed_Toy_dpe = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_value*(delta_perp_value>=0.)*(delta_perp_value<=(2*3.14159))*(delta_s_value>=0.)*(delta_s_value<=(2*3.14159))", fixed_toys );
		  vector<Double_t>* Fixed_Toy_dper = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_error", fixed_toys );
		  vector<Double_t>* Fixed_Toy_dpa = TTree_Processing::Buffer_Branch( input_tree, "delta_para_value", fixed_toys );
		  vector<Double_t>* Fixed_Toy_dpar = TTree_Processing::Buffer_Branch( input_tree, "delta_para_error", fixed_toys );
		  vector<Double_t>* Fixed_Toy_dsp = TTree_Processing::Buffer_Branch( input_tree, "delta_s_value*(delta_s_value>=0.)*(delta_s_value<=(2*3.14159))*(delta_perp_value>=0.)*(delta_perp_value<=(2*3.14159))", fixed_toys );
		  vector<Double_t>* Fixed_Toy_dsr = TTree_Processing::Buffer_Branch( input_tree, "delta_s_error", fixed_toys );
		  vector<Double_t>* Fixed_Toy_phi = TTree_Processing::Buffer_Branch( input_tree, "Phi_s_value", fixed_toys );
		  vector<Double_t>* Fixed_Toy_phe = TTree_Processing::Buffer_Branch( input_tree, "Phi_s_error", fixed_toys );
		  vector<Double_t>* Fixed_Toy_gae = TTree_Processing::Buffer_Branch( input_tree, "gamma_error", fixed_toys );
		  vector<Double_t>* Fixed_Toy_Aser = TTree_Processing::Buffer_Branch( input_tree, "As_sq_error", fixed_toys );
		  vector<Double_t>* Fixed_Toy_trer = TTree_Processing::Buffer_Branch( input_tree, "timeResolutionScale_error", fixed_toys );
		  vector<Double_t>* Fixed_Toy_dMe = TTree_Processing::Buffer_Branch( input_tree, "deltaM_error", fixed_toys );
		  vector<Double_t>* Fixed_Toy_mp0 = TTree_Processing::Buffer_Branch( input_tree, "mistagP0_error", fixed_toys );
		  vector<Double_t>* Fixed_Toy_mp1 = TTree_Processing::Buffer_Branch( input_tree, "mistagP1_error", fixed_toys );
		  vector<Double_t>* Fixed_Toy_dGe = TTree_Processing::Buffer_Branch( input_tree, "deltaGamma_error", fixed_toys );
		 */

		vector<Double_t>* Free_Toy_NLL = TTree_Processing::Buffer_Branch( input_tree, "NLL", free_toys );

		/*vector<Double_t>* Free_Toy_As = TTree_Processing::Buffer_Branch( input_tree, "As_sq_value", free_toys );
		  vector<Double_t>* Free_Toy_dpe = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_value*(delta_perp_value>=0.)*(delta_perp_value<=(2*3.14159))*(delta_s_value>=0.)*(delta_s_value<=(2*3.14159))", free_toys );
		  vector<Double_t>* Free_Toy_dpa = TTree_Processing::Buffer_Branch( input_tree, "delta_para_value", free_toys );
		  vector<Double_t>* Free_Toy_dsp = TTree_Processing::Buffer_Branch( input_tree, "delta_s_value*(delta_s_value>=0.)*(delta_s_value<=(2*3.14159))*(delta_perp_value>=0.)*(delta_perp_value<=(2*3.14159))", free_toys );
		  vector<Double_t>* Free_Toy_phi = TTree_Processing::Buffer_Branch( input_tree, "Phi_s_value", free_toys );
		  vector<Double_t>* Free_Toy_phe = TTree_Processing::Buffer_Branch( input_tree, "Phi_s_error", free_toys );
		  vector<Double_t>* Free_Toy_gae = TTree_Processing::Buffer_Branch( input_tree, "gamma_error", free_toys );
		  vector<Double_t>* Free_Toy_dper = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_error", free_toys );
		  vector<Double_t>* Free_Toy_dpar = TTree_Processing::Buffer_Branch( input_tree, "delta_para_error", free_toys );
		  vector<Double_t>* Free_Toy_dsr = TTree_Processing::Buffer_Branch( input_tree, "delta_s_error", free_toys );
		  vector<Double_t>* Free_Toy_Aser = TTree_Processing::Buffer_Branch( input_tree, "As_sq_error", free_toys );
		  vector<Double_t>* Free_Toy_trer = TTree_Processing::Buffer_Branch( input_tree, "timeResolutionScale_error", free_toys );
		  vector<Double_t>* Free_Toy_dMe = TTree_Processing::Buffer_Branch( input_tree, "deltaM_error", free_toys );
		  vector<Double_t>* Free_Toy_mp0 = TTree_Processing::Buffer_Branch( input_tree, "mistagP0_error", free_toys );
		  vector<Double_t>* Free_Toy_mp1 = TTree_Processing::Buffer_Branch( input_tree, "mistagP1_error", free_toys );
		  vector<Double_t>* Free_Toy_dGe = TTree_Processing::Buffer_Branch( input_tree, "deltaGamma_error", free_toys );

		  vector<Double_t>* Free_Toy_diff = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_value-delta_s_value", free_toys );
		  vector<Double_t>* Fixed_Toy_diff = TTree_Processing::Buffer_Branch( input_tree, "delta_perp_value-delta_s_value", fixed_toys );

		  vector<Double_t>* Free_Toy_ds_cut = TTree_Processing::Buffer_Branch( input_tree, "(delta_s_value)*(abs(delta_perp_value-delta_s_value)<=1)", free_toys );
		  vector<Double_t>* Free_Toy_dp_cut = TTree_Processing::Buffer_Branch( input_tree, "(delta_perp_value)*(abs(delta_perp_value-delta_s_value)<=1)", free_toys );
		  vector<Double_t>* Fixed_Toy_ds_cut = TTree_Processing::Buffer_Branch( input_tree, "(delta_s_value)*(abs(delta_perp_value-delta_s_value)<=1)", fixed_toys );
		  vector<Double_t>* Fixed_Toy_dp_cut = TTree_Processing::Buffer_Branch( input_tree, "(delta_perp_value)*(abs(delta_perp_value-delta_s_value)<=1)", fixed_toys );
		 */
		gStyle->SetOptStat(0);
		TString grid_name = Flat_Coord_String( controlled_parameter_name, *grid_i );

		/*
		TString colz("colz"), e("");
		TH2* th1 = Histogram_Processing::Plot_2D( *Free_Toy_dsp, *Free_Toy_dpe, "scat_free_"+grid_name+".pdf", e, rand );
		TH2* th2 = Histogram_Processing::Plot_2D( *Free_Toy_dsp, *Free_Toy_dpe, "scat_free_"+grid_name+"_contz.pdf", colz, rand );

		TH2* th3 = Histogram_Processing::Plot_2D( *Fixed_Toy_dsp, *Fixed_Toy_dpe, "scat_fix_"+grid_name+".pdf", e, rand );

		TH2* th4 = Histogram_Processing::Plot_2D( *Fixed_Toy_dsp, *Fixed_Toy_dpe, "scat_fix_"+grid_name+"_contz.pdf", colz, rand );
		TH2* th5 = Histogram_Processing::Plot_2D( *Free_Toy_ds_cut, *Free_Toy_dp_cut, "scat_free_cut"+grid_name+".pdf", e, rand );
		TH2* th6 = Histogram_Processing::Plot_2D( *Free_Toy_ds_cut, *Free_Toy_dp_cut, "scat_free_cut_"+grid_name+"_contz.pdf", colz, rand );

		TH2* th7 = Histogram_Processing::Plot_2D( *Fixed_Toy_ds_cut, *Fixed_Toy_dp_cut, "scat_fix_cut_"+grid_name+".pdf", e, rand );
		TH2* th8 = Histogram_Processing::Plot_2D( *Fixed_Toy_ds_cut, *Fixed_Toy_dp_cut, "scat_fix_cut_"+grid_name+"_contz.pdf", colz, rand );

		TH1* diff_th1 = Histogram_Processing::Get_TH1( *Free_Toy_diff, rand );
		TH1* diff2_th1 = Histogram_Processing::Get_TH1( *Fixed_Toy_diff, rand );
		TString canv_nam_p("name");canv_nam_p+=rand->Rndm();
                string caonvas_string( canvas_nam_p.Data() );
                replace( canvas_string.begin(), canvas_string.end(), '.', '_' );
		canv_nam_p = canvas_string.c_str();
		here->cd();
		TCanvas*  canv_p = EdStyle::RapidFitCanvas( canv_nam_p, canv_nam_p );
		canv_p->SetLogy();
		diff_th1->Draw();
		diff2_th1->Draw("SAME");
		canv_p->Update();
		canv_p->Write("Overlay_diff_"+grid_name+".pdf");
		 */

		//	Draw the NLL distributions at this point
		//NLL_dists( *Free_Toy_NLL, *Fixed_Toy_NLL, controlled_parameter_name, *grid_i, rand );

		//	Print some information on the toys at this grid point
		//Print_Toy_Info( controlled_parameter_name, *grid_i, (int)Fixed_Toy_NLL->size(), (int)Free_Toy_NLL->size() );

		cout << "here" << endl;

		//      Record and output the DLL from data at this point
		cout << Data_Coordinate.size() << "\t" << Data_Coordinate[0].size() << endl;
		double LOCAL_DATA_DLL = Data_Coordinate[0][0] - GLOBAL_BEST_NLL;
		cout << Data_Coordinate[0][0] << " - " << GLOBAL_BEST_NLL << " = " << LOCAL_DATA_DLL << endl;
		DATA_DLL.push_back( LOCAL_DATA_DLL );

		//	Make some plots of the distributions of all free parameters fluctuating over the fit
		//Plot_Nuisance_Parameters( input_tree, free_toys, fixed_toys, Free_Params, controlled_parameter_name, *grid_i, rand, LOCAL_DATA_DLL );

		//	Check the total number of 'good' toys which we can use at this grid point
		unsigned int toys_to_test = (unsigned)Free_Toy_NLL->size();
		if( Free_Toy_NLL->size() != Fixed_Toy_NLL->size() )
		{
			cerr << "\tWARNING: DIFFERENT NUMBERS OF TOYS BETWEEN FIXED/FREE" << endl;
			cerr << Free_Toy_NLL->size() << "\t" << Fixed_Toy_NLL->size() << endl << endl;
			//	This if statement is to protect the CL from errors in the data
			//exit(-5);
			continue;
		}
		if( toys_to_test == 0 ) continue;

		cout << Free_Toy_NLL->size() << "\t" << Fixed_Toy_NLL->size() << endl << endl;

		//	Store the DLL from the toys in an array
		vector<double> Toy_DLL_Dist, other_dist, total;
		for( unsigned int i=0; i != toys_to_test; ++i )
		{
			bool added = false;
			//if( (*Fixed_Toy_As)[i] > 0.0012 && (*Free_Toy_As)[i] > 0.0012 )
			//{
			//if( fabs((*Free_Toy_dpa)[i]) < strong_up_lim && fabs((*Fixed_Toy_dpa)[i]) < strong_up_lim )
			//{
			/*
			   if( fabs((*Free_Toy_phi)[i]) < 1 && fabs((*Fixed_Toy_phi)[i]) < 1 )
			   {*/
			//if( fabs((*Free_Toy_diff)[i]) < 1. && fabs((*Fixed_Toy_diff)[i]) < 1. )
			//{
			//if( (*Free_Toy_As)[i] > 0.0012 && (*Fixed_Toy_As)[i] > 0.0012 )
			//{
			//if( (*Free_Toy_dsp)[i] > 0. && (*Fixed_Toy_dsp)[i] > 0. )
			//{
			//if( (*Free_Toy_dsp)[i] > strong_lo_lim && (*Fixed_Toy_dsp)[i] > strong_lo_lim )
			//{
			//if( (*Free_Toy_dpe)[i] > 0. && (*Fixed_Toy_dpe)[i] > 0. )
			//{
			//if( (*Free_Toy_dpa)[i] > 0 && (*Fixed_Toy_dpa)[i] > 0. )
			//{
			//if( (*Free_Toy_dpe)[i] > strong_lo_lim && (*Fixed_Toy_dpe)[i] > strong_lo_lim )
			//{
			//if( (*Free_Toy_dpa)[i] > 3.14159 && (*Fixed_Toy_dpa)[i] > 3.14159 )
			//{
			/*
			   if( fabs((*Free_Toy_phi)[i]) < 1.5 && fabs((*Fixed_Toy_phi)[i]) < 1.5 )
			   {*/
			/*if( (*Free_Toy_gae)[i] < 0.00575 && (*Fixed_Toy_gae)[i] < 0.00575  )
			  {
			  if( (*Fixed_Toy_dper)[i] > 0.15 && (*Free_Toy_dper)[i] > 0.15 )
			  {
			  if( (*Fixed_Toy_dper)[i] < 0.4 && (*Free_Toy_dper)[i] < 0.4 )
			  {
			  if( (*Fixed_Toy_dpar)[i] < 0.35 && (*Free_Toy_dpar)[i] < 0.35 )
			  {
			  if( (*Fixed_Toy_dpar)[i] > 0.05 && (*Free_Toy_dpar)[i] > 0.05 )
			  {
			  if( (*Fixed_Toy_dsr)[i] < 0.35 && fabs((*Free_Toy_dsr)[i]) < 0.35 )
			  {
			  if( (*Fixed_Toy_dsr)[i] > 0.15 && fabs((*Free_Toy_dsr)[i]) > 0.15 )
			  {
			  if( (*Fixed_Toy_Aser)[i] > 0.0075 && (*Free_Toy_Aser)[i] > 0.0075 )
			  {
			  if( (*Fixed_Toy_Aser)[i] < 0.0125 && (*Free_Toy_Aser)[i] < 0.0125 )
			  {
			  if( (*Fixed_Toy_dMe)[i] < 0.125 && (*Free_Toy_dMe)[i] < 0.125 )
			  {
			  if( (*Fixed_Toy_dMe)[i] > 0.05 && (*Free_Toy_dMe)[i] > 0.05 )
			  {				*/	/*
			   if( (*Fixed_Toy_phe)[i] < 0.15 && (*Free_Toy_phe)[i] < 0.15 )
			   {
			   if( (*Fixed_Toy_phe)[i] > 0.05 && (*Free_Toy_phe)[i] > 0.05 )
			   {
			   if( (*Fixed_Toy_trer)[i] < 0.04 && (*Free_Toy_trer)[i] < 0.04 )
			   {
			   if( (*Fixed_Toy_mp0)[i] < 0.0088 && (*Free_Toy_mp0)[i] < 0.0088 )
			   {
			   if( (*Fixed_Toy_mp0)[i] > 0.0084 && (*Free_Toy_mp0)[i] > 0.0084 )
			   {
			   if( (*Fixed_Toy_mp1)[i] < 0.02425 && (*Free_Toy_mp1)[i] < 0.02425 )
			   {
			   if( (*Fixed_Toy_mp1)[i] > 0.02375 && (*Free_Toy_mp1)[i] > 0.02375 )
			   {
			   if( fabs((*Fixed_Toy_dGe)[i])  <= 0.03 && fabs((*Free_Toy_dGe)[i]) <= 0.03 )
			   {*/
			//cout << Fixed_Toy_NLL[i] << "\t" << Free_Toy_NLL[i] << "\t\t" << Fixed_Toy_NLL[i] - Free_Toy_NLL[i] << endl;
			double toy_dll = (*Fixed_Toy_NLL)[i] - (*Free_Toy_NLL)[i];
			/*if( toy_dll < 0. )*/ Toy_DLL_Dist.push_back( toy_dll );
			added=true;
			//}
			//}}}}}}}}}}}}}}}}}}}}}}}}}
			//}}}}}}}}}}}//}}
			//}//}}}}}}}}}}//}

			if( !added )
			{
				other_dist.push_back( (*Fixed_Toy_NLL)[i] - (*Free_Toy_NLL)[i] );
			}
			total.push_back( (*Fixed_Toy_NLL)[i] - (*Free_Toy_NLL)[i] );
		}

		TH1* total_th1 = Histogram_Processing::Get_TH1( total, rand , 100, get_minimum(total), get_maximum(total) );
		TH1* acc_th1 = Histogram_Processing::Get_TH1( Toy_DLL_Dist, rand , 100, get_minimum(total), get_maximum(total) );
		//TH1* rej_th1 = Histogram_Processing::Get_TH1( other_dist, rand , 100, get_minimum(total), get_maximum(total) );

		total_th1->SetMarkerColor( 1 ); total_th1->SetLineColor( 1 );
		//total_th1->SetLineWidth( 4 );
		acc_th1->SetMarkerColor( 4 ); acc_th1->SetLineColor( 4 );
		//rej_th1->SetMarkerColor( 2 ); rej_th1->SetLineColor( 2 );

		//TGraph* total_graph = new TGraph( total_th1 );
		//TGraph* accepted_graph = new TGraph( acc_th1 ); accepted_graph->SetMarkerColor( 1 ); accepted_graph->SetLineColor( 1 );
		//TGraph* rejected_graph = new TGraph( rej_th1 ); rejected_graph->SetMarkerColor( 2 ); rejected_graph->SetLineColor( 2 );
	
		TString TMulti_Name("TMultiGraph_"); TMulti_Name+=rand->Rndm();
                string MultiStr( TMulti_Name.Data() );
                replace( MultiStr.begin(), MultiStr.end(), '.', '_' );
		TMulti_Name = MultiStr.c_str();
		//TMultiGraph* acc_rej = new TMultiGraph( TMulti_Name, TMulti_Name );
		TString m_canvas_name("TCanvas_"); m_canvas_name+=rand->Rndm();
                string CanvasStr( m_canvas_name.Data() );
                replace( CanvasStr.begin(), CanvasStr.end(), '.', '_' );
		m_canvas_name = CanvasStr.c_str();
		TCanvas* DLL_Overlay = EdStyle::RapidFitCanvas( m_canvas_name, m_canvas_name );
		DLL_Overlay->SetLogy();
		total_th1->Draw("lp9");
		DLL_Overlay->Update();
		acc_th1->Draw("SAME lp9");
		DLL_Overlay->Update();
		//rej_th1->Draw("SAME lp9");
		DLL_Overlay->Update();
		total_th1->GetXaxis()->SetTitle("DLL");
		total_th1->GetYaxis()->SetTitle("freq");
		TLegend* tleg = new TLegend( 0.75, 0.65, 0.9, 0.9 );
		tleg->AddEntry( total_th1, "total", "lp" );
		tleg->AddEntry( acc_th1, "accepted", "lp" );
		//tleg->AddEntry( rej_th1, "rejected", "lp" );
		TLine* data_l = new TLine( LOCAL_DATA_DLL, 0, LOCAL_DATA_DLL, total_th1->GetMaximum() );
		data_l->SetLineColor(3);
		data_l->Draw();
		tleg->AddEntry( data_l, "DLL from Data", "l" );
		tleg->Draw();
		DLL_Overlay->Update();
		DLL_Overlay->Print("DLL_Overlay"+grid_name+".pdf");
		DLL_Overlay->Print("DLL_Overlay"+grid_name+".png");
		DLL_Overlay->Print("DLL_Overlay"+grid_name+".C");

		//sort( Toy_DLL_Dist.begin(), Toy_DLL_Dist.end() );
		//print( Toy_DLL_Dist );


		//	Plot the output of the DLL from toys
		Output_GridPoint( controlled_parameter_name, *grid_i, (double)LOCAL_DATA_DLL, Toy_DLL_Dist );


		//	calculate the CL by comparing the DLL from toys to the DLL from data
		unsigned int toy_dll_smaller = 0;
		for( vector<double>::iterator dll_i = Toy_DLL_Dist.begin(); dll_i != Toy_DLL_Dist.end(); ++dll_i )
		{
			if( *dll_i < LOCAL_DATA_DLL ) ++toy_dll_smaller;
		}
		double LOCAL_CL_FROM_FC = double(toy_dll_smaller)/double(Toy_DLL_Dist.size());

		cout << "Data DLL: " << LOCAL_DATA_DLL << endl;
		cout << "\tC.L.:\t" << setprecision(10) << LOCAL_CL_FROM_FC << "\tNumber:\t" << Toy_DLL_Dist.size() << "\tEfficiency:\t" << double(double(Toy_DLL_Dist.size()) / double(Fixed_Toy_NLL->size())) << endl;

		//	Print some more output for this Grid Point
		//cout << "DATA: " << LOCAL_DATA_DLL << "\tTOY: " << LOCAL_CL_FROM_FC << endl;
		ALL_CL_FROM_FC.push_back( make_pair( *grid_i, LOCAL_CL_FROM_FC ) );

		//exit(1);
	}

	/*
	if( controlled_parameter_name.size() == 1 )
	{
		vector<vector<double> > FC_X_vals_temp = return_first( ALL_CL_FROM_FC );
		vector<double> FC_X_vals = request_element( FC_X_vals_temp, 0 );
		vector<double> FC_Y_vals = return_second( ALL_CL_FROM_FC );

		pair<vector<double>,vector<double> > FC_data = make_pair( FC_X_vals, FC_Y_vals );

		vector<OutputPlots*> _1D_Outputs = this->Plot_1D( reparam( FC_data ), DATA_DLL, (*Global_CV)[0], (*Global_CV_err)[0], rand );
	}
	else if( controlled_parameter_name.size() == 2 )
	{
		vector<vector<double> > DLL_X_vals = return_first( ALL_CL_FROM_FC );

		vector<pair<vector<double>, double> > DLL_Vals;

		vector<double>::iterator dll_val_i = DLL_Vals.begin();
		vector<vector<double> >::iterator dll_x_i = DLL_X_vals.begin();
		for( ; dll_x_i != DLL_X_vals.end(); ++dll_x_i, ++dll_val_i )
		{
			DLL_Vals.push_back( make_pair( (*dll_i), (*dll_val_i) ) );
		}

		vector<OutputPlots*> _2D_Outputs = this->Plot_2D( ALL_CL_FROM_FC, DLL_Vals, Global_CV, Global_CV_err, rand );
	}
	*/

	vector<double> X_data, X_data2, X_data3, X_1_data, X_1_data2, real_1_data, Y_data, Y_data2, Y_data3, Y_1_data, Y_1_data2;

	vector<double> sqrt_1_minus_Data_X, sqrt_1_minus_Data_Y, sqrt_1_minus_Theory_X, sqrt_1_minus_Theory_Y, sqrt_1_minus_FC_X, sqrt_1_minus_FC_Y;

	for( unsigned int i=0; i< ALL_CL_FROM_FC.size(); ++i )
	{
		cout << ALL_CL_FROM_FC[i].first[0] << "\t" << ALL_CL_FROM_FC[i].second << "\t" << DATA_DLL[i] << endl;

		double x =  ALL_CL_FROM_FC[i].first[0];
		double x_err = erf( fabs((x-(*Global_CV)[0])/(*Global_CV_err)[0])/sqrt(2.) );

		//	X_data and Y_data contain phi_s and FC-CL respectivley
		//	X_data2 and Y_data2 contain the theorectical CL based on the CV fit result
		X_data.push_back( x );
		X_data2.push_back( x );
		X_data3.push_back( x );
		Y_data.push_back( ALL_CL_FROM_FC[i].second );
		Y_data2.push_back( x_err );
		Y_data3.push_back( erf( sqrt(DATA_DLL[i]) ) );

		//	A common parameter
		double n = (x-(*Global_CV)[0])/(*Global_CV_err)[0];


		X_1_data.push_back( n*n );
		X_1_data2.push_back( n*n );

		sqrt_1_minus_Data_X.push_back( sqrt( X_1_data.back() ) );
		sqrt_1_minus_Theory_X.push_back( sqrt( X_1_data.back() ) );
		sqrt_1_minus_FC_X.push_back( sqrt( X_1_data.back() ) );

		if( x < (*Global_CV)[0] ) sqrt_1_minus_Data_X.back() = -sqrt_1_minus_Data_X.back();
		if( x < (*Global_CV)[0] ) sqrt_1_minus_Theory_X.back() = -sqrt_1_minus_Theory_X.back();
		if( x < (*Global_CV)[0] ) sqrt_1_minus_FC_X.back() = -sqrt_1_minus_FC_X.back();

		if( x < (*Global_CV)[0] ) X_1_data.back() = -X_1_data.back();
		if( x < (*Global_CV)[0] ) X_1_data2.back() = -X_1_data2.back();

		//	Y_1_data contains log(1-FC_CL)
		Y_1_data.push_back( 1.-ALL_CL_FROM_FC[i].second );
		//	Y_1_data2 contains log(1-Theory_CL)
		Y_1_data2.push_back( 1.-x_err );
		//	real_1_data cotains log( 1-erf(data_DLL) )
		real_1_data.push_back( 1.-erf( sqrt(DATA_DLL[i]) ) );

		sqrt_1_minus_Data_Y.push_back( real_1_data.back() );
		sqrt_1_minus_Theory_Y.push_back( Y_1_data2.back() );
		sqrt_1_minus_FC_Y.push_back( Y_1_data.back() );

		if( (i+1) != ALL_CL_FROM_FC.size() )
		{
			if( ( x < (*Global_CV)[0] && ALL_CL_FROM_FC[i+1].first[0] > (*Global_CV)[0] ) || ( x > (*Global_CV)[0] && ALL_CL_FROM_FC[i+1].first[0] < (*Global_CV)[0] ) )
			{
				X_data2.push_back( (*Global_CV)[0] );
				X_1_data2.push_back( 0. );
				Y_data2.push_back( 0. );
				Y_1_data2.push_back( 1. );
				real_1_data.push_back( 1. );
				Y_data3.push_back( 0. );
				X_data3.push_back( (*Global_CV)[0] );

				sqrt_1_minus_Data_X.push_back( 0. );
				sqrt_1_minus_Data_Y.push_back( 1. );
				sqrt_1_minus_Theory_X.push_back( 0. );
				sqrt_1_minus_Theory_Y.push_back( 1. );
			}
		}
	}



	//	Make CL plot

	TGraph* final_graph = new TGraph( (int)Y_data.size(), &(X_data[0]), &(Y_data[0]) );
	final_graph->SetTitle("");final_graph->SetName("final");
	final_graph->SetLineColor(2);
	final_graph->SetMarkerColor(2);

	TGraph* final_cl = new TGraph( (int)Y_data2.size(), &(X_data2[0]), &(Y_data2[0]) );
	final_cl->SetTitle( "" );final_graph->SetName("C_L_");
	final_cl->SetLineColor(4);
	final_cl->SetMarkerColor(4);

	TGraph* data1 = new TGraph( (int)Y_data3.size(), &(X_data3[0]), &(Y_data3[0]) );
	data1->SetTitle("");data1->SetName("data1");
	data1->SetLineColor(3);
	data1->SetMarkerColor(3);

	TMultiGraph* mult = new TMultiGraph( "multi_f", "" );
	mult->Add( final_graph );
	mult->Add( final_cl );
	mult->Add( data1 );

	TCanvas* cf = EdStyle::RapidFitCanvas( "cf", "" );
	cf->SetTitle("");
	mult->Draw("APC");
	cf->Update();

	TPaveText* text_1 = Histogram_Processing::addLHCbLabel("1.03 fb^{-1}",true);
	text_1->SetFillStyle(0);
	text_1->Draw("SAME");

	double max = mult->GetXaxis()->GetXmax();
	double min = mult->GetXaxis()->GetXmin();
	double range = fabs(max-min);
	if( max<min ) min=max;

	double min_line_1=min+range*0.25;
	double min_line_2=min+range*0.2;
	double min_line_3=min+range*0.1;

	if( fabs( min - (*Global_CV)[0] ) < 1E-1 )
	{
		min_line_1 = min;
		min_line_2 = min;
		min_line_3 = min;
	}

	TLine* _68_CL = new TLine( min_line_1, 0.68, min+0.75*range, 0.68 );
	_68_CL->SetLineColor( 6 );
	TLine* _90_CL = new TLine( min_line_1, 0.90, min+0.75*range, 0.90 );
	_90_CL->SetLineColor( 7 );
	TLine* _95_CL = new TLine( min_line_2, 0.95, min+0.8*range, 0.95 );
	_95_CL->SetLineColor( 8 );
	TLine* _99_CL = new TLine( min_line_3, 0.99, min+0.9*range, 0.99 );
	_99_CL->SetLineColor( 9 );

	_68_CL->Draw();
	_90_CL->Draw();
	_95_CL->Draw();
	_99_CL->Draw();
	cf->Update();

	TLegend *leg = new TLegend(0.77,0.9,0.97,0.6);
	//leg->SetFillStyle(1001);
	leg->SetFillStyle( 0 );
	leg->SetBorderSize(0);

	leg->AddEntry(final_graph,"FelmanCousins CL","lp");
	leg->AddEntry(data1, "CL from Data", "lp" );
	leg->AddEntry(final_cl,"Theoretical CL from Data", "lp" );
	leg->AddEntry(_68_CL, "68\% CL", "l" );
	leg->AddEntry(_90_CL, "90\% CL", "l" );
	leg->AddEntry(_95_CL, "95\% CL", "l" );
	leg->AddEntry(_99_CL, "99\% CL", "l" );

	leg->Draw();

	mult->GetYaxis()->SetTitle( "C.L." );
	mult->GetXaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter_name[0]) + " " + EdStyle::GetParamRootUnit(controlled_parameter_name[0]) );

	Histogram_Processing::Silent_Print( cf, "Final.pdf" );
	Histogram_Processing::Silent_Print( cf, "Final.png" );
	Histogram_Processing::Silent_Print( cf, "Final.C" );

	mult->GetYaxis()->SetRangeUser( 0.5, 1. );
	cf->Update();

	Histogram_Processing::Silent_Print( cf, "Final_zoom.pdf" );
	Histogram_Processing::Silent_Print( cf, "Final_zoom.png" );
	Histogram_Processing::Silent_Print( cf, "Final_zoom.C" );

	TCanvas* cfz = EdStyle::RapidFitCanvas( "cfz", "" );

	TMultiGraph* mult_z = new TMultiGraph( "cfz", "" );
	mult_z->Add( final_graph );
	mult_z->Add( final_cl );
	mult_z->Add( data1 );
	cfz->SetLogy();
	mult_z->Draw("APC");
	mult_z->GetYaxis()->SetTitle( "C.L." );
	mult_z->GetXaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter_name[0]) + " " + EdStyle::GetParamRootUnit(controlled_parameter_name[0]) );
	Histogram_Processing::Silent_Print( cfz, "logFinal.pdf" );
	Histogram_Processing::Silent_Print( cfz, "logFinal.png" );
	Histogram_Processing::Silent_Print( cfz, "logFinal.C" );



	//	Make 1-CL Plot

	TGraph* final_graph2 = new TGraph( (int)Y_1_data.size(), &(X_1_data[0]), &(Y_1_data[0]) );
	final_graph2->SetTitle("");
	final_graph2->SetLineColor(2);
	final_graph2->SetMarkerColor(2);

	TGraph* final_cl2 = new TGraph( (int)Y_1_data2.size(), &(X_1_data2[0]), &(Y_1_data2[0]) );
	final_cl2->SetTitle( "" );
	final_cl2->SetLineColor(4);
	final_cl2->SetMarkerColor(4);

	TGraph* real_cl2 = new TGraph( (int)real_1_data.size(), &(X_1_data2[0]), &(real_1_data[0]) );
	real_cl2->SetTitle( "" );
	real_cl2->SetLineColor(3);
	real_cl2->SetMarkerColor(3);

	TMultiGraph* mult2 = new TMultiGraph( "multi_f2", "" );
	mult2->Add( final_graph2 );
	mult2->Add( final_cl2 );
	mult2->Add( real_cl2 );

	TCanvas* cf2 = EdStyle::RapidFitCanvas( "cf2", "" );
	cf2->SetTitle( "" );
	cf2->SetLogy();
	mult2->Draw("APC");
	cf2->Update();


	TPaveText* text_2 = Histogram_Processing::addLHCbLabel("1.03 fb^-1",true);
	text_2->SetFillStyle(0);
	text_2->Draw("SAME");
	cf2->Update();

	min = mult2->GetXaxis()->GetXmin();
	max = mult2->GetXaxis()->GetXmax();
	range = fabs(max-min);
	if( max<min ) min=max;

	min_line_1=min+range*0.25;
	min_line_2=min+range*0.2;
	min_line_3=min+range*0.1;

	if( fabs( min/range ) < .1  )
	{
		min_line_1 = min;
		min_line_2 = min;
		min_line_3 = min;
	}

	TLine* _1_68_CL = new TLine( min_line_1, 0.32, min+range*0.75, 0.32 );
	_1_68_CL->SetLineColor( 6 );
	TLine* _1_90_CL = new TLine( min_line_1, 0.10, min+range*0.75, 0.10 );
	_1_90_CL->SetLineColor( 7 );
	TLine* _1_95_CL = new TLine( min_line_2, 0.05, min+range*0.8, 0.05 );
	_1_95_CL->SetLineColor( 8 );
	TLine* _1_99_CL = new TLine( min_line_3, 0.01, min+0.9*range, 0.01 );
	_1_99_CL->SetLineColor( 9 );

	_1_68_CL->Draw();
	_1_90_CL->Draw();
	_1_95_CL->Draw();
	_1_99_CL->Draw();
	cf2->SetLogy();
	cf2->Update();

	TLegend *leg2 = new TLegend(0.77,0.9,0.97,0.6);
	//leg2->SetFillStyle(1001);
	leg2->SetFillStyle( 0 );
	leg2->SetBorderSize(0);

	leg2->AddEntry(final_graph2,"FelmanCousins CL","lp" );
	leg2->AddEntry(real_cl2, "CL from Data", "lp" );
	leg2->AddEntry(final_cl2,"Theoretical CL from Data", "lp" );
	leg2->AddEntry(_1_68_CL, "68\% CL", "l" );
	leg2->AddEntry(_1_90_CL, "90\% CL", "l" );
	leg2->AddEntry(_1_95_CL, "95\% CL", "l" );
	leg2->AddEntry(_1_99_CL, "99\% CL", "l" );
	leg2->Draw();

	//mult2->GetYaxis()->SetRangeUser( 0.001, 1. );

	mult2->GetYaxis()->SetTitle( "1 - log( CL )" );
	mult2->GetXaxis()->SetTitle( "2 #Delta LL" );

	cf2->Update();

	Histogram_Processing::Silent_Print( cf2, "1-Final.pdf" );
	Histogram_Processing::Silent_Print( cf2, "1-Final.png" );
	Histogram_Processing::Silent_Print( cf2, "1-Final.C" );



	//	make 1-sqrt(CL) Plot

	TGraph* FC_graph = new TGraph( sqrt_1_minus_FC_X.size(), &(sqrt_1_minus_FC_X[0]), &(sqrt_1_minus_FC_Y[0]) );
	FC_graph->SetTitle("");
	FC_graph->SetLineColor(2);
	FC_graph->SetMarkerColor(2);

	TGraph* Data_Graph = new TGraph( sqrt_1_minus_Theory_X.size(), &(sqrt_1_minus_Theory_X[0]), &(sqrt_1_minus_Theory_Y[0]) );
	Data_Graph->SetTitle( "" );
	Data_Graph->SetLineColor(4);
	Data_Graph->SetMarkerColor(4);

	TGraph* Theory_Graph = new TGraph( sqrt_1_minus_Data_X.size(), &(sqrt_1_minus_Data_X[0]), &(sqrt_1_minus_Data_Y[0]) );
	Theory_Graph->SetTitle( "" );
	Theory_Graph->SetLineColor(3);
	Theory_Graph->SetMarkerColor(3);

	TMultiGraph* mult3 = new TMultiGraph( "multi_f3", "" );
	mult3->Add( FC_graph );
	mult3->Add( Theory_Graph );
	mult3->Add( Data_Graph );

	TCanvas* cf3 = EdStyle::RapidFitCanvas( "cf3", "" );
	cf3->SetTitle( "" );
	cf3->SetLogy();
	mult3->Draw("APC");
	cf3->Update();


	TPaveText* text_3 = Histogram_Processing::addLHCbLabel("1.03 fb^-1",true);
	text_3->SetFillStyle(0);
	text_3->Draw("SAME");
	cf3->Update();

	min = mult3->GetXaxis()->GetXmin();
	max = mult3->GetXaxis()->GetXmax();
	range = fabs(max-min);
	if( max<min ) min=max;

	min_line_1=min+range*0.25;
	min_line_2=min+range*0.2;
	min_line_3=min+range*0.1;

	if( fabs( min/range ) < .1 )
	{
		min_line_1 = min;
		min_line_2 = min;
		min_line_3 = min;
	}

	TLine* _1_68_CL_3 = new TLine( min_line_1, 0.32, min+range*0.75, 0.32 );
	_1_68_CL_3->SetLineColor( 6 );
	TLine* _1_90_CL_3 = new TLine( min_line_1, 0.10, min+range*0.75, 0.10 );
	_1_90_CL_3->SetLineColor( 7 );
	TLine* _1_95_CL_3 = new TLine( min_line_2, 0.05, min+range*0.8, 0.05 );
	_1_95_CL_3->SetLineColor( 8 );
	TLine* _1_99_CL_3 = new TLine( min_line_3, 0.01, min+0.9*range, 0.01 );
	_1_99_CL_3->SetLineColor( 9 );

	_1_68_CL_3->Draw();
	_1_90_CL_3->Draw();
	_1_95_CL_3->Draw();
	_1_99_CL_3->Draw();
	cf3->SetLogy();
	cf3->Update();


	TLegend *leg3 = new TLegend(0.77,0.9,0.97,0.6);
	//leg3->SetFillStyle(1001);
	leg3->SetFillStyle( 0 );
	leg3->SetBorderSize(0);

	leg3->AddEntry(FC_graph,"FelmanCousins CL","lp" );
	leg3->AddEntry(Theory_Graph, "CL from Data", "lp" );
	leg3->AddEntry(Data_Graph,"Theoretical CL from Data", "lp" );
	leg3->AddEntry(_1_68_CL_3, "68\% CL", "l" );
	leg3->AddEntry(_1_90_CL_3, "90\% CL", "l" );
	leg3->AddEntry(_1_95_CL_3, "95\% CL", "l" );
	leg3->AddEntry(_1_99_CL_3, "99\% CL", "l" );

	leg3->Draw();

	mult3->GetYaxis()->SetRangeUser( 0.001, 1. );

	mult3->GetYaxis()->SetTitle( "1 - log( CL )" );
	mult3->GetXaxis()->SetTitle( "#sqrt{2 #Delta LL}" );

	cf3->Update();

	Histogram_Processing::Silent_Print( cf3, "1-Final_sqrt.pdf" );
	Histogram_Processing::Silent_Print( cf3, "1-Final_sqrt.png" );
	Histogram_Processing::Silent_Print( cf3, "1-Final_sqrt.C" );

}

vector<OutputPlots*> FeldmanCousinsAnalysis::Plot_2D( vector<pair<vector<double>,double> > ALL_CL_FROM_FC, vector<pair<vector<double>,double> > DATA_DLL, vector<double> Global_CV, vector<double> Global_CV_err, TRandom3* rand )
{
	(void) ALL_CL_FROM_FC;
	(void) DATA_DLL;
	(void) Global_CV;
	(void) Global_CV_err;


	vector<OutputPlots*> output;
	return output;
}

vector<OutputPlots*> FeldmanCousinsAnalysis::Plot_1D( vector<pair<double, double> > ALL_CL_FROM_FC, vector<double> DATA_DLL, double Global_CV, double Global_CV_err, TRandom3* rand )
{
	vector<double> FC_X, Theory_X, Data_X, _1_minus_FC_X, _1_minus_Theory_X, _1_minus_Data_X;
	vector<double> FC_Y, Theory_Y, Data_Y, _1_minus_FC_Y, _1_minus_Theory_Y, _1_minus_Data_Y;

	for( unsigned int i=0; i< ALL_CL_FROM_FC.size(); ++i )
	{

		cout << ALL_CL_FROM_FC[i].first << "\t" << ALL_CL_FROM_FC[i].second << "\t" << DATA_DLL[i] << endl;

		double x =  ALL_CL_FROM_FC[i].first;
		double x_err = erf( fabs((x-(Global_CV/Global_CV_err))/sqrt(2.) ) );

		//      X_data and Y_data contain phi_s and FC-CL respectivley
		//      X_data2 and Y_data2 contain the theorectical CL based on the CV fit result
		FC_X.push_back( x );
		Theory_X.push_back( x );
		Data_X.push_back( x );

		FC_Y.push_back( ALL_CL_FROM_FC[i].second );
		Theory_Y.push_back( x_err );
		Data_Y.push_back( erf( sqrt(DATA_DLL[i]) ) );

		//      A common parameter
		double n = (x-Global_CV)/Global_CV_err;

		_1_minus_FC_X.push_back( n*n );
		_1_minus_Theory_X.push_back( n*n );
		if( x < Global_CV ) _1_minus_FC_X.back() = -_1_minus_FC_X.back();
		if( x < Global_CV ) _1_minus_Theory_X.back() = -_1_minus_Theory_X.back();

		//      Y_1_data contains log(1-FC_CL)
		_1_minus_FC_Y.push_back( 1.-ALL_CL_FROM_FC[i].second );
		//      Y_1_data2 contains log(1-Theory_CL)
		_1_minus_Theory_Y.push_back( 1.-x_err );
		//      real_1_data cotains log( 1-erf(data_DLL) )
		_1_minus_Data_Y.push_back( 1.-erf( sqrt(DATA_DLL[i]) ) );

		if( (i+1) != ALL_CL_FROM_FC.size() )
		{
			if( x < Global_CV && ALL_CL_FROM_FC[i+1].first > Global_CV )
			{
				Theory_X.push_back( Global_CV );
				Theory_Y.push_back( 0. );
				_1_minus_Theory_X.push_back( 0. );
				_1_minus_Theory_Y.push_back( 1. );

				Data_X.push_back( Global_CV );
				Data_Y.push_back( 0. );
				_1_minus_Data_X.push_back( 0. );
				_1_minus_Data_Y.push_back( 1. );
			}
		}
	}

	vector<OutputPlots*> returnable_plots;
	OutputPlots* plot1 = Plot_1D_Vectors( make_pair( Theory_X, Theory_Y ), make_pair( Data_X, Data_Y ), make_pair( FC_X, FC_Y ), rand );
	OutputPlots* plot2 = Plot_1D_Vectors( make_pair( _1_minus_Theory_X, _1_minus_Theory_Y), make_pair( _1_minus_Data_X, _1_minus_Data_Y ), make_pair( _1_minus_FC_X, _1_minus_FC_Y ), rand );

	returnable_plots.push_back( plot1 );
	returnable_plots.push_back( plot2 );
	return returnable_plots;
}

OutputPlots* FeldmanCousinsAnalysis::Plot_1D_Vectors( pair<vector<double>,vector<double> > theory, pair<vector<double>,vector<double> > data, pair<vector<double>,vector<double> > FC, TRandom3* rand )
{
	TGraph* FC_graph = new TGraph( (int)FC.first.size(), &(FC.first[0]), &(FC.second[0]) );
	FC_graph->SetTitle("");FC_graph->SetName("FC");
	FC_graph->SetLineColor(2);
	FC_graph->SetMarkerColor(2);

	TGraph* theory_graph = new TGraph( (int)theory.first.size(), &(theory.first[0]), &(theory.second[0]) );
	theory_graph->SetTitle( "" );theory_graph->SetName("C_L_");
	theory_graph->SetLineColor(4);
	theory_graph->SetMarkerColor(4);

	TGraph* data_graph = new TGraph( (int)data.first.size(), &(data.first[0]), &(data.second[0]) );
	data_graph->SetTitle("");data_graph->SetName("data1");
	data_graph->SetLineColor(3);
	data_graph->SetMarkerColor(3);


	TMultiGraph* mult = new TMultiGraph( "multi_f", "" );
	mult->Add( FC_graph );
	mult->Add( theory_graph );
	mult->Add( data_graph );


	TString rand_num; rand_num+=rand->Rndm();
	string rand_string( rand_num.Data() );
	replace( rand_string.begin(), rand_string.end(), '.', '_' );
	rand_num = rand_string.c_str();
	TCanvas* cf = EdStyle::RapidFitCanvas( "canv_"+rand_num, "" );
	mult->Draw("APC");
	cf->Update();

	TPaveText* text_1 = Histogram_Processing::addLHCbLabel("1.03 fb^{-1}",true);
	text_1->SetFillStyle(1001);
	text_1->Draw("SAME");

	double max = mult->GetXaxis()->GetXmax();
	double min = mult->GetXaxis()->GetXmin();
	double range = fabs(max-min);
	if( max<min ) min=max;

	TLine* _68_CL = new TLine( min+0.25*range, 0.68, min+0.75*range, 0.68 );
	_68_CL->SetLineColor( 6 );
	TLine* _90_CL = new TLine( min+0.25*range, 0.90, min+0.75*range, 0.90 );
	_90_CL->SetLineColor( 7 );
	TLine* _95_CL = new TLine( min+0.25*range, 0.95, min+0.75*range, 0.95 );
	_95_CL->SetLineColor( 8 );
	TLine* _99_CL = new TLine( min+0.10*range, 0.99, min+0.90*range, 0.99 );
	_99_CL->SetLineColor( 9 );

	_68_CL->Draw();
	_90_CL->Draw();
	_95_CL->Draw();
	_99_CL->Draw();
	cf->Update();

	TLegend *leg = new TLegend( 0.75, 0.8, 0.95, 0.6 );
	leg->SetFillStyle( 0 );
	leg->SetBorderSize( 0 );

	leg->AddEntry( FC_graph, "FelmanCousins CL", "lp" );
	leg->AddEntry( data_graph, "CL from Data", "lp" );
	leg->AddEntry( theory_graph, "Theoretical CL from Data", "lp" );
	leg->AddEntry( _68_CL, "68\% CL", "l" );
	leg->AddEntry( _90_CL, "90\% CL", "l" );
	leg->AddEntry( _95_CL, "95\% CL", "l" );
	leg->AddEntry( _99_CL, "99\% CL", "l" );

	leg->Draw();

	cf->Update();

	vector<TObject*> other_objects;
	other_objects.push_back( _68_CL );
	other_objects.push_back( _90_CL );
	other_objects.push_back( _95_CL );
	other_objects.push_back( _99_CL );

	OutputPlots* returnable_plot = new OutputPlots( "FC-Plot", cf, (TObject*) mult, leg, text_1, "APL" );

	return returnable_plot;
}


/*
   {

   vector<double> X_data, X_data2, X_data3, X_1_data, X_1_data2, real_1_data, Y_data, Y_data2, Y_data3, Y_1_data, Y_1_data2;

   for( unsigned int i=0; i< ALL_CL_FROM_FC.size(); ++i )
   {
   cout << ALL_CL_FROM_FC[i].first[0] << "\t" << ALL_CL_FROM_FC[i].second << "\t" << DATA_DLL[i] << endl;

   double x =  ALL_CL_FROM_FC[i].first[0];
   double x_err = erf( fabs((x-(*Global_CV)[0])/(*Global_CV_err)[0])/sqrt(2.) );

//	X_data and Y_data contain phi_s and FC-CL respectivley
//	X_data2 and Y_data2 contain the theorectical CL based on the CV fit result
X_data.push_back( x );
X_data2.push_back( x );
X_data3.push_back( x );
Y_data.push_back( ALL_CL_FROM_FC[i].second );
Y_data2.push_back( x_err );
Y_data3.push_back( erf( sqrt(DATA_DLL[i]) ) );

//	A common parameter
double n = (x-(*Global_CV)[0])/(*Global_CV_err)[0];


X_1_data.push_back( n*n );
X_1_data2.push_back( n*n );		
if( x < (*Global_CV)[0] ) X_1_data.back() = -X_1_data.back();
if( x < (*Global_CV)[0] ) X_1_data2.back() = -X_1_data2.back();

//	Y_1_data contains log(1-FC_CL)
Y_1_data.push_back( 1.-ALL_CL_FROM_FC[i].second );
//	Y_1_data2 contains log(1-Theory_CL)
Y_1_data2.push_back( 1.-x_err );
//	real_1_data cotains log( 1-erf(data_DLL) )
real_1_data.push_back( 1.-erf( sqrt(DATA_DLL[i]) ) );

if( (i+1) != ALL_CL_FROM_FC.size() )
{
if( x < (*Global_CV)[0] && ALL_CL_FROM_FC[i+1].first[0] > (*Global_CV)[0] )
{
X_data2.push_back( (*Global_CV)[0] );
X_1_data2.push_back( 0. );
Y_data2.push_back( 0. );
Y_1_data2.push_back( 1. );
real_1_data.push_back( 1. );
Y_data3.push_back( 0. );
X_data3.push_back( (*Global_CV)[0] );
}
}
}

TGraph* final_graph = new TGraph( (int)Y_data.size(), &(X_data[0]), &(Y_data[0]) );
final_graph->SetTitle("");final_graph->SetName("final");
final_graph->SetLineColor(2);
final_graph->SetMarkerColor(2);

TGraph* final_cl = new TGraph( (int)Y_data2.size(), &(X_data2[0]), &(Y_data2[0]) );
final_cl->SetTitle( "" );final_graph->SetName("C.L.");
final_cl->SetLineColor(4);
final_cl->SetMarkerColor(4);

TGraph* data1 = new TGraph( (int)Y_data3.size(), &(X_data3[0]), &(Y_data3[0]) );
data1->SetTitle("");data1->SetName("data1");
data1->SetLineColor(3);
data1->SetMarkerColor(3);

TMultiGraph* mult = new TMultiGraph( "multi_f", "" );
mult->Add( final_graph );
mult->Add( final_cl );
mult->Add( data1 );

TCanvas* cf = EdStyle::RapidFitCanvas( "cf", "" );
cf->SetTitle("");
mult->Draw("APC");
cf->Update();

TPaveText* text_1 = Histogram_Processing::addLHCbLabel("1.03 fb^{-1}",true);
text_1->SetFillStyle(1001);
text_1->Draw("SAME");

double max = mult->GetXaxis()->GetXmax();
double min = mult->GetXaxis()->GetXmin();
double range = fabs(max-min);
if( max<min ) min=max;

TLine* _68_CL = new TLine( min+0.25*range, 0.68, min+0.75*range, 0.68 );
_68_CL->SetLineColor( 6 );
TLine* _90_CL = new TLine( min+0.25*range, 0.90, min+0.75*range, 0.90 );
_90_CL->SetLineColor( 7 );
TLine* _95_CL = new TLine( min+0.25*range, 0.95, min+0.75*range, 0.95 );
_95_CL->SetLineColor( 8 );
TLine* _99_CL = new TLine( min+0.1*range, 0.99, min+0.9*range, 0.99 );
_99_CL->SetLineColor( 9 );

_68_CL->Draw();
_90_CL->Draw();
_95_CL->Draw();
_99_CL->Draw();
cf->Update();

TLegend *leg = new TLegend(0.75,0.8,0.95,0.6);
//leg->SetFillStyle(1001);
leg->SetFillStyle( 0 );
leg->SetBorderSize(0);

leg->AddEntry(final_graph,"FelmanCousins CL","lp");
leg->AddEntry(data1, "CL from Data", "lp" );
leg->AddEntry(final_cl,"Theoretical CL from Data", "lp" );
leg->AddEntry(_68_CL, "68\% CL", "l" );
leg->AddEntry(_90_CL, "90\% CL", "l" );
leg->AddEntry(_95_CL, "95\% CL", "l" );
leg->AddEntry(_99_CL, "99\% CL", "l" );

leg->Draw();

mult->GetYaxis()->SetTitle( "C.L." );
mult->GetXaxis()->SetTitle( EdStyle::GetParamRootName(controlled_parameter_name[0]) + " " + EdStyle::GetParamRootUnit(controlled_parameter_name[0]) );

Histogram_Processing::Silent_Print( cf, "Final.pdf" );
Histogram_Processing::Silent_Print( cf, "Final.png" );
Histogram_Processing::Silent_Print( cf, "Final.C" );


TGraph* final_graph2 = new TGraph( (int)Y_1_data.size(), &(X_1_data[0]), &(Y_1_data[0]) );
final_graph2->SetTitle("");
final_graph2->SetLineColor(2);
final_graph2->SetMarkerColor(2);

TGraph* final_cl2 = new TGraph( (int)Y_1_data2.size(), &(X_1_data2[0]), &(Y_1_data2[0]) );
final_cl2->SetTitle( "" );
final_cl2->SetLineColor(4);
final_cl2->SetMarkerColor(4);

TGraph* real_cl2 = new TGraph( (int)real_1_data.size(), &(X_1_data2[0]), &(real_1_data[0]) );
real_cl2->SetTitle( "" );
real_cl2->SetLineColor(3);
real_cl2->SetMarkerColor(3);

TMultiGraph* mult2 = new TMultiGraph( "multi_f2", "" );
mult2->Add( final_graph2 );
mult2->Add( final_cl2 );
mult2->Add( real_cl2 );

TCanvas* cf2 = EdStyle::RapidFitCanvas( "cf2", "" );
cf2->SetTitle( "" );
cf2->SetLogy();
mult2->Draw("APC");
cf2->Update();


TPaveText* text_2 = Histogram_Processing::addLHCbLabel("1.03 fb^-1",true);
text_2->SetFillStyle(1001);
text_2->Draw("SAME");
cf2->Update();

min = mult2->GetXaxis()->GetXmin();
max = mult2->GetXaxis()->GetXmax();
range = fabs(max-min);
if( max<min ) min=max;

TLine* _1_68_CL = new TLine( min+range*0.25, 0.32, min+range*0.75, 0.32 );
_1_68_CL->SetLineColor( 6 );
TLine* _1_90_CL = new TLine( min+range*0.25, 0.10, min+range*0.75, 0.10 );
_1_90_CL->SetLineColor( 7 );
TLine* _1_95_CL = new TLine( min+range*0.25, 0.05, min+range*0.75, 0.05 );
_1_95_CL->SetLineColor( 8 );
TLine* _1_99_CL = new TLine( min+0.1*range, 0.01, min+0.9*range, 0.01 );
_1_99_CL->SetLineColor( 9 );

_1_68_CL->Draw();
_1_90_CL->Draw();
_1_95_CL->Draw();
_1_99_CL->Draw();
cf2->SetLogy();
cf2->Update();

TLegend *leg2 = new TLegend(0.75,0.8,0.95,0.6);
//leg2->SetFillStyle(1001);
leg2->SetFillStyle( 0 );
leg2->SetBorderSize(0);

leg2->AddEntry(final_graph2,"FelmanCousins CL","lp" );
leg2->AddEntry(real_cl2, "CL from Data", "lp" );
leg2->AddEntry(final_cl2,"Theoretical CL from Data", "lp" );
leg2->AddEntry(_1_68_CL, "68\% CL", "l" );
leg2->AddEntry(_1_90_CL, "90\% CL", "l" );
leg2->AddEntry(_1_95_CL, "95\% CL", "l" );
leg2->AddEntry(_1_99_CL, "99\% CL", "l" );

leg2->Draw();

mult2->GetYaxis()->SetRangeUser( 0.001, 1. );

mult2->GetYaxis()->SetTitle( "1 - log( CL )" );
mult2->GetXaxis()->SetTitle( "2 #Delta LL" );

cf2->Update();

Histogram_Processing::Silent_Print( cf2, "1-Final.pdf" );

}
*/

