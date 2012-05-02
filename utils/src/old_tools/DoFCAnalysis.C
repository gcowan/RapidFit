
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

#include "DoFCAnalysis.h"
#include "NTuple_Processing.h"
#include "Histo_Processing.h"
#include "Template_Functions.h"
#include "RapidFitParameters.h"

#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace::std;

void Print_Info_At_Grid_ij( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, vector<vector<Double_t> > this_Grid_Coordinate )
{
	cout << "AT: ";
	vector<string>::iterator param_i = controlled_parameter_name.begin();
	vector<Double_t>::iterator grid_ij = coordinate.begin();
	for( ; param_i != controlled_parameter_name.end(); ++param_i )
	{
		cout << *param_i << ": " << *grid_ij << "\t";
	}
	cout << "\t\t" << "THERE ARE: " << this_Grid_Coordinate.size() << " TOYS!" << endl;
}

void Print_Toy_Info( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, int Fixed, int Free )
{
	cout << "AT: ";
	vector<string>::iterator param_i = controlled_parameter_name.begin();
	vector<Double_t>::iterator grid_ij = coordinate.begin();
	for( ; param_i != controlled_parameter_name.end(); ++param_i )
	{
		cout << *param_i << ": " << *grid_ij << "\t";
	}
	cout << "\t\t" << "THERE ARE:\t" << Fixed << " Fixed Toys and " << Free << " Free Toys." << endl;
}

TString Flat_Coord_String( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate )
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

void Plot_Nuisance_Parameters( TTree* input_tree, TString& free_toys, TString& fixed_toys, vector<TString>& Free_Params, vector<string> controlled_parameter_name, vector<Double_t> coordinate, TRandom3* rand, double DLL_Data )
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
			//vector<vector<Double_t> > Free_this_Grid_Coordinate = Plotter_Data( input_tree, *Param_i+*var_i, free_toys, rand );
			//vector<vector<Double_t> > Fixed_this_Grid_Coordinate = Plotter_Data( input_tree, *Param_i+*var_i, fixed_toys, rand );

			vector<Double_t>* param_data_free = Buffer_Branch( input_tree, *Param_i+*var_i, free_toys ); //Free_this_Grid_Coordinate[0];
			vector<Double_t>* nll_free = Buffer_Branch( input_tree, "NLL", free_toys );
			vector<Double_t>* nll_free_As = Buffer_Branch( input_tree, "As_sq_value", free_toys );
			vector<Double_t>* nll_free_dpa = Buffer_Branch( input_tree, "delta_para_value", free_toys );
			vector<Double_t>* nll_free_dpe = Buffer_Branch( input_tree, "delta_perp_value", free_toys );
			vector<Double_t>* nll_free_ds = Buffer_Branch( input_tree, "delta_s_value", free_toys );
			vector<Double_t>* param_data_fixed = Buffer_Branch( input_tree, *Param_i+*var_i, fixed_toys ); //Fixed_this_Grid_Coordinate[0];
			vector<Double_t>* nll_fixed_As = Buffer_Branch( input_tree, "As_sq_value", fixed_toys );
			vector<Double_t>* nll_fixed_dpa = Buffer_Branch( input_tree, "delta_para_value", fixed_toys );
			vector<Double_t>* nll_fixed_dpe = Buffer_Branch( input_tree, "delta_perp_value", fixed_toys );
			vector<Double_t>* nll_fixed = Buffer_Branch( input_tree, "NLL", fixed_toys );
			vector<Double_t>* nll_fixed_ds = Buffer_Branch( input_tree, "delta_s_value", fixed_toys );
			for( unsigned int i=0; i != param_data_free->size(); ++i )
			{
				if( ((*nll_fixed)[i]-(*nll_free)[i]) < 0 )
				{
					cout << "BAD:\t" << *Param_i+*var_i << "\tFREE: " << (*param_data_free)[i] << "\tFIXED: " << (*param_data_fixed)[i] << "\t\tNLL: " << (*nll_free)[i] << "\t"  << (*nll_fixed)[i] << "\tDLL: " << ((*nll_fixed)[i]-(*nll_free)[i]) << endl;
				}
				if( ((*nll_fixed)[i]-(*nll_free)[i]) > 2*DLL_Data && ((*nll_free_As)[i] > 0.0001 && (*nll_fixed_As)[i] > 0.0001 ) )
				{
					if( fabs( (*nll_free_dpa)[i] ) < 6. && fabs((*nll_fixed_dpa)[i]) < 6. )
					{
					if( fabs( (*nll_free_dpe)[i] ) < 6. && fabs((*nll_fixed_dpe)[i]) < 6. )
					{
					if( fabs( (*nll_free_ds)[i] ) < 6. && fabs((*nll_fixed_ds)[i]) < 6. )
					{
					if( (*nll_free_dpa)[i] > 0. && (*nll_fixed_dpa)[i] > 0. )
					{
					if( (*nll_free_dpe)[i] > 0. && (*nll_fixed_dpe)[i] > 0. )
					{
					if( (*nll_free_ds)[i] > 0. && (*nll_fixed_ds)[i] > 0. )
					{
					cout << "BAD2:\t" << *Param_i+*var_i << "\tFREE: " << (*param_data_free)[i] << "\tFIXED: " << (*param_data_fixed)[i] << "\t\tNLL: " << (*nll_free)[i] << "\t"  << (*nll_fixed)[i] << "\tDLL: " << ((*nll_fixed)[i]-(*nll_free)[i]) << endl;
				}	}}} }}}
			}

			TCanvas* free_c1 = new TCanvas( "TCanv_free_"+rand_str, "TCanv_free_"+rand_str, 1680, 1050 );
			TH1* free_data_th1 = Get_TH1( *param_data_free, rand );//, get_optimal_histo_bins(param_data_free) );
			free_data_th1->Draw();
			free_c1->Update();
			free_c1->Print( *Param_i+"_"+*var_i+"_"+grid_name+"_free.pdf" );

			TCanvas* fixed_c1 = new TCanvas( "TCanv_fixed_"+rand_str, "TCanv_fixed_"+rand_str, 1680, 1050 );
			TH1* fixed_data_th1 = Get_TH1( *param_data_fixed, rand );//, get_optimal_histo_bins(param_data_fixed) );
			fixed_data_th1->Draw();
			fixed_c1->Update();
			fixed_c1->Print( *Param_i+"_"+*var_i+"_"+grid_name+"_fixed.pdf" );

			TGraph* free_graph = new TGraph( free_data_th1 );
			TGraph* fixed_graph = new TGraph( fixed_data_th1 );

			TMultiGraph* multi_graph = new TMultiGraph( "TMulti_"+rand_str, "TMulti_"+rand_str );

			multi_graph->Add( free_graph );
			multi_graph->Add( fixed_graph );

			TCanvas* c1 = new TCanvas( "TCanv_"+rand_str, "TCanv_"+rand_str, 1680, 1050 );

			multi_graph->Draw( "APL" );

			c1->Update();
			c1->Print( *Param_i+"_"+*var_i+"_"+grid_name+"_multi.pdf" );
		}
	}
}

void NLL_dists( vector<Double_t>& param_data_free, vector<Double_t>& param_data_fixed, vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, TRandom3* rand )
{
	TString grid_name = Flat_Coord_String( controlled_parameter_name, coordinate );
	TString rand_str; rand_str += rand->Rndm();
	TCanvas* free_c1 = new TCanvas( "TCanv_free_"+rand_str, "TCanv_free_"+rand_str, 1680, 1050 );
	TH1* free_data_th1 = Get_TH1( param_data_free, rand );//, get_optimal_histo_bins(param_data_free) );
	free_data_th1->Draw();
	free_c1->Update();
	free_c1->Print( "NLL_"+grid_name+"_free.pdf" );

	TCanvas* fixed_c1 = new TCanvas( "TCanv_fixed_"+rand_str, "TCanv_fixed_"+rand_str, 1680, 1050 );
	TH1* fixed_data_th1 = Get_TH1( param_data_fixed, rand );//, get_optimal_histo_bins(param_data_fixed) );
	fixed_data_th1->Draw();
	fixed_c1->Update();
	fixed_c1->Print( "NLL_"+grid_name+"_fixed.pdf" );

	TGraph* free_graph = new TGraph( free_data_th1 );
	TGraph* fixed_graph = new TGraph( fixed_data_th1 );

	TMultiGraph* multi_graph = new TMultiGraph( "TMulti_"+rand_str, "TMulti_"+rand_str );

	free_graph->SetLineColor(1);
	multi_graph->Add( free_graph );
	multi_graph->Add( fixed_graph );

	TCanvas* c1 = new TCanvas( "TCanv_"+rand_str, "TCanv_"+rand_str, 1680, 1050 );

	multi_graph->Draw("APL");

	c1->Update();
	c1->Print( "NLL_"+grid_name+"_multi.pdf" );
}

void Output_GridPoint( vector<string>& controlled_parameter_name, vector<Double_t>& coordinate, double LOCAL_DATA_DLL, vector<Double_t>& Toy_DLL_Dist )
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

	TCanvas* c1 = new TCanvas( "TCanvas_"+grid_name, "TCanvas_"+grid_name, 1680, 1050 );

	grid_th1->Draw("PE");
	c1->Update();           //      STUPID ROOT

	double y_max = grid_th1->GetBinContent( grid_th1->GetMaximumBin() );

	TLine* data_DLL = new TLine( (Double_t)LOCAL_DATA_DLL, 0., (Double_t)LOCAL_DATA_DLL, (Double_t)y_max );

	data_DLL->SetLineWidth(8);

	data_DLL->Draw( "SAME" );
	c1->Update();

	c1->Print("DLL_Dist"+grid_name+".pdf");
	c1->Print("DLL_Dist"+grid_name+".png");

	grid_th1->Write("",TObject::kOverwrite);
	c1->Write("",TObject::kOverwrite);

	grid_point->Close();
}

pair<vector<double>*, vector<double>* > Get_Best_CV( TTree* input_tree, vector<string>& controlled_param_name )
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

void DoFCAnalysis( TTree* input_tree, vector<string>& controlled_parameter_name, TRandom3* rand )
{
	//	Just to make sure we never lose any informatio
	input_tree->SetEstimate( input_tree->GetEntries() );

	//	DLL coordinates have the same _gen as the global CV, but all have unique fixed _value
	vector<TString> Free_Params = get_free_non_scanned_parameters( input_tree, controlled_parameter_name );

	//	This Draw String gives the returned control parameters from the cut string
	TString Grid_Draw_String = Construct_Draw_String( controlled_parameter_name );

	//	This Cut String returns only the fit results which are 
	TString Grid_Cut_String = Construct_Cut_String( input_tree, controlled_parameter_name );

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
	vector<vector<Double_t> > Grid_Coordinates = Plotter_Data( input_tree, Grid_Draw_String, Grid_Cut_String, rand );
	Grid_Coordinates = rotate( Grid_Coordinates );


	//	Output for the user
	cout << "DIMENTION: " << Grid_Coordinates[0].size() << "\t" << "COORDINATES: " << Grid_Coordinates.size() << endl;



	//	Vectors to store the output from the analysis
	vector<double> DATA_DLL;
	vector<pair<vector<double>, double> > ALL_CL_FROM_FC;

	//	Make an output Directory for the FC analysis tool
	gSystem->mkdir("FC_DLL_distributions");
	gSystem->cd("FC_DLL_distributions");

	vector<vector<Double_t> >::iterator grid_i = Grid_Coordinates.begin();
	for( ; grid_i != Grid_Coordinates.end(); ++grid_i )
	{
		//	Cut String to select only a SINGLE fit result for the LL at this point
		TString data_cut = Grid_Cut_String + "&&" + Data_At_Grid_ij( controlled_parameter_name, *grid_i ); 

		//	Slect only toy Results which were generated at this grid point
		//	By definition the _gen value is equal to the fit minima from the LL scan result at this point
		//	Could even expand on this to use the _gen of all nuisence too, but thata extreme belt&braces
		TString toys_cut = Fits_At_Grid_ij( controlled_parameter_name, *grid_i, true );

		//	Cut string which will return only free toys
		TString free_condition = Construct_Fixed_condition( controlled_parameter_name, false );

		//	Cut string which will return only fixed toys
		TString fixed_condition = Construct_Fixed_condition( controlled_parameter_name, true );

		//	Only Fixed/Free toys at this coordinate
		TString free_toys = toys_cut + "&&" + free_condition + "&&(NLL>0)";
		TString fixed_toys = toys_cut + "&&" + fixed_condition + "&&(NLL>0)";

		//	Some output for the user
		cout << "Free CUT:\t" << free_toys << endl;
		cout << "Fixed CUT:\t" << fixed_toys << endl;

		//	this data will be mangled, however we only want the single unique result so don't care
		vector<vector<Double_t> > Data_Coordinate = Plotter_Data( input_tree, "NLL", data_cut, rand );

		//	Get the NLL for the fixed/free toys in a completely UNMANGLED way
		vector<Double_t>* Fixed_Toy_NLL = Buffer_Branch( input_tree, "NLL", fixed_toys );
		vector<Double_t>* Fixed_Toy_As = Buffer_Branch( input_tree, "As_sq_value", fixed_toys );
		vector<Double_t>* Fixed_Toy_dpa = Buffer_Branch( input_tree, "delta_para_value", fixed_toys );
		vector<Double_t>* Fixed_Toy_dpe = Buffer_Branch( input_tree, "delta_perp_value", fixed_toys );
		vector<Double_t>* Fixed_Toy_dsp = Buffer_Branch( input_tree, "delta_s_pull", fixed_toys );
		vector<Double_t>* Free_Toy_NLL = Buffer_Branch( input_tree, "NLL", free_toys );
		vector<Double_t>* Free_Toy_As = Buffer_Branch( input_tree, "As_sq_value", free_toys );
		vector<Double_t>* Free_Toy_dpa = Buffer_Branch( input_tree, "delta_para_value", free_toys );
		vector<Double_t>* Free_Toy_dpe = Buffer_Branch( input_tree, "delta_perp_value", free_toys );
		vector<Double_t>* Free_Toy_dsp = Buffer_Branch( input_tree, "delta_s_pull", free_toys );

		//	Draw the NLL distributions at this point
		NLL_dists( *Free_Toy_NLL, *Fixed_Toy_NLL, controlled_parameter_name, *grid_i, rand );

		//	Print some information on the toys at this grid point
		Print_Toy_Info( controlled_parameter_name, *grid_i, (int)Fixed_Toy_NLL->size(), (int)Free_Toy_NLL->size() );

                //      Record and output the DLL from data at this point
                double LOCAL_DATA_DLL = Data_Coordinate[0][0] - GLOBAL_BEST_NLL;
                cout << LOCAL_DATA_DLL << endl;
                DATA_DLL.push_back( LOCAL_DATA_DLL );

		//	Make some plots of the distributions of all free parameters fluctuating over the fit
		//Plot_Nuisance_Parameters( input_tree, free_toys, fixed_toys, Free_Params, controlled_parameter_name, *grid_i, rand, LOCAL_DATA_DLL );




		//	Check the total number of 'good' toys which we can use at this grid point
		unsigned int toys_to_test = Free_Toy_NLL->size();
		if( Free_Toy_NLL->size() != Fixed_Toy_NLL->size() )
		{
			cerr << "\tWARNING: DIFFERENT NUMBERS OF TOYS BETWEEN FIXED/FREE" << endl;
			//	This if statement is to protect the CL from errors in the data
			exit(-5);
			//continue;
		}
		if( toys_to_test == 0 ) continue;


		//	Store the DLL from the toys in an array
		vector<double> Toy_DLL_Dist;
		for( unsigned int i=0; i != toys_to_test; ++i )
		{
			if( (*Fixed_Toy_As)[i] > 0.0001 && (*Free_Toy_As)[i] > 0.0001 )
			{
				if( fabs((*Free_Toy_dpa)[i]) < 6 && fabs((*Fixed_Toy_dpa)[i]) < 6 )
				{	
				if( fabs((*Free_Toy_dpe)[i]) < 6 && fabs((*Fixed_Toy_dpe)[i]) < 6 )
				{
				if( fabs((*Free_Toy_dsp)[i]) < 6. && fabs((*Fixed_Toy_dsp)[i]) < 6. )
				{

				if( (*Free_Toy_dpa)[i] > 0. && (*Fixed_Toy_dpa)[i] > 0. )
				{
				if( (*Free_Toy_dpe)[i] > 0. && (*Fixed_Toy_dpe)[i] > 0. )
				{
				if( (*Free_Toy_dsp)[i] > 0. && (*Fixed_Toy_dsp)[i] > 0. )
				{

				//cout << Fixed_Toy_NLL[i] << "\t" << Free_Toy_NLL[i] << "\t\t" << Fixed_Toy_NLL[i] - Free_Toy_NLL[i] << endl;
				double toy_dll = (*Fixed_Toy_NLL)[i] - (*Free_Toy_NLL)[i];
				if( toy_dll >= 0. ) Toy_DLL_Dist.push_back( toy_dll );
			}	}}}	}}}
		}
		sort( Toy_DLL_Dist.begin(), Toy_DLL_Dist.end() );

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


		//	Print some more output for this Grid Point
		cout << "DATA: " << LOCAL_DATA_DLL << "\tTOY: " << LOCAL_CL_FROM_FC << endl;
		ALL_CL_FROM_FC.push_back( make_pair( *grid_i, LOCAL_CL_FROM_FC ) );
	}



	vector<double> X_data, X_data2, Y_data, Y_data2, Y_1_data, Y_1_data2;

	for( unsigned int i=0; i< ALL_CL_FROM_FC.size(); ++i )
	{
		cout << ALL_CL_FROM_FC[i].first[0] << "\t" << ALL_CL_FROM_FC[i].second << "\t" << DATA_DLL[i] << endl;
		X_data.push_back( ALL_CL_FROM_FC[i].first[0] );
		X_data2.push_back( ALL_CL_FROM_FC[i].first[0] );
		Y_data2.push_back( erf( fabs((X_data2.back()-(*Global_CV)[0])/(*Global_CV_err)[0])/sqrt(2.) ) );
		Y_1_data2.push_back( 1.-erf( fabs((X_data2.back()-(*Global_CV)[0])/(*Global_CV_err)[0])/sqrt(2.) ) );
		if( (i+1) != ALL_CL_FROM_FC.size() )
		{
			if( ALL_CL_FROM_FC[i].first[0] < (*Global_CV)[0] && ALL_CL_FROM_FC[i+1].first[0] > (*Global_CV)[0] )
			{
				X_data2.push_back( (*Global_CV)[0] );
				Y_data2.push_back( 0. );
				Y_1_data2.push_back( 1. );
			}
		}
		Y_data.push_back( ALL_CL_FROM_FC[i].second );
		Y_1_data.push_back( 1.-ALL_CL_FROM_FC[i].second );
	}

	TGraph* final_graph = new TGraph( Y_data.size(), &(X_data[0]), &(Y_data[0]) );
	final_graph->SetTitle("final");final_graph->SetName("final");
	final_graph->SetLineColor(2);
	final_graph->SetMarkerColor(2);

	TGraph* final_cl = new TGraph( Y_data2.size(), &(X_data2[0]), &(Y_data2[0]) );
	final_cl->SetTitle( "CL" );final_graph->SetName("CL");
	final_cl->SetLineColor(4);
	final_cl->SetMarkerColor(4);

	TMultiGraph* mult = new TMultiGraph( "multi_f", "multi_f" );
	mult->Add( final_graph );
	mult->Add( final_cl );

	TCanvas* cf = new TCanvas( "cf", "cf", 1680, 1050 );
	mult->Draw("APC");
	cf->Update();
	cf->Print("Final.pdf");

	//cf->SetLogy();
	//cf->Update();

        TGraph* final_graph2 = new TGraph( Y_1_data.size(), &(X_data[0]), &(Y_1_data[0]) );
        final_graph2->SetTitle("final2");final_graph->SetName("final2");
        final_graph2->SetLineColor(2);
        final_graph2->SetMarkerColor(2);

	TGraph* final_cl2 = new TGraph( Y_1_data2.size(), &(X_data2[0]), &(Y_1_data2[0]) );
        final_cl2->SetTitle( "CL2" );final_graph->SetName("CL2");
        final_cl2->SetLineColor(4);
        final_cl2->SetMarkerColor(4);

	TMultiGraph* mult2 = new TMultiGraph( "multi_f2", "multi_f2" );
	mult2->Add( final_graph2 );
	mult2->Add( final_cl2 );

	TCanvas* cf2 = new TCanvas( "cf2", "cf2", 1680, 1050 );

	cf2->SetLogy();
	mult2->Draw("APC");
	cf2->Update();

	cf2->SetLogy();
	cf2->Update();
	cf2->Print("1-Final.pdf");
}

/*
   void Plot_Output( vector<pair<vector<double>,double> >& input, vector<string>& controlled_parameter_name )
   {
   if( vector[0].first.size() == 1 )
   Plot_Output_1D( input, controlled_parameter_name );
   else if( vector[0].first.size() == 2 )
   Plot_Output_2D( input, controlled_parameter_name );
   return
   }

   void Plot_Output_1D( vector<pair<vector<double>,double> >& input, vector<string>& controlled_parameter_name )
   {
   vector<pair<vector<double>,double> >::iterator grid_point_i = input.begin();

   TString param_name = controlled_parameter_name[0];

   vector<double> X_data;
   vector<double> Y_data;

   for( ; grid_point_i != input.end(); ++grid_point_i )
   {
   double fc_CL = grid_point_i->second;
   double param_val = grid_point_i->first[0];
   X_data.push_back( param_val );
   Y_data.push_bacl( fc_CL );
   }

   TGraph* final_graph = new TGraph( (int)input.size(), X_data, Y_data );
   final_graph->SetTitle( param_name+"_TGraph" );
   }

   void Plot_Output_2D( vector<pair<vector<double>,double> >& input, vector<string>& controlled_parameter_name )
   {

   }
 */

