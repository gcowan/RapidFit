//	This program is intended to give you an overlay of 1D LL plots formatted in the RapidFit style
//	It DOES work, however due to ROOT bugs it requires a later version to work correctly

//	ROOT Headers
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TLine.h"
#include "TColor.h"
#include "TH2.h"
//	RapidFit Headers
#include "EdStyle.h"
#include "NTuple_Processing.h"
#include "TString_Processing.h"
//	System Headers
#include <vector>
#include <iostream>
#include <cstdlib>

using namespace::std;

int main( int argc, char* argv[] )
{
	cout <<" .----------------.  .----------------.  .----------------.  .----------------.  .----------------."<<endl;
	cout <<" | .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |"<<endl;
	cout <<" | |  _______     | || |      __      | || |   ______     | || |     _____    | || |  ________    | |"<<endl;
	cout <<" | | |_   __ \\    | || |     /  \\     | || |  |_   __ \\   | || |    |_   _|   | || | |_   ___ `.  | |"<<endl;
	cout <<" | |   | |__) |   | || |    / /\\ \\    | || |    | |__) |  | || |      | |     | || |   | |   `. \\ | |"<<endl;
	cout <<" | |   |  __ /    | || |   / ____ \\   | || |    |  ___/   | || |      | |     | || |   | |    | | | |"<<endl;
	cout <<" | |  _| |  \\ \\_  | || | _/ /    \\ \\_ | || |   _| |_      | || |     _| |_    | || |  _| |___.' / | |"<<endl;
	cout <<" | | |____| |___| | || ||____|  |____|| || |  |_____|     | || |    |_____|   | || | |________.'  | |"<<endl;
	cout <<" | |              | || |              | || |              | || |              | || |              | |"<<endl;
	cout <<" | '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |"<<endl;
	cout <<"  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'"<<endl;
	cout <<"  .----------------.  .----------------."<<endl;
	cout <<" | .--------------. || .--------------. |"<<endl;
	cout <<" | |   _____      | || |   _____      | |"<<endl;
	cout <<" | |  |_   _|     | || |  |_   _|     | |"<<endl;
	cout <<" | |    | |       | || |    | |       | |"<<endl;
	cout <<" | |    | |   _   | || |    | |   _   | |"<<endl;
	cout <<" | |   _| |__/ |  | || |   _| |__/ |  | |"<<endl;
	cout <<" | |  |________|  | || |  |________|  | |"<<endl;
	cout <<" | |              | || |              | |"<<endl;
	cout <<" | '--------------' || '--------------' |"<<endl;
	cout <<"  '----------------'  '----------------'"<<endl;
	cout <<endl;
	cout << "Usage:"<<endl;
	cout << endl << argv[0] << "\t" << "Param_to_plot" << "\t" << "File1.root" << "\t" << "File2.root" << "\t...\t" << "FileN.root" << endl;
	cout <<endl;

	if( argc < 3 ) exit(-5);
	//	Do this before any Canvas is constructed
	//	ROOT is a BITCH for not allowing you to easily change internal crap
	//      Setup the Canvas and such
	EdStyle* RapidFit_Style = new EdStyle();
	RapidFit_Style->SetStyle();

	//	The name of output tree's from RapidFit assuming Flat NTuple format
	TString TREE_NAME( "RapidFitResult" );

	//	Parameter Passed from the command line
	TString Param_Of_Choice( argv[1] );

	//	switch for plotting the _error plots for each parameter of interest
	bool want_ROADMAP_ERROR=true;


	//	Input Object Holders
	vector<TString> Input_File_Names;
	vector<TFile*> Input_Files;
	vector<TTree*> Input_Tree_per_File;
	//	All args>=3 are assumed to be files due to design
	for( int i=2; i < argc ; ++i )
	{
		Input_File_Names.push_back( TString( argv[i] ) );

		Input_Files.push_back( new TFile( Input_File_Names.back(), "READ" ) );

		Input_Tree_per_File.push_back( (TTree*) gDirectory->Get( TREE_NAME ) );

		TString new_name = TREE_NAME;
		new_name+=i;

		//	ROOT does NOT like objects having the same name, keep it from throwing a hissy fit!
		Input_Tree_per_File.back()->SetName( new_name );
	}

	cout << endl <<"Opened files..." << endl<<endl;


	TString value_suffix = "_value";

	//      Get a list of all branches in allresults
	vector<TString> all_parameters = get_branch_names( Input_Tree_per_File[0] );
	//      Get a list of all branches in allresults with '_value' in their name
	vector<TString> all_parameter_values = filter_names( all_parameters, string(value_suffix.Data()) );


	//	Create a vector to hold the overlay graphs for all of the parameters in the fit
	vector<TMultiGraph*> all_params_drift_multi;
	//all_params_drift_multi.resize( all_parameters.size() );

	for( unsigned int i=0; i< all_parameters.size(); ++i )
	{
		all_params_drift_multi.push_back( new TMultiGraph() );
	}

	TMultiGraph *mg = new TMultiGraph();

	TString NLL="NLL";


	vector<double> minimum_NLL;
	vector<TString> ALL_Draw_Strings;
	vector<TTree*> temp_trees;

	cout << "Calculating DLL..." << endl<<endl;

	//	For all Input Files, derrive the string that will allow the DeltaLL to start at minima
	for( unsigned int i=0; i< Input_Tree_per_File.size(); ++i )
	{

		TString _1D_Draw_String, NLL_Draw_String;

		TTree* temp_tree = Input_Tree_per_File[i]->CopyTree( "NLL>0","fast",Input_Tree_per_File[i]->GetEntries(),0 );

		TString Name="TEMP_ROOT";
		Name+=i;
		temp_tree->SetName( Name );

		temp_trees.push_back( temp_tree );
		minimum_NLL.push_back( temp_tree->GetMinimum( NLL ) );

		NLL_Draw_String = "(" + NLL + "-";
		NLL_Draw_String+=minimum_NLL.back();
		NLL_Draw_String.Append(")");

		_1D_Draw_String = NLL_Draw_String + ":" + Param_Of_Choice + value_suffix;

		ALL_Draw_Strings.push_back( _1D_Draw_String );
		cout << ALL_Draw_Strings.back() << endl;
	}

	TCanvas* new_canvas = new TCanvas("Output", "Output", 1680, 1050 );

	vector<TGraph*> temp_graphs;
	cout << "Constructing Graphs:" << endl;


	//	For thowing away, but certain objects only exist after Canvas creation
	TCanvas* bad_c = new TCanvas("new2", "new2", 1680, 1050);

	for( unsigned int i=0; i< Input_Files.size(); ++i )
	{
		bad_c->cd();
		Input_Tree_per_File[i]->Draw( ALL_Draw_Strings[i], "NLL>0" );
		TGraph* temp_graph = new TGraph( int(Input_Tree_per_File[i]->GetSelectedRows()), Input_Tree_per_File[i]->GetV2(), Input_Tree_per_File[i]->GetV1() );
		TString Name( Param_Of_Choice );
		Name.Append("_");Name+=i;
		temp_graph->SetName( Name );
		temp_graph->SetLineColor( i+1 );
		temp_graph->SetMarkerColor( i+1 );
		new_canvas->cd();
		temp_graphs.push_back( temp_graph );
		mg->Add(temp_graph);
	}

	vector<TString> Drift_Param_Draw_String;

	for( unsigned int i=0; i< all_parameter_values.size(); ++i )
	{
		TString temp_str = all_parameter_values[i];
		temp_str.Append( ":" );
		temp_str.Append( Param_Of_Choice );
		temp_str.Append( value_suffix );
		Drift_Param_Draw_String.push_back( temp_str );
		cout << Drift_Param_Draw_String.back() << endl;
	}

	cout << endl<<"Constructing Nuisence Parameter Plots..." << endl;

	TCanvas* temp_c = new TCanvas( "tmp2", "tmp2", 1680, 1050 );
	//	For all files
	for( unsigned int i=0; i< Input_Files.size(); ++i )
	{
		//	For all fit parameters
		for( unsigned int j=0; j< Drift_Param_Draw_String.size(); ++j )
		{
			Input_Tree_per_File[i]->Draw( Drift_Param_Draw_String[j] , "NLL>0" );
			TGraph* temp_graph = new TGraph( int( Input_Tree_per_File[i]->GetSelectedRows() ), Input_Tree_per_File[i]->GetV2(), Input_Tree_per_File[i]->GetV1() );
			temp_graph->SetLineColor( i+1 );
			temp_graph->SetMarkerColor( i+1 );
			TString Name("File_");
			Name+=j;Name.Append("_");Name+=i;
			temp_graph->SetName( Name );
			all_params_drift_multi[j]->Add( temp_graph );
		}
	}
	delete temp_c;

	new_canvas->cd();
	cout << endl <<"Overlaying And Plotting..." << endl;
	mg->Draw("AC*");

	//	Now that the axis exist we can worry about labeling them
	mg->GetXaxis()->SetTitle( EdStyle::GetParamRootName( Param_Of_Choice ) );
	mg->GetYaxis()->SetTitle( "#Delta LL" );
	new_canvas->Update();

	TString Output_Graph_Name = "Output_"+Param_Of_Choice;
	new_canvas->Print( Output_Graph_Name+".png" );
	new_canvas->Print( Output_Graph_Name+".pdf" );


	if(	want_ROADMAP_ERROR	)
	{
		//	Now mark and draw the 1 sigma error
		double xmax=0, xmin=0;
		double error=0.5;
		//	Zoom in so you can see it
		mg->GetYaxis()->SetRangeUser( 0., 5. );
		//	This was an attempt to get the range to auto adust
		new_canvas->Update();
		xmax = mg->GetXaxis()->GetXmax();
		xmin = mg->GetXaxis()->GetXmin();
		TLine *error_line = new TLine( xmin, error, xmax, error);
		error_line->SetLineColor( Color_t(3) );
		TLine *param_error = NULL;

		//	Definition of errors around CV
		//	Would be nice to be able to read this from somewhere...
		if( Param_Of_Choice == "gamma" )
		{
			param_error = new TLine( 0.6807-0.02, error, 0.6807+.02, error);
		} else if( Param_Of_Choice == "deltaGamma" )
		{
			param_error = new TLine( 0.06-0.06, error, 0.06+0.06, error);
		} else if( Param_Of_Choice == "Azero_sq" )
		{
			param_error = new TLine( 0.6-0.015, error, 0.6+0.015, error);
		} else if( Param_Of_Choice == "Aperp_sq" )
		{
			param_error = new TLine( 0.16-0.02, error, 0.16+0.02, error);
		} else if( Param_Of_Choice == "delta_para" )
		{
			param_error = new TLine( 2.5-0.1, error, 2.5+0.1, error);
		} else if( Param_Of_Choice == "delta_perp" )
		{
			param_error = new TLine( -0.17-0.5, error, -0.17+0.5, error);
		} else if( Param_Of_Choice == "Phi_s" )
		{
			param_error = new TLine( -0.7-0.25, error, -0.7+0.25, error);
		}
		if( param_error != NULL )
		{
			param_error->SetLineColor( Color_t(2) );
			param_error->SetLineWidth( 10 );
			error_line->SetLineWidth( 10 );
			error_line->Draw();
			param_error->Draw();
			new_canvas->Update();
			new_canvas->Print( Output_Graph_Name+"_error.png" );
			new_canvas->Print( Output_Graph_Name+"_error.pdf" );
		}
	}

	for( unsigned int i=0; i < all_parameter_values.size(); ++i )
	{
		TString canvas_name="Canvas_";
		canvas_name+=i;
		TCanvas* out_canvas = new TCanvas(canvas_name, canvas_name, 1680, 1050 );
		all_params_drift_multi[i]->Draw("AC*");
		out_canvas->Update();
		all_params_drift_multi[i]->GetXaxis()->SetTitle( EdStyle::GetParamRootName( Param_Of_Choice ) );
		all_params_drift_multi[i]->GetYaxis()->SetTitle( EdStyle::GetParamRootName( all_parameter_values[i] ) );
		TString Output_Graph_Name = "Output_"+Param_Of_Choice+"_"+all_parameter_values[i];
		out_canvas->Print( Output_Graph_Name+".png");
		out_canvas->Print( Output_Graph_Name+".pdf");
	}

	cout << endl;
	return 0;
}
