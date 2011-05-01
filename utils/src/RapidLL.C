//	ROOT Headers
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TAxis.h"
//	RapidFit Headers
#include "EdStyle.h"
//	System Headers
#include <vector>
#include <iostream>
#include <cstdlib>

using namespace::std;

//      Pass this the pointer to the ntuple you're handing and it will give you a list
//      of branches that it contains in an easy to handle manner
vector<TString> get_branch_names( TTree* local_tree )
{
        //      To be populated and returned to the user
        vector<TString> temp_branch_names;

        //      Get the list of branches within the TTree
        TObjArray* branch_obj_array = local_tree->GetListOfBranches();

        //      Loop over all found branch objects and request their names
        for( unsigned short int i=0; i < branch_obj_array->GetEntries() ; ++i )
        {
                TObject* branch_object = (*branch_obj_array)[i];
                temp_branch_names.push_back((const char*) branch_object->GetName());
        }

        //      Return the vector of names I have found
        return temp_branch_names;
}

//      Pass this the name of the ntuple you're handing and it will give you a list
//      of branches that it contains in an easy to handle manner
vector<TString> get_branch_names( TString absolute_path )
{
        //      Get the tree that has been defined by the user
        TTree* local_tree = (TTree*) gDirectory->Get((const char*) absolute_path);
        //      Return a vector of Branch names as defined by the user
        return get_branch_names( local_tree );
}

//      Pass this an array of strings and it will find all strings matching a substring
vector<TString> filter_names( vector<TString> all_names, string substring )
{
        //      To be populated and returned to the user
        vector<TString> returnable_names;

        //      Address of the end of the string object
        size_t found=string::npos;

        //      Loop over all strings in the vector
        for( unsigned short int i=0; i<all_names.size(); ++i )
        {       //      Again using STL functions :)
                string temp_str = (all_names[i].Data());
                //      Attempt to find the coordinate of the substring within the string
                found = temp_str.find( substring );
                if( found!=string::npos )       //      If the substring is found
                {
                        returnable_names.push_back(all_names[i]);
                }
        }
        //      Return all strings found containing the substring
        return returnable_names;
}

int main( int argc, char* argv[] )
{
	cout << argv[0] << "\t" << "Param_to_plot" << "\t" << "File1.root" << "\t" << "File2.root" << "\t...\t" << "FileN.root" << endl;

	//	Do this before any Canvas is constructed
	//	ROOT is a BITCH for not allowing you to easily change internal crap
	//      Setup the Canvas and such
	EdStyle* RapidFit_Style = new EdStyle();
	RapidFit_Style->SetStyle();

	TString TREE_NAME( "RapidFitResult" );

	TString Param_Of_Choice( argv[1] );

	TString* Input_File_Names = new TString[size_t(argc-2)];
	vector<TFile*> Input_Files;
	vector<TTree*> Input_Tree_per_File;
	for( int i=2; i < argc ; ++i )
	{
		Input_File_Names[i-2] = TString( argv[i] );

		Input_Files.push_back( new TFile( Input_File_Names[i-2], "READ" ) );

		Input_Tree_per_File.push_back( (TTree*) gDirectory->Get( TREE_NAME ) );

		TString new_name = TREE_NAME;
		new_name+=i;

		//	ROOT does NOT like objects having the same name, keep it from throwing a hissy fit!
		Input_Tree_per_File.back()->SetName( new_name );
	}

	cout << "opened files..." << endl;


	TString value_suffix = "_value";

	//      Get a list of all branches in allresults
	vector<TString> all_parameters = get_branch_names( Input_Tree_per_File[0] );
	//      Get a list of all branches in allresults with '_value' in their name
	vector<TString> all_parameter_values = filter_names( all_parameters, string(value_suffix.Data()) );


	//	Create a vector to hold the overlay graphs for all of the parameters in the fit
	vector<TMultiGraph*> all_params_drift_multi;
	all_params_drift_multi.resize( all_parameters.size() );

	for( unsigned int i=0; i< all_parameters.size(); ++i )
	{
		all_params_drift_multi[i] = new TMultiGraph();
	}

	TMultiGraph *mg = new TMultiGraph();

	TString NLL="NLL";


	vector<double> minimum_NLL( Input_Files.size() );
	vector<TString> ALL_Draw_Strings( Input_Files.size() );

	cout << "Calculating DLL..." << endl;

	TString _1D_Draw_String, NLL_Draw_String;
	for( unsigned int i=0; i< Input_Files.size(); ++i )
	{
		TTree* temp_tree = Input_Tree_per_File[i]->CopyTree( "NLL>0","fast",Input_Tree_per_File[i]->GetEntries(),0 );
		minimum_NLL[i] = temp_tree->GetMinimum( NLL );

		NLL_Draw_String = "(" + NLL + "-";
		NLL_Draw_String+=minimum_NLL[i];
		NLL_Draw_String.Append(")");

		_1D_Draw_String = NLL_Draw_String + ":" + Param_Of_Choice + value_suffix;

		ALL_Draw_Strings[i] = _1D_Draw_String;
	}

	TFile* working_file = new TFile( "/tmp/work.root", "RECREATE" );

	TCanvas* new_canvas = new TCanvas("Output", "Output", 1680, 1050 );

	vector<TGraph*> temp_graphs;
	cout << "Constructing Graphs..." << endl;

	TCanvas* bad_c = new TCanvas("new2", "new2", 1680, 1050);
	for( unsigned int i=0; i< Input_Files.size(); ++i )
	{
//		TCanvas* bad_c = new TCanvas("new2", "new2", 1680, 1050);
		bad_c->cd();
		Input_Tree_per_File[i]->Draw( ALL_Draw_Strings[i], "NLL>0" );
		//TGraph2D* temp_graph = new TGraph2D( Input_Tree_per_File[i]->GetHistogram() );
		TGraph* temp_graph = new TGraph( int(Input_Tree_per_File[i]->GetSelectedRows()), Input_Tree_per_File[i]->GetV2(), Input_Tree_per_File[i]->GetV1() );
		new_canvas->cd();
		temp_graphs.push_back( temp_graph );
		mg->Add(temp_graph);
//		new_canvas->Update();
//		delete bad_c;
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

	cout << "Constructing Nuisence Parameter Plots..." << endl;

	TCanvas* temp_c = new TCanvas( "tmp2", "tmp2", 1680, 1050 );
	//	For all files
	for( unsigned int i=0; i< Input_Files.size(); ++i )
	{
		//	For all fit parameters
		for( unsigned int j=0; j< Drift_Param_Draw_String.size(); ++j )
		{
			Input_Tree_per_File[i]->Draw( Drift_Param_Draw_String[j] , "NLL>0" );
			TGraph* temp_graph = new TGraph( int( Input_Tree_per_File[i]->GetSelectedRows() ), Input_Tree_per_File[i]->GetV2(), Input_Tree_per_File[i]->GetV1() );
			all_params_drift_multi[j]->Add( temp_graph );
		}
	}
//	delete temp_c;

	new_canvas->cd();
	cout << "Overlaying And Plotting..." << endl;
	mg->Draw("AC*");

	//	Now that the axis exist we can worry about labeling them
	mg->GetXaxis()->SetTitle( EdStyle::GetParamRootName( Param_Of_Choice ) );
	mg->GetYaxis()->SetTitle( "#Delta LL" );

	TString Output_Graph_Name = "Output_"+Param_Of_Choice;
	new_canvas->Print( Output_Graph_Name+".png" );
	new_canvas->Print( Output_Graph_Name+".pdf" );


	for( unsigned int i=0; i < all_params_drift_multi.size(); ++i )
	{
		TString canvas_name="Canvas_";
		canvas_name+=i;
		TCanvas* out_canvas = new TCanvas(canvas_name, canvas_name, 1680, 1050 );
		all_params_drift_multi[i]->Draw("AC*");
		all_params_drift_multi[i]->GetXaxis()->SetTitle( EdStyle::GetParamRootName( Param_Of_Choice ) );
		all_params_drift_multi[i]->GetYaxis()->SetTitle( EdStyle::GetParamRootName( all_parameter_values[i] ) );
		TString Output_Graph_Name = "Output_"+Param_Of_Choice+"_"+all_parameter_values[i];
		out_canvas->Print( Output_Graph_Name+".png");
		out_canvas->Print( Output_Graph_Name+".pdf");
//		delete out_canvas;
	}


	//	Clean UP
//	while( !all_params_drift_multi.empty() )
//	{
//		delete all_params_drift_multi.back();
//		all_params_drift_multi.pop_back();
//	}
//	delete mg;
//	delete Input_File_Names;
//	while( !Input_Files.empty() )
//	{
//		delete Input_Files.back();
//		Input_Files.pop_back();
//	}
//	delete new_canvas;

	return 0;
}
