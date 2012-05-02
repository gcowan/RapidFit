//	ROOT Headers
#include "TFile.h"
#include "TTree.h"
//	RapidFit Utils Headers
#include "NTuple_Processing.h"
#include "TString_Processing.h"
//	System Headers
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

using namespace::std;

int main( int argc, char* argv[] )
{

	cout << argv[0] << "\tThis is a tool for adding an index into an ntuple to allow for RapidFit to apply 'cuts' to run over different parts of the data." << endl << endl;

	if( argc != 2 ) exit(-55);

	cout << argv[0] << "\t" << "Input_File.root" << endl;

	cout << "Opening Input File: " << argv[1] <<" and requesting a list of the Trees within it" << endl;

	TFile* in_file = new TFile( TString( argv[1] ), "READ" );

        vector<pair<string,string> > input_trees;

        get_TTree_list_here( &input_trees );

	cout << "Found: " << input_trees.size() << " TTree objects" << endl;

	TString Out_File = filename_append( TString(argv[1]), TString("_indexed") );

	TFile* out_file = new TFile( Out_File, "RECREATE" );

	out_file->SetCompressionLevel( -9 );

	TString Output_File = gDirectory->GetPath();

	cout << "Looping over all found objects and indexing" << endl << endl;

	for( unsigned int i=0; i< input_trees.size(); ++i )
	{
		cout << "Processing:\t" << input_trees[i].second ;
		
		//	goto the directory containing the TTree object
		gDirectory->cd( input_trees[i].first.c_str() );

		TTree* input_tree = (TTree*) gDirectory->Get( input_trees[i].second.c_str() );

		gDirectory->cd( Output_File );

		TTree* output_tree = input_tree->CloneTree( -1 );

		Int_t number=-1;

		TBranch* index_branch = output_tree->Branch( "index", &number, "index/I" );

		for( int j=0; j<( input_tree->GetEntries() ); ++j )
		{
			number = j;
			index_branch->Fill();
		}

		output_tree->Write();

		cout << "\t(Done)" << endl;
	}

	out_file->Close();
	in_file->Close();
	cout << endl << "Finished, Goobye :D" << endl;

}

