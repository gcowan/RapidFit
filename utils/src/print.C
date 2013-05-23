///	ROOT Headers
#include "TString.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
///	RapidFit Utils Headers
#include "TTree_Processing.h"
#include "ROOT_File_Processing.h"
#include "Template_Functions.h"
///	System Headers
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

using namespace::std;

int main( int argc, char* argv[] )
{
	(void) argc;
	cout << argv[0] << "\t" << "File.root" ;
	cout << endl << endl << "Tool for listing the contents of a .root file" << endl << endl;


	//	Setup the variables that I need
	vector<pair<string,string> > found_names;

	for( unsigned int fileNum =1; fileNum < argc; ++fileNum )
	{
		TFile* input = new TFile( TString( argv[fileNum] ), "READ" );


		//	Empty the vector, populate it with all TTree objects in the file,
		//	and, print the results.
		clear( &found_names );
		ROOT_File_Processing::get_TTree_list_here( &found_names );
		cout << "FOUND: " << found_names.size() << " TTREE OBJECTS:" << endl;
		print( found_names );
		for( unsigned int i=0; i< found_names.size(); ++i )
		{
			cout << endl << "For Tree: " << found_names[i].first <<found_names[i].second << ":" << endl;
			gDirectory->cd( found_names[i].first.c_str() );
			TTree* local_tree = (TTree*) gDirectory->Get( found_names[i].second.c_str() );
			cout << "Entries: " << local_tree->GetEntries() << endl;
			vector<TString> branchNames = TTree_Processing::get_branch_names( local_tree );
			vector<TString> branchTitles = TTree_Processing::get_branch_titles( local_tree );
			vector<TString> branchTypes = TTree_Processing::get_branch_types( local_tree );
			cout << "Branches:" << endl;
			for( unsigned int j=0; j< branchNames.size(); ++j )
			{
				cout << "\t\tName: " << setw(20) << branchNames[j] << "\tTitle: " << setw(20) << branchTitles[j] << "\tType: " << branchTypes[j] << endl;
			}
			cout << endl;
		}


		clear( &found_names );
		ROOT_File_Processing::get_TNtuple_list_here( &found_names );
		cout << endl <<"FOUND: " << found_names.size() << " TNTUPLE OBJECTS:" << endl;
		print( found_names );

		clear( &found_names );
		ROOT_File_Processing::get_TH1_list_here( &found_names );
		cout << endl <<"FOUND: " << found_names.size() << " TH1 OBJECTS:" << endl;
		print( found_names );


	        clear( &found_names );
		ROOT_File_Processing::get_TH2_list_here( &found_names );
	        cout << endl <<"FOUND: " << found_names.size() << " TH2 OBJECTS:" << endl;
		print( found_names );


		clear( &found_names );
		ROOT_File_Processing::get_TH3_list_here( &found_names );
	        cout << endl <<"FOUND: " << found_names.size() << " TH3 OBJECTS:" << endl;
		print( found_names );


		clear( &found_names );
	        ROOT_File_Processing::get_TGraph_list_here( &found_names );
	        cout << endl <<"FOUND: " << found_names.size() << " TGraph OBJECTS:" << endl;
		print( found_names );


		clear( &found_names );
	        ROOT_File_Processing::get_TGraph2D_list_here( &found_names );
	        cout << endl <<"FOUND: " << found_names.size() << " TGraph2D OBJECTS:" << endl;
		print( found_names );

		cout << endl;

		input->Close();

		if( fileNum != argc-1 ) cout << endl << "--------------------------------------------------------" << endl << endl;
	}

	return 0;
}
