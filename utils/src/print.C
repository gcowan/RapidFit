//	ROOT Headers
#include "TString.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
//	RapidFit Utils Headers
#include "TString_Processing.h"
#include "NTuple_Processing.h"
//	System Headers
#include <vector>
#include <string>
#include <iostream>

using namespace::std;

int main( int argc, char* argv[] )
{
	(void) argc;
	cout << argv[0] << "\t" << "File.root" ;
	cout << endl << endl << "Tool for listing the contents of a .root file" << endl << endl;


	//	Setup the variables that I need
	vector<pair<string,string> > found_names;
	TFile* input = new TFile( TString( argv[1] ), "READ" );


	//	Empty the vector, populate it with all TTree objects in the file,
	//	and, print the results.
	clear( &found_names );
	get_TTree_list_here( &found_names );
	cout << "FOUND: " << found_names.size() << " TTREE OBJECTS:" << endl;
	print( found_names );


	clear( &found_names );
	get_TNtuple_list_here( &found_names );
	cout << endl <<"FOUND: " << found_names.size() << " TNTUPLE OBJECTS:" << endl;
	print( found_names );


	clear( &found_names );
	get_TH1_list_here( &found_names );
	cout << endl <<"FOUND: " << found_names.size() << " TH1 OBJECTS:" << endl;
	print( found_names );


        clear( &found_names );
	get_TH2_list_here( &found_names );
        cout << endl <<"FOUND: " << found_names.size() << " TH2 OBJECTS:" << endl;
	print( found_names );


	clear( &found_names );
	get_TH3_list_here( &found_names );
        cout << endl <<"FOUND: " << found_names.size() << " TH3 OBJECTS:" << endl;
	print( found_names );


	clear( &found_names );
        get_TGraph_list_here( &found_names );
        cout << endl <<"FOUND: " << found_names.size() << " TGraph OBJECTS:" << endl;
	print( found_names );


	clear( &found_names );
        get_TGraph2D_list_here( &found_names );
        cout << endl <<"FOUND: " << found_names.size() << " TGraph2D OBJECTS:" << endl;
	print( found_names );

	cout << endl;

	input->Close();

	return 0;
}
