//	ROOT Headers
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
//	RapidFit Headers
#include "StringProcessing.h"
//	RapidFit Utils Headers
#include "NTuple_Processing.h"
#include "String_Processing.h"
#include "Histo_Processing.h"
//	System Headers
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>

using namespace::std;

int rand_gen( int _num )
{
	TRandom3* gen = new TRandom3(0);
	return gen->Rndm() * _num;
}

int main( int argc, char* argv[] )
{

	TFile* input_data = new TFile( TString( argv[1] ), "READ" );

	TString file_root = gDirectory->GetPath();

	vector<pair<string,string> > found_names;

	get_TTree_list( file_root, &found_names );

	//	Vector found_names now contains a series of pairs constructed such that
	//	it has the path and name of all the ttree objects in a .root file

	//	Some string which is in the ttree name that identifies the tree as being what I want
	string DataString="Data";

	//	Find all the trees which have this string in their name
	vector<pair<string,string> > filtered = return_pair_check_second( &found_names, DataString );

	//	Nice stdout template function I wrote for a vector of pairs :D
	print( filtered );


	vector<TH1D*> all_TH1 = make_weighted_TH1_Histos( filtered, "X_data", "Y_data" );
	//	Remember that if you call this more than once you will need to rename the objects returned from here


	for( unsigned int i=0; i< all_TH1.size(); ++i )
	{
		TString Canvas_Name="Canvas_Name"; Canvas_Name+=i;
		TCanvas* c1 = EdStyle::RapidFitCanvas( Canvas_Name, Canvas_Name );
		all_TH1[i]->Draw();
		OptimallyRebin( all_TH1[i], 1 );
		c1->Print( Canvas_Name+".pdf" );
		c1->Print( Canvas_Name+".png" );
		c1->Print( Canvas_Name+".C" );
	}


	vector<double> vec_num;
	for( unsigned int i=0; i< 10; ++i ) vec_num.push_back( (double)i );

	vector<double> new_vec = vec_num;

	print( vec_num );

	random_shuffle( new_vec.begin(), new_vec.end(), rand_gen );

	print( vec_num );
	print( new_vec );

	return 0;
}

