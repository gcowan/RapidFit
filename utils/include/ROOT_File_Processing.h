#ifndef _ROOTFILE_PROCESS
#define _ROOTFILE_PROCESS

//	ROOT Headers
#include "TFile.h"
#include "TString.h"
#include "TClass.h"
#include "TRandom3.h"
//	System Headers
#include <vector>

using namespace::std;

class ROOT_File_Processing
{
	public:

		//	Open a .root file safely i.e. don't crash on failure
		static TFile* OpenFile( string, string ReadWrite="" );

		//	Open Multiple .root files safely
		static vector<TFile*> OpenMultipleFiles( vector<string> );

		//	Open a file and return the first TTree object (or object derrived from TTree)
		static TTree* GetFirstTree( string filename, string ReadWrite="" );

		//	Open a file ane get the requested TTree from the file
		static TTree* GetTree( string filename, string tuplename );

		//	Open Multiple files and get the requested TTrees from the files
		static vector<TTree*> GetMultipleTrees( vector<string> all_filenames, vector<string> all_tuplenames );

		//	Open Multiple files and get the same TTree from each file
		static vector<TTree*> GetMultipleTrees( vector<string> all_filenames, string tuplename );
		static vector<TTree*> GetMultipleTrees( vector<TString> all_filenames, TString tuplename );

		//  From the root current path look for all keys (objects in root file) and loop over them
		//  For each one that is actually an object of inherit_type store it's name and the number of events it has
		//  Slightly stupid and ineffient as it always has to decend into the file stuture
		//  This shouldn't be a huge problem as I never anticipate finding files with hundereds of files
		static void get_object_list( TString current_path, vector<pair<string,string> > *found_names, TClass* inherit_type, TString *relative_path );


		//	Wrapper functions for get_object_list which provides some more clearly named functions with easy to interperate results

		static void get_TDirectory_list( TString current_path, vector<pair<string,string> > *found_names );

		static void get_TTree_list( TString current_path, vector<pair<string,string> > *found_names );

		static void get_TH1_list( TString current_path, vector<pair<string,string> > *found_names );

		static void get_TH2_list( TString current_path, vector<pair<string,string> > *found_names );

		static void get_TH3_list( TString current_path,vector<pair<string,string> > *found_names );

		static void get_TGraph_list( TString current_path, vector<pair<string,string> > *found_names );

		static void get_TGraph2D_list( TString current_path, vector<pair<string,string> > *found_names );

		static void get_TNtuple_list( TString current_path, vector<pair<string,string> > *found_names );


		//	Wrapper functions to return a list of all parameters 
		static void get_TTree_list_here( vector<pair<string,string> > *found_names );

		static void get_TH1_list_here( vector<pair<string,string> > *found_names );

		static void get_TH2_list_here( vector<pair<string,string> > *found_names );

		static void get_TH3_list_here( vector<pair<string,string> > *found_names );

		static void get_TGraph_list_here( vector<pair<string,string> > *found_names );

		static void get_TGraph2D_list_here( vector<pair<string,string> > *found_names );

		static void get_TNtuple_list_here( vector<pair<string,string> > *found_names );


		//	Open an object which is at pair.first with name pair.second

		static TObject* open_object( pair<string,string> );

};

#endif

