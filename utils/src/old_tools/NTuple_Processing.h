#ifndef _NTUPLE_PROCESS
#define _NTUPLE_PROCESS

//	ROOT Headers
#include "TFile.h"
#include "TString.h"
#include "TClass.h"
#include "TRandom3.h"
//      Utils Headers
#include "Template_Functions.h"
//	System Headers
#include <vector>

using namespace::std;

//	Open a .root file safely i.e. don't crash on failure
TFile* OpenFile( string );

//	Open Multiple .root files safely
vector<TFile*> OpenMultipleFiles( vector<string> );

//	Open a file and return the first TTree object (or object derrived from TTree)
TTree* GetFirstTree( string filename );

//	Open a file ane get the requested TTree from the file
TTree* GetTree( string filename, string tuplename );

//	Open Multiple files and get the requested TTrees from the files
vector<TTree*> GetMultipleTrees( vector<string> all_filenames, vector<string> all_tuplenames );

//	Open Multiple files and get the same TTree from each file
vector<TTree*> GetMultipleTrees( vector<string> all_filenames, string tuplename );

//  From the root current path look for all keys (objects in root file) and loop over them
//  For each one that is actually an object of inherit_type store it's name and the number of events it has
//  Slightly stupid and ineffient as it always has to decend into the file stuture
//  This shouldn't be a huge problem as I never anticipate finding files with hundereds of files
void get_object_list( TString current_path, vector<pair<string,string> > *found_names, TClass* inherit_type, TString *relative_path );


//	Wrapper functions for get_object_list which provides some more clearly named functions with easy to interperate results

void get_TDirectory_list( TString current_path, vector<pair<string,string> > *found_names );

void get_TTree_list( TString current_path, vector<pair<string,string> > *found_names );

void get_TH1_list( TString current_path, vector<pair<string,string> > *found_names );

void get_TH2_list( TString current_path, vector<pair<string,string> > *found_names );

void get_TH3_list( TString current_path,vector<pair<string,string> > *found_names );

void get_TGraph_list( TString current_path, vector<pair<string,string> > *found_names );

void get_TGraph2D_list( TString current_path, vector<pair<string,string> > *found_names );

void get_TNtuple_list( TString current_path, vector<pair<string,string> > *found_names );


//	Wrapper functions to return a list of all parameters 

void get_TTree_list_here( vector<pair<string,string> > *found_names );

void get_TH1_list_here( vector<pair<string,string> > *found_names );

void get_TH2_list_here( vector<pair<string,string> > *found_names );

void get_TH3_list_here( vector<pair<string,string> > *found_names );

void get_TGraph_list_here( vector<pair<string,string> > *found_names );

void get_TGraph2D_list_here( vector<pair<string,string> > *found_names );

void get_TNtuple_list_here( vector<pair<string,string> > *found_names );


//	Open an object which is at pair.first with name pair.second

TObject* open_object( pair<string,string> );



//      Pass this the pointer to the ntuple you're handing and it will give you a list
//      of branches that it contains in an easy to handle manner
vector<TString> get_branch_names( TTree* local_tree );

//      Pass this the pointer to the ntuple you're handing and it will give you a list
//      of branches that it contains in an easy to handle manner
vector<TString> get_branch_titles( TTree* local_tree );

//      Pass this the pointer to the ntuple you're handing and it will give you a list
//      of branches that it contains in an easy to handle manner
vector<TString> get_branch_types( TTree* local_tree );

//      Pass this the name of the ntuple you're handing and it will give you a list
//      of branches that it contains in an easy to handle manner
vector<TString> get_branch_names( TString absolute_path );

//	Create a new tree using a copytree method with a unique name and with the applied Cut_String
TTree* Cut( TTree* input_tree, TString Cut_String, TRandom3* random );



//      Check that the Global Minima as defined at entry 0 within the file is actually the best minima or not
void Check_Minima( TTree* input_tree, TString Cut_String, Float_t* Global_Best_NLL, TString NLL, TString param1_val, TString param2_val );

//	Check for any 'toy studies' exiting within the input file
bool Toys_Check( TTree* input_tree, TString Cut_String, TTree* Toy_Tree );


//	These functions pass pointer to potentially large vectors to reduce the amount of copy operations

//	Load a whole branch from a TTree into system memory (a Memory intesive operation, but allows the I/O for read/write to be done in 'bursts')
vector<Double_t>* Buffer_Branch( TTree* input_tree, TString input_branch, TString cut_string="" );
//	Load multiple branches and store them in memory
vector<vector<Double_t> >* Buffer_Multiple_Branches( TTree* input_tree, TString input_branch, TString cut_string="" );


#endif

