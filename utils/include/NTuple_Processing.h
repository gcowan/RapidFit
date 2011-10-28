#ifndef _NTUPLE_PROCESS
#define _NTUPLE_PROCESS

//	ROOT Headers
#include "TFile.h"
#include "TString.h"
#include "TClass.h"
#include "TRandom3.h"
//	System Headers
#include <vector>

using namespace::std;

//	Open a .root file safely i.e. don't crash on failure
TFile* OpenFile( string );

//	Open Multiple .root files safely
vector<TFile*> OpenMultipleFiles( vector<string> );

//	Open a file ane get the requested TTree from the file
TTree* GetTree( string filename, string tuplename );

//	Open Multiple files and get the requested TTrees from the files
vector<TTree*> GetMultipleTrees( vector<string> all_filenames, vector<string> all_tuplenames );

//	Open Multiple files and get the same TTree from each file
vector<TTree*> GetMultipleTrees( vector<string> all_filenames, vector<string> all_tuplenames );

//  From the root current path look for all keys (objects in root file) and loop over them
//  For each one that is actually an object of inherit_type store it's name and the number of events it has
//  Slightly stupid and ineffient as it always has to decend into the file stuture
//  This shouldn't be a huge problem as I never anticipate finding files with hundereds of files
void get_object_list( TString current_path, vector<pair<string,string> > *found_names, TClass* inherit_type, TString *relative_path );


//	Wrapper functions for get_object_list which provides some more clearly named functions with easy to interperate results

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

//      Pass this the name of the ntuple you're handing and it will give you a list
//      of branches that it contains in an easy to handle manner
vector<TString> get_branch_names( TString absolute_path );

//	Create a new tree using a copytree method with a unique name and with the applied Cut_String
TTree* Cut( TTree* input_tree, TString Cut_String, TRandom3* random );




//      Get the data stored in one branch 'Data_String' after applying the cut string
//      This method employs TTree->Draw() which is the fastest method for accessing data in a TTree
vector<Double_t> Get_Data( TTree* input_tree, TString Cut_String, TString Data_String );

//      A wrapper for the vector<double> Get_Data2D for instances where it's easier to write code to
//      Loop over array objects between 1/2/3 D
vector<vector<Double_t> > Get_Datav( TTree* input_tree, TString Cut_String, TString Data_String );

//      2D version of Get_Data
vector<pair<Double_t,Double_t> > Get_Data2D( TTree* input_tree, TString Cut_String, TString Data_String );

//      A wrapper for the vector<pair<double,double> > Get_Data2D for instances where it's easier to write code to
//      Loop over array objects between 1/2/3 D
vector<vector<Double_t> > Get_Data2Dv( TTree* input_tree, TString Cut_String, TString Data_String );

//      3D version of Get_Data
vector<vector<Double_t> > Get_Data3D( TTree* input_tree, TString Cut_String, TString Data_String );


//      Check that the Global Minima as defined at entry 0 within the file is actually the best minima or not
void Check_Minima( TTree* input_tree, TString Cut_String, Float_t* Global_Best_NLL, TString NLL, TString Double_Tolerance, TString param1_val, TString param2_val );

bool Toys_Check( TTree* input_tree, TString Cut_String, TTree* Toy_Tree );

#endif
