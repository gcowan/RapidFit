//	ROOT Headers
#include "TObject.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TDirectory.h"
#include "TString.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLeaf.h"
//	RapidFit Utils Header
#include "EdStyle.h"
#include "TString_Processing.h"
#include "NTuple_Processing.h"
//	System Headers
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using namespace::std;

//	Open a .root file without barfing if it doesn't exist, really this should be handled by ROOT SILENTLY && INTERNALLY
TFile* OpenFile( string filename )
{

	ifstream input_file;

	input_file.open( filename.c_str(), ifstream::in );
	input_file.close();

	if( !input_file.fail() )
	{
		//      Open the File
		return TFile::Open( filename.c_str() );
	} else {
		return NULL;
	}

	return NULL;
}

vector<TFile*> OpenMultipleFiles( vector<string> all_filenames )
{
	vector<TFile*> output;
	for( vector<string>::iterator file_i=all_filenames.begin(); file_i != all_filenames.end(); ++file_i )
	{
		output.push_back( OpenFile( *file_i ) );
	}
	return output;
}

TTree* GetFirstTree( string filename )
{
	TFile* input_file = OpenFile( filename ); (void) input_file;

	vector<pair<string,string> > ttree_names;

	get_TTree_list_here( &ttree_names );

	if( ttree_names.size() != 1 )
	{
		cerr << "WARNING: Mutliple Trees found in file: " << filename << ", you're probably expecting just 1 tree when using the GetFirstTree function" << endl;
	}

	TTree* input_tree = NULL;

	gDirectory->cd( ttree_names[0].first.c_str() );

	input_tree = (TTree*) gDirectory->Get( ttree_names[0].second.c_str() );

	return input_tree;
}

TTree* GetTree( string filename, string tuplename )
{
	TFile* input_file = OpenFile( filename );
	TTree* input_tree = NULL;

	//cout << input_file << endl;

	vector<pair<string,string> > all_Trees;

	get_TTree_list_here( &all_Trees );

	vector<string> tree_names = return_second(all_Trees);

	//cout << tuplename << "\t" << all_Trees.size() << endl;
	//print( all_Trees );

	int tuple_index = VectorContains( &tree_names, &tuplename );

	if( tuple_index == -1 ) return NULL;

	string object_path = all_Trees[(unsigned)tuple_index].first;
	object_path+="/"+tuplename;

	input_file->GetObject( object_path.c_str(), input_tree );

	return input_tree;
}

vector<TTree*> GetMultipleTrees( vector<string> all_filenames, vector<string> all_tuplenames )
{
	vector<TTree*> output;
	if( all_filenames.size() != all_tuplenames.size() )
	{
		cerr << "ERROR:\tReuested number of files different to number of tuples wanted!!!" << endl;
		return output;
	}
	vector<string>::iterator file_i=all_filenames.begin();
	vector<string>::iterator tuple_i=all_tuplenames.begin();
	for( ; file_i != all_filenames.end(); ++file_i, ++tuple_i )
	{
		output.push_back( GetTree( *file_i, *tuple_i ) );
	}
	return output;
}

vector<TTree*> GetMultipleTrees( vector<string> all_filenames, string tuplename )
{
	vector<TTree*> output;
	for( vector<string>::iterator file_i = all_filenames.begin(); file_i != all_filenames.end(); ++file_i )
	{
		output.push_back( GetTree( *file_i, tuplename ) );
	}
	return output;
}

vector<TTree*> GetMultipleTrees( vector<TString> all_filenames, TString tuplename )
{

	return GetMultipleTrees( TString2string( all_filenames ), string( tuplename.Data() ) );
}

//  From the root current path look for all keys (objects in root file) and loop over them
//  For each one that is actually an object of inherit_type store it's name and the number of events it has
//  Slightly stupid and ineffient as it always has to decend into the file stuture
//  This shouldn't be a huge problem as I never anticipate finding files with hundereds of files
void get_object_list( TString current_path, vector<pair<string,string> > *found_names, TClass* inherit_type, TString *relative_path )
{
	//  goto current search path
	gDirectory->cd( current_path );
	TDirectory *current_path_directory = gDirectory;
	TString current_path_str( current_path );

	//  get all keys present in current path
	TIter nextkey( current_path_directory->GetListOfKeys() );
	TKey *key, *oldkey=0;

	//cout << current_path << endl;
	while ( ( key = (TKey*) nextkey() ) )
	{
		//  Loop over all keys (objects) in the current path
		//  Leave on the last key

		if ( oldkey && !strcmp( oldkey->GetName(), key->GetName() ) ) continue;
		//  Should already be here?
		//gDirectory->cd( current_path );
		//  Read in the object
		TObject *obj = key->ReadObj();

		//  Get the name of the object in current directory
		string obj_name( (char*) obj->GetName() );
		TString temp( obj_name );

		//  If this is a Tuple object we want to store some info about it
		if ( obj->IsA()->InheritsFrom( inherit_type ) )
		{
			//  Store the full name of the object and it's full path to determine if we've already found it
			string full_name( relative_path->Data() );
			full_name.append( obj_name );

			vector<string> all_found_names;

			for( vector<pair<string,string> >::iterator found_i=found_names->begin(); found_i!=found_names->end(); ++found_i )
			{
				string full_name = found_i->first;
				full_name.append( string("/") );
				full_name.append( found_i->second );
				all_found_names.push_back( full_name );
			}

			//  If this is an object of the wanted type we foud ignore it
			if( VectorContains( &all_found_names, &full_name ) == -1 )
			{
				//  If this is a new class object save it's name and event_number
				//relative_path->Append( obj_name );
				//found_names->push_back( (*relative_path).Data() );
				pair<string,string> path_and_name;
				path_and_name.first = string( (*relative_path).Data() );
				path_and_name.second = obj_name;
				found_names->push_back( path_and_name );
				continue;
			}
		} else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {

			//  If this is a directory we want to decend down into it
			current_path_str.Append( obj_name );
			current_path_str.Append( "/" );
			relative_path->Append( obj_name );
			relative_path->Append( "/" );
			get_object_list( current_path_str, found_names, inherit_type, relative_path );
		}
	}
}

void get_TDirectory_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TDirectory::Class(), new TString(current_path) );
}

//	Wrapper to find TH1 objects
void get_TH1_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TH1::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void get_TH2_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TH2::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void get_TH3_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TH3::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void get_TNtuple_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TNtuple::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void get_TTree_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TTree::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void get_TGraph_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TGraph::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void get_TGraph2D_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TGraph2D::Class(), new TString(current_path) );
}

void get_TTree_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TTree_list( current_path, found_names );
}

void get_TH1_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TH1_list( current_path, found_names );
}

void get_TH2_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TH2_list( current_path, found_names );
}

void get_TH3_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TH3_list( current_path, found_names );
}

void get_TGraph_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TGraph_list( current_path, found_names );
}

void get_TGraph2D_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TGraph2D_list( current_path, found_names );
}

void get_TNtuple_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TNtuple_list( current_path, found_names );
}

vector<int> get_Event_Number( vector<pair<string,string> > *found_names )
{
	vector<int> returnable_list;
	TString current_path = gDirectory->GetPath();
	get_TTree_list( current_path, found_names );
	for( vector<pair<string,string> >::iterator found_i = found_names->begin(); found_i != found_names->end(); ++found_i )
	{
		gDirectory->cd( found_i->first.c_str() );

		TTree* temp_tree = (TTree*) gDirectory->Get( found_i->second.c_str() );

		returnable_list.push_back( (int) temp_tree->GetEntries() );
	}
	return returnable_list;
}

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

//      Pass this the pointer to the ntuple you're handing and it will give you a list
//      of branches that it contains in an easy to handle manner
vector<TString> get_branch_titles( TTree* local_tree )
{
	//      To be populated and returned to the user
	vector<TString> temp_branch_names;

	//      Get the list of branches within the TTree
	TObjArray* branch_obj_array = local_tree->GetListOfBranches();

	//      Loop over all found branch objects and request their names
	for( unsigned short int i=0; i < branch_obj_array->GetEntries() ; ++i )
	{
		TObject* branch_object = (*branch_obj_array)[i];
		temp_branch_names.push_back((const char*) branch_object->GetTitle());
	}

	//      Return the vector of names I have found
	return temp_branch_names;
}

//      Pass this the pointer to the ntuple you're handing and it will give you a list
//      of branches that it contains in an easy to handle manner
vector<TString> get_branch_types( TTree* local_tree )
{
	//      To be populated and returned to the user
	vector<TString> temp_branch_names;

	//      Get the list of branches within the TTree
	TObjArray* branch_obj_array = local_tree->GetListOfBranches();

	//      Loop over all found branch objects and request their names
	for( unsigned short int i=0; i < branch_obj_array->GetEntries() ; ++i )
	{
		TBranch* branch_object = (TBranch*)(*branch_obj_array)[i];
		if( branch_object != NULL )
		{
			TLeaf* leaf_object = (TLeaf*) branch_object->GetListOfLeaves()->First();
			temp_branch_names.push_back( leaf_object->GetTypeName() );
		}
		else
		{
			temp_branch_names.push_back("NULL");
		}
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


TObject* open_object( pair<string,string> details )
{
	TString Original_Path = gDirectory->GetPath();

	gDirectory->cd( details.first.c_str() );

	TObject* object = gDirectory->Get( details.second.c_str() );

	gDirectory->cd( Original_Path );

	return object;
}


//      Perform a cut on a TTree object and return a new TTree with a unique name
TTree* Cut( TTree* input_tree, TString Cut_String, TRandom3* random )
{
	//      Tree to be returned
	TTree* wanted_tree = NULL;
	//      Perform The Cut
	wanted_tree = input_tree->CopyTree( Cut_String, "fast", input_tree->GetEntries() );

	//      Generate a Unique Name
	TString Name="Tree_";
	double rand = random->Rndm();
	Name+=rand;

	//      Set Unique Name
	wanted_tree->SetName(Name);
	wanted_tree->SetTitle(Name);

	//      Return the TTree
	return wanted_tree;
}


//      Check that the Global Minima as defined at entry 0 within the file is actually the best minima or not
void Check_Minima( TTree* input_tree, TString Cut_String, Float_t* Global_Best_NLL, TString NLL, TString param1_val, TString param2_val )
{
	TTree* wanted = input_tree->CopyTree( Cut_String, "fast", input_tree->GetEntries() );

	Float_t Global_Min_NLL = Float_t( wanted->GetMinimum( NLL ) );

	//      Strings that are universally defined
	TString Double_Tolerance="0.000001";

	//      Only if we have a better minima Minuit did NOT find
	if( (Global_Min_NLL-*Global_Best_NLL) < 0 )
	{
		TString Global_Min_NLL_Str;
		Global_Min_NLL_Str+=Global_Min_NLL;
		TString Catch( "abs(" + NLL + "-" + Global_Min_NLL_Str + ")<" + Double_Tolerance);
		TTree* local_best = wanted->CopyTree( Catch, "fast", wanted->GetEntries() );
		//      GetMinimum == GetMaximum == Get 0th event
		double true_X = local_best->GetMinimum(param1_val);
		double true_Y = local_best->GetMinimum(param2_val);
		cout << "\n\t\tWARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING\n"<<endl;
		cout << "\t\tTRUE MINIMUM NLL = " << Global_Min_NLL << "\t\tAt:\t\tX:" << true_X << "\tY:\t" << true_Y <<endl<<endl;
		//              cout << "\tNEW MINIMA FOUND AT :\tX:\t" << X_true_min << "\tY:\t" << Y_true_min << endl<<endl;
		cout << "\t\tWARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING!!!WARNING\n"<<endl;
		*Global_Best_NLL=Global_Min_NLL;
	}

}


//	Not sure why I wrote this I think it can be safely removed in future as I appear not to use it
//vector<vector<Double_t> > Get_All_Data( TTree* input_tree, TString Cut_String, TString Data_String )
//{
//	int dimention_num = 1;
//	string input_str = Data_String.Data();
//	for( size_t found=0; found!=string::npos; )
//	{
//		found = input_str.find( ":", found+1, string::npos );
//		if( found!=string::npos ) ++dimention_num;
//	}
//	cout << dimention_num << endl;
//	input_tree->SetEstimate(input_tree->GetEntries());
//	input_tree->Draw( Data_String, Cut_String );
//	int Data_Num = int( input_tree->GetSelectedRows() );
//	Double_t* temp_pointer = input_tree->GetV1();
//	vector<vector<Double_t> > Returnable_Data;
//	for( int i=0; i < Data_Num; ++i )
//	{
//		vector<Double_t> temp_vec;
//		temp_vec.push_back( double(temp_pointer[i]) );
//		Returnable_Data.push_back( temp_vec );
//	}
//	return Returnable_Data;
//}

bool Toys_Check( TTree* input_tree, TString Cut_String, TTree* Toy_Tree )
{
	Toy_Tree = input_tree->CopyTree( Cut_String, "fast", input_tree->GetEntries() );
	if( Toy_Tree->GetEntries() > 0 ) return true;
	return false;
}



vector<Double_t>* Buffer_Branch( TTree* input_tree, TString input_branch, TString cut_string )
{
	vector<Double_t>* output_branch_data = new vector<Double_t>();

	//      By default ROOT limits this to 1M entries, this has been written explicitly for some cases when this is insufficient
	input_tree->SetEstimate( input_tree->GetEntries() );

	//      Use TSelector to accesst the data of interest
	input_tree->Draw( input_branch, cut_string, "goff" );

	Double_t* temp_pointer = input_tree->GetV1();

	for( unsigned int entry_i=0; entry_i< input_tree->GetSelectedRows(); ++entry_i )
	{
		output_branch_data->push_back( temp_pointer[entry_i] );
	}

	return output_branch_data;
}

vector<vector<Double_t> >* Buffer_Multiple_Branches( TTree* input_tree, TString input_branch, TString cut_string )
{
	vector<vector<Double_t> >* output_branch_data = new vector<vector<Double_t> >();

	//      By default ROOT limits this to 1M entries, this has been written explicitly for some cases when this is insufficient
	input_tree->SetEstimate( input_tree->GetEntries() );

	//      Use TSelector to accesst the data of interest
	input_tree->Draw( input_branch, cut_string, "goff" );

	//	Draw actually done inside Player object constucted by the TTree
	//	TPlayer creates a Selector instance which stored the number of dimensions in the 'input_branch' statement
	int numberofinputColumns = input_tree->GetPlayer()->GetDimension();

	//	Move the objects into our existance
	for( int i=0; i< numberofinputColumns; ++i )
	{
		vector<Double_t>* thisColumn = new vector<Double_t>();

		Double_t* temp_pointer = input_tree->GetVal( i );

		for( unsigned int entry_i=0; entry_i< input_tree->GetSelectedRows(); ++entry_i )
		{
			thisColumn->push_back( temp_pointer[entry_i] );
		}

		output_branch_data->push_back( *thisColumn );
	}
	return output_branch_data;
}

/*
vector<Double_t>* GetUnsortedData( TTree* input_tree, TString input_branch, TString cut_string, int wanted_entries, int start )
{
	if( wanted_entries = -1 ) wanted_entries = input_tree->GetEntries();
	TTree* temp_tree = input_tree->CopyTree( cut_string, "", wanted_entries, start );

	int new_total = temp_tree->GetEntries();

	vector<Double_t>*

	TBranch* wanted_branch = temp_tree->

	TString branch_type = 

	for( int i=0; i != new_total; ++i )
	{
		
	}

}
*/

