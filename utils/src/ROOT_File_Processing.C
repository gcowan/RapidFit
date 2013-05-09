
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TIterator.h"
#include "TKey.h"

#include "ROOT_File_Processing.h"
#include "Template_Functions.h"
#include "StringOperations.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace::std;

//	Open a .root file without barfing if it doesn't exist, really this should be handled by ROOT SILENTLY && INTERNALLY
TFile* ROOT_File_Processing::OpenFile( string filename, string ReadWrite )
{
	string status = ReadWrite;
	if( ReadWrite.empty() ) status = "READ";

	ifstream input_file;

	input_file.open( filename.c_str(), ifstream::in );
	input_file.close();

	if( !input_file.fail() )
	{
		//      Open the File
		return TFile::Open( filename.c_str(), status.c_str() );
	} else {
		return NULL;
	}

	return NULL;
}

vector<TFile*> ROOT_File_Processing::OpenMultipleFiles( vector<string> all_filenames )
{
	vector<TFile*> output;
	for( vector<string>::iterator file_i=all_filenames.begin(); file_i != all_filenames.end(); ++file_i )
	{
		output.push_back( OpenFile( *file_i ) );
	}
	return output;
}

TTree* ROOT_File_Processing::GetFirstTree( string filename, string ReadWrite )
{
	TFile* input_file = OpenFile( filename, ReadWrite ); (void) input_file;

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

TTree* ROOT_File_Processing::GetTree( string filename, string tuplename )
{
	TFile* input_file = OpenFile( filename );
	TTree* input_tree = NULL;

	//cout << input_file << endl;

	vector<pair<string,string> > all_Trees;

	get_TTree_list_here( &all_Trees );

	vector<string> tree_names = return_second(all_Trees);

	//cout << tuplename << "\t" << all_Trees.size() << endl;
	//print( all_Trees );

	int tuple_index = StringOperations::VectorContains( &tree_names, &tuplename );

	if( tuple_index == -1 ) return NULL;

	string object_path = all_Trees[(unsigned)tuple_index].first;
	object_path+="/"+tuplename;

	input_file->GetObject( object_path.c_str(), input_tree );

	return input_tree;
}

vector<TTree*> ROOT_File_Processing::GetMultipleTrees( vector<string> all_filenames, vector<string> all_tuplenames )
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

vector<TTree*> ROOT_File_Processing::GetMultipleTrees( vector<string> all_filenames, string tuplename )
{
	vector<TTree*> output;
	for( vector<string>::iterator file_i = all_filenames.begin(); file_i != all_filenames.end(); ++file_i )
	{
		output.push_back( GetTree( *file_i, tuplename ) );
	}
	return output;
}

vector<TTree*> ROOT_File_Processing::GetMultipleTrees( vector<TString> all_filenames, TString tuplename )
{

	return GetMultipleTrees( StringOperations::TString2string( all_filenames ), string( tuplename.Data() ) );
}

//  From the root current path look for all keys (objects in root file) and loop over them
//  For each one that is actually an object of inherit_type store it's name and the number of events it has
//  Slightly stupid and ineffient as it always has to decend into the file stuture
//  This shouldn't be a huge problem as I never anticipate finding files with hundereds of files
void ROOT_File_Processing::get_object_list( TString current_path, vector<pair<string,string> > *found_names, TClass* inherit_type, TString *relative_path )
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
			if( StringOperations::VectorContains( &all_found_names, &full_name ) == -1 )
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

void ROOT_File_Processing::get_TDirectory_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TDirectory::Class(), new TString(current_path) );
}

//	Wrapper to find TH1 objects
void ROOT_File_Processing::get_TH1_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TH1::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void ROOT_File_Processing::get_TH2_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TH2::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void ROOT_File_Processing::get_TH3_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TH3::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void ROOT_File_Processing::get_TNtuple_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TNtuple::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void ROOT_File_Processing::get_TTree_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TTree::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void ROOT_File_Processing::get_TGraph_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TGraph::Class(), new TString(current_path) );
}

//      Wrapper to find TH1 objects
void ROOT_File_Processing::get_TGraph2D_list( TString current_path, vector<pair<string,string> > *found_names )
{
	get_object_list( current_path, found_names, TGraph2D::Class(), new TString(current_path) );
}

void ROOT_File_Processing::get_TTree_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TTree_list( current_path, found_names );
}

void ROOT_File_Processing::get_TH1_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TH1_list( current_path, found_names );
}

void ROOT_File_Processing::get_TH2_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TH2_list( current_path, found_names );
}

void ROOT_File_Processing::get_TH3_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TH3_list( current_path, found_names );
}

void ROOT_File_Processing::get_TGraph_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TGraph_list( current_path, found_names );
}

void ROOT_File_Processing::get_TGraph2D_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TGraph2D_list( current_path, found_names );
}

void ROOT_File_Processing::get_TNtuple_list_here( vector<pair<string,string> > *found_names )
{
	TString current_path = gDirectory->GetPath();
	get_TNtuple_list( current_path, found_names );
}

TObject* ROOT_File_Processing::open_object( pair<string,string> details )
{
	TString Original_Path = gDirectory->GetPath();

	gDirectory->cd( details.first.c_str() );

	TObject* object = gDirectory->Get( details.second.c_str() );

	gDirectory->cd( Original_Path );

	return object;
}

