
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "ROOT_File_Processing.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "Mathematics.h"

#include "EdStyle.h"
#include "TTree_Processing.h"

#include <vector>
#include <string>
#include <iostream>

using namespace::std;

vector<int> TTree_Processing::get_Event_Number( vector<pair<string,string> > *found_names )
{
	vector<int> returnable_list;
	TString current_path = gDirectory->GetPath();
	ROOT_File_Processing::get_TTree_list( current_path, found_names );
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
vector<TString> TTree_Processing::get_branch_names( TTree* local_tree )
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
vector<TString> TTree_Processing::get_branch_titles( TTree* local_tree )
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
vector<TString> TTree_Processing::get_branch_types( TTree* local_tree )
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
vector<TString> TTree_Processing::get_branch_names( TString absolute_path )
{
	//      Get the tree that has been defined by the user
	TTree* local_tree = (TTree*) gDirectory->Get((const char*) absolute_path);
	//      Return a vector of Branch names as defined by the user
	return get_branch_names( local_tree );
}

//      Perform a cut on a TTree object and return a new TTree with a unique name
TTree* TTree_Processing::TreeFromCut( TTree* input_tree, TString Cut_String, TRandom* random )
{
	if( random == NULL ) random=gRandom;
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

vector<Double_t>* TTree_Processing::Buffer_Branch( TTree* input_tree, TString input_branch, TString cut_string )
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

vector<vector<Double_t> >* TTree_Processing::Buffer_Multiple_Branches( TTree* input_tree, TString input_branch, TString cut_string )
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

vector<vector<Double_t> > TTree_Processing::Plotter_Data( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom* random, int upper_lim )
{
	if( random == NULL ) random = gRandom;
	if( upper_lim == -1 ) upper_lim = input_tree->GetEntries();
	//      Plot the graph using TTree->Draw()
	//      The resulting graph (TH3) contains empty axis and a TPolyMarker3D object
	input_tree->SetEstimate(input_tree->GetEntries());  // Fix the size of the array of doubles to be created (There will never be more than this
	TDirectory* temp_path = gDirectory;
	TString canv_name = "Plot_Data_"; canv_name += random->Rndm()*1000;
	TCanvas* temp_canv = EdStyle::RapidFitCanvas( canv_name, canv_name );
	//	Cannot implement goff in the case of 3D data that we wish to be in TPolyMarker3D form
	input_tree->Draw( Draw_String, Cut_String, "", upper_lim, 0 );

	int numberofinputColumns = input_tree->GetPlayer()->GetDimension();

	if( numberofinputColumns == 3 )
	{
		//      Get the Points that have been plotted (TPolyMarker3D object named "TPolyMarker3D", see ROOTtalk)
		TPolyMarker3D *pm = new TPolyMarker3D( *(TPolyMarker3D*)gPad->FindObject("TPolyMarker3D") );
		double temp = random->Rndm();
		//temp_canv->Close();
		//temp_path->cd();
		TString Name = "TPoly3_";
		Name+=temp;

		if( pm!= NULL )
		{
			pm->SetName(Name);

			//      Get a list of ONLY unique coordinates due to the short comings of the interpolation within TGraph2D
			vector<vector<Double_t> > returnable_data = Mathematics::Unique_Coords( pm );

			return returnable_data;
		}
		cout << "Didn't pick up any points in Cut" << endl;
		vector<vector<Double_t> > dummy_return;
		return dummy_return;
	}
	else if ( numberofinputColumns == 2 )
	{
		Double_t* first = input_tree->GetV1();
		Double_t* second = input_tree->GetV2();

		int rows = (int)input_tree->GetSelectedRows();

		vector<pair<Double_t,Double_t> > to_be_sorted;

		for( int i=0; i< rows; ++i )
		{
			to_be_sorted.push_back( make_pair( first[(unsigned)i], second[(unsigned)i] ) );
		}

		vector<vector<Double_t> > returnable_data = Mathematics::Unique_Coords( to_be_sorted );
		return returnable_data;
	}
	else
	{
		Double_t* first = input_tree->GetV1();

		int rows = (int)input_tree->GetSelectedRows();

		vector<Double_t> to_be_sorted;

		for( int i=0; i< rows; ++i )
		{
			to_be_sorted.push_back( first[(unsigned)i] );
		}

		vector<vector<Double_t> > returnable_data = Mathematics::Unique_Coords( to_be_sorted );
		return returnable_data;
	}
}

//	Return a TTree object composed from a vector of vector of data objects, I provide a stupid branch naming scheme
TTree* TTree_Processing::vecvec2TTree( vector<vector<Float_t> > input_vec )
{
	TTree* new_tree = new TTree( "tree2Draw", "tree2Drw" );

	Float_t* Float_data = new Float_t[ input_vec[0].size() ];
	for( unsigned int i=0; i<input_vec[0].size(); ++i )
	{
		TString br_name("Branch_");
		br_name+=i;
		TString br_title(br_name);
		br_title.Append("/F");
		new_tree->Branch( br_name, &Float_data[i], br_title );
	}

	for( unsigned int i=0; i<input_vec.size(); ++i)
	{
		for( unsigned int j=0; j<input_vec[i].size(); ++j )
		{
			Float_data[j] = input_vec[i][j];
			new_tree->Fill();
		}
	}

	return new_tree;
}

void TTree_Processing::AddBranch( TTree* input_tree, string title, vector<double> data )
{
	if( (input_tree->GetEntries() != (int)data.size()) && (input_tree->GetEntries() != 0) )
	{
		cout << "ERROR: CANNOT ADD BRANCH " << title << " DATA SIZE DID NOT MATCH TREE." << endl;
		return;
	}

	unsigned int size = (unsigned)data.size();

	double thisValue=0.;
	TString name=title.c_str();
	name.Append("/D");
	TBranch* output_branch_obj = input_tree->Branch( title.c_str(), &thisValue, name );

	input_tree->SetEntries( (int)size );

	for( unsigned int i=0; i< size; ++i )
	{
		thisValue = data[i];
		output_branch_obj->Fill();
	}
}

void TTree_Processing::AddBranch( TTree* input_tree, string title, vector<int> data )
{
	if( (input_tree->GetEntries() != (int)data.size()) && (input_tree->GetEntries() != 0) )
	{
		cout << "ERROR: CANNOT ADD BRANCH " << title << " DATA SIZE DID NOT MATCH TREE." << endl;
		return;
	}

	unsigned int size = (unsigned)data.size();

	int thisValue=0;
	TString name=title.c_str();
	name.Append("/I");
	TBranch* output_branch_obj = input_tree->Branch( title.c_str(), &thisValue, name );

	input_tree->SetEntries( (int)size );

	for( unsigned int i=0; i< size; ++i )
	{
		thisValue = data[i];
		output_branch_obj->Fill();
	}
}

void TTree_Processing::AddBranch( TTree* input_tree, string title, vector<bool> data )
{
	if( (input_tree->GetEntries() != (int)data.size()) && (input_tree->GetEntries() != 0) )
	{
		cout << "ERROR: CANNOT ADD BRANCH " << title << " DATA SIZE DID NOT MATCH TREE." << endl;
		return;
	}

	unsigned int size = (unsigned)data.size();

	bool thisValue=false;
	TString name=title.c_str();
	name.Append("/B");
	TBranch* output_branch_obj = input_tree->Branch( title.c_str(), &thisValue, name );

	input_tree->SetEntries( (int)size );

	for( unsigned int i=0; i< size; ++i )
	{
		thisValue = data[i];
		output_branch_obj->Fill();
	}
}

