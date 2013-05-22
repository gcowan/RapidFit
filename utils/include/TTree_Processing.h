#ifndef _TTree_Processing
#define _TTree_Processing

#include "TString.h"
#include "TTree.h"
#include "TRandom3.h"

#include <vector>
#include <string>

using namespace::std;

class TTree_Processing
{
	public:

		/*!
		 * @brief Pass this the pointer to the ntuple you're handing and it will give you a list
		 *        of branches that it contains in an easy to handle manner
		 *
		 * @return Returns a vector of TStrings which contain the names of the branches in the TTree
		 */
		static vector<TString> get_branch_names( TTree* local_tree );

		//      Pass this the pointer to the ntuple you're handing and it will give you a list
		//      of branches that it contains in an easy to handle manner
		static vector<TString> get_branch_titles( TTree* local_tree );

		//      Pass this the pointer to the ntuple you're handing and it will give you a list
		//      of branches that it contains in an easy to handle manner
		static vector<TString> get_branch_types( TTree* local_tree );

		//      Pass this the name of the ntuple you're handing and it will give you a list
		//      of branches that it contains in an easy to handle manner
		static vector<TString> get_branch_names( TString absolute_path );

		//	Create a new tree using a copytree method with a unique name and with the applied Cut_String
		static TTree* TreeFromCut( TTree* input_tree, TString Cut_String, TRandom* random=NULL );

		//	These functions pass pointer to potentially large vectors to reduce the amount of copy operations I pass pointers, don't like it, don't use it!

		//	Load a whole branch from a TTree into system memory (a Memory intesive operation, but allows the I/O for read/write to be done in 'bursts')
		static vector<Double_t>* Buffer_Branch( TTree* input_tree, TString input_branch, TString cut_string="" );

		//	Load multiple branches and store them in memory
		static vector<vector<Double_t> >* Buffer_Multiple_Branches( TTree* input_tree, TString input_branch, TString cut_string="" );

		//	Return a vector of the number of entries in each ntuple in the list
		static vector<int> get_Event_Number( vector<pair<string,string> > *found_names );

		//	Create a ttree from a vector of vectors
		static TTree* vecvec2TTree( vector<vector<Float_t> > input_vec );

		//	return the *UNIQUE* corrdinates contained in the input_tree from the Draw_String after applying the Cut_String
		//
		//	ATTENTION: THIS WILL MANGLE DATA COMPARED TO HOW IT APPEARS IN THE INPUT TREE, BUT IT GIVES UNIQUE SORTED DATA POINTS SO I DONT CARE!
		static vector<vector<Double_t> > Plotter_Data( TTree* input_tree, TString Draw_String, TString Cut_String, TRandom* random=NULL, int upper_lim=-1 );

		static void AddBranch( TTree* input_tree, string title, vector<double> data );
		static void AddBranch( TTree* input_tree, string title, vector<int> data );
		static void AddBranch( TTree* input_tree, string title, vector<bool> data );
};

#endif

