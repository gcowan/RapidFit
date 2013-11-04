//#define isMC

//#define SSpion // uses SSpion tagger (Bu and Bd) ELSE SSkaon (Bs only)
//	To Compile:
//	
//	g++ -Wall -Wextra -Wsign-compare -Wmissing-noreturn -msse -msse2 -m3dnow -g -ansi -O3 -funroll-all-loops `root-config --cflags --libs --ldflags` -lTreePlayer -Wall -g CreateNtupleB.C -o CreateNtupleB
//	./CreateNtupleB ntuplaA.root B_s0
//	
//	(( to Run in CINT:
//	(( root -q -b "CreateNtupleB.C(\"ntuplaA.root", \"B_s0\")"
//	
//     WARNING, M. Dorigo and O. Leroy 8 October 2013:
//     tagdecision and tagomega (os, ss and total) contains by default the nnetKaon
//     	
//
//	ROOT Headers
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFormula.h"
#include "TString.h"
#include "TTreeFormula.h"
//	System Headers
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <cstdlib>

using namespace::std;

template <typename T> vector<T> Buffer_Branch( TTree* input_tree, TString input_branch, TString cut_string );

void DumpBranch( TString input_file, TString input_branch, TString cut_string, TString output_file );

void DumpBranch( TString input_file, TString input_tree, TString input_branch, TString cut_string, TString output_file )
{
	TFile* thisFile = TFile::Open( input_file, "READ" );

	if( thisFile == NULL )
	{
		cerr << "File Not Found" << endl;
	}

	TTree* thisTree = (TTree*) thisFile->Get( input_tree );

	if( thisTree == NULL )
	{
		cerr << "Tree Not Found" << endl;
		exit(0);
	}

	vector<double> thisBranchData = Buffer_Branch<double>( thisTree, input_branch, cut_string );

	if( thisBranchData.empty() )
	{
		cerr << "No Branch Data Found" << endl;
		exit(0);
	}

	stringstream outputFileContent;

	for( unsigned int i=0; i< thisBranchData.size(); ++i )
	{
		outputFileContent << thisBranchData[i] << endl;
	}

	ofstream outFile;
	outFile.open( output_file );
	outFile << outputFileContent.str();
	outFile.close();
}

//#	Window to the code when compiled and run externally to ROOT rather than in CINT
#ifndef __CINT__
int main( int argc, char* argv[] )
{
	//	Give standard Usage information and exit if we have the wrong number of inputs
	cout << "Usage:" << " " << argv[0] << " " << "Input_File.root  inputTree  inputBranch  CutString  OutputFile.txt" << endl;
	if( argc != 6 ) {
		return(-1);
	}

	cout << "Dumping: " << argv[1] << endl;

	string file = argv[1];
	string tree = argv[2];
	string branch = argv[3];
	string cut = argv[4];
	string output = argv[5];
	DumpBranch( file.c_str(), tree.c_str(), branch.c_str(), cut.c_str(), output.c_str() );

}
#endif

	template <typename T>
vector<T> Buffer_Branch( TTree* input_tree, TString input_branch, TString cut_string )
{
	vector<T> output_branch_data;

	//      By default ROOT limits this to 1M entries, this has been written explicitly for some cases when this is insufficient
	input_tree->SetEstimate( input_tree->GetEntries() );

	//      Use TSelector to accesst the data of interest
	input_tree->Draw( input_branch, cut_string, "goff" );

	Double_t* temp_pointer = input_tree->GetV1();

	for(Long64_t entry_i = 0; entry_i < input_tree->GetSelectedRows(); ++entry_i) {
		output_branch_data.push_back( T(temp_pointer[entry_i]) );
	}
	return output_branch_data;
}

