#ifndef _TSTRING_PROCESS
#define _TSTRING_PROCESS

//	ROOT Headers
#include "TString.h"
#include "TPaveText.h"
//	System Headers
#include <string>
#include <vector>
#include <iostream>

using namespace::std;

//	Add the LHCb label to plots
//	use false at the end to simulation
TPaveText* addLHCbLabel(TString footer, bool DATA=true);

TString remove_suffix( TString whole_string, TString deliminator );

string remove_suffix( string whole_string, string deliminator );

//      Is a TString empty ?
bool is_empty( TString input );

//      Compare a string to a TString
int comp_TString( string test, TString input );

//      Compare a TString to a string
int TString_comp( TString test, string input );

//  Useful for appendding a filename
void filename_append( TString input, TString *output, TString ext );

//      Wrapper
TString filename_append( TString input, TString ext );

string filename_append( string input, string ext );

//	Print a vector of variables which are compatible with cout with a single print command
template<class T> void print ( vector<T> output, bool new_line=true )
{
	cout << "T:" <<output[0]<< endl;
	typename std::vector<T>::iterator output_i;
	for( output_i = output.begin(); output_i != output.end(); ++output_i )
	{
		cout << *output_i ;
		if( new_line ) { cout << endl; }
		else { cout << ",\t" << endl; }
	}
}

//	Print a vector of pairs of variables which are compatible with cout with a single command
template<class T, class U> void print ( vector<pair<T,U> > output, bool new_line=true )
{
        typename std::vector<std::pair<T,U> >::iterator output_i;
        for( output_i = output.begin(); output_i != output.end(); ++output_i )
        {
                cout << output_i->first << "  ,  " << output_i->second ;
                if( new_line ) { cout << endl; }
                else { cout << " ::  " << endl; }
        }
}

//	empty a vector
template<class T> void clear ( vector<T>* input )
{
	while( !input->empty() ) { input->pop_back(); }
}

//	empty a vector of pairs
template<class T, class U> void clear ( vector<pair<T,U> >* input )
{
	while( !input->empty() ) { input->pop_back(); }
}

//      Pass this an array of strings and it will find all strings matching a substring
vector<TString> filter_names( vector<TString> all_names, TString substring );

//	Remove the extention from an array of strings
void strip_strings( vector<string>* list, string ext );

//	Create a new string composed of the 
vector<string> strip_all_strings( vector<string>* list, string ext );

#endif
