#ifndef _STRING_PROCESS
#define _STRING_PROCESS

//	ROOT Headers
#include "TString.h"
#include "TPaveText.h"
//	System Headers
#include <string>
#include <vector>
#include <iostream>

using namespace::std;

vector<TString> postpend( vector<TString> input, TString suffix );

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

//      Print a vector of pairs of variables which are compatible with cout with a single command
template<class T, class U> std::vector<T> return_first( vector<pair<T,U> > output )
{
	typename std::vector<T> output_vec;
	typename std::vector<std::pair<T,U> >::iterator output_i;
	for( output_i = output.begin(); output_i != output.end(); ++output_i )
	{
		output_vec.push_back( output_i->first );
	}
	return output_vec;
}

//      Print a vector of pairs of variables which are compatible with cout with a single command
template<class T, class U> std::vector<U> return_second( vector<pair<T,U> > output )
{
	typename std::vector<U> output_vec;
	typename std::vector<std::pair<T,U> >::iterator output_i;
	for( output_i = output.begin(); output_i != output.end(); ++output_i )
	{
		output_vec.push_back( output_i->second );
	}
	return output_vec;
}

//      Print a vector of pairs of variables which are compatible with cout with a single command
template<class T, class U> std::vector<U> swap( vector<pair<T,U> > input )
{
	typename std::vector<std::pair<U,T> > output_vec;
	typename std::vector<std::pair<T,U> >::iterator input_i;
	for( input_i = input.begin(); input_i != input.end(); ++input_i )
	{
		std::pair<T,U> new_pair;
		new_pair.first = input_i->second;
		new_pair.second = input_i->first;
		output_vec.push_back( new_pair );
	}
	return output_vec;
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

vector<string> strip_all_first_strings( vector<pair<string,string> >* list, string ext );

vector<string> return_all_second_strings( vector<pair<string,string> >* list, string ext );

vector<string> return_all_second_strings( vector<pair<string,string> >* list, string ext );

vector<pair<string,string> > return_pair_check( vector<pair<string,string> >* list, vector<string> new_list, vector<string> final_list, string ext );

vector<pair<string,string> > return_pair_check_first( vector<pair<string,string> >* list, string ext );

vector<pair<string,string> > return_pair_check_second( vector<pair<string,string> >* list, string ext );

#endif
