#ifndef _TSTRING_PROCESS
#define _TSTRING_PROCESS

//	ROOT Headers
#include "TString.h"
#include "TPaveText.h"
//      Utils Headers
#include "Template_Functions.h"
#include "RapidFitParameters.h"
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

//	Append the filename with the deisred extention
string filename_append( string input, string ext );

//      Pass this an array of strings and it will find all strings matching a substring
vector<TString> filter_names( vector<TString> all_names, TString substring );

//	Remove the extention from an array of strings
void strip_strings( vector<string>* list, string ext );

//	Create a new string composed of the 
vector<string> strip_all_strings( vector<string>* list, string ext );

//	Copied functions from within RapidFit StringProcessing to try and guarantee that any changes I make don't break RapidFit
vector<string> SplitString( string, char );
int CharacterPosition( string, char );
vector<int> StringPositions( string, string );
void RemoveCharacter( string&, char );
string ReplaceString( string&, string, string );
void RemoveWhiteSpace( vector<string>& );
vector<string> CombineUniques( vector<string>, vector<string> );
int VectorContains( vector<TString> const&, TString const& );
int VectorContains( vector<string> const*, string const* );
TString CondenseStrings( vector<string>, int, int );

//	Filter the vector TStrings and return a sublist which contains the second argument
vector<TString> GetStringContaining( vector<TString>, TString );
vector<TString> StripStrings( vector<TString>, TString );

//	convert a full vector to the most convenient object to have at the time
vector<string> TString2string( vector<TString> input );
vector<TString> string2TString( vector<string> input );

//	Remove the suffix from the input string and return the string
TString RemoveSuffix( TString input, TString suffix );

TString getFileName( TString input );

TString Clean( TString input );

#endif
