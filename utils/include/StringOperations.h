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

class StringOperations
{
	public:

		//	Return a vector which contains all of the input TStrings postpended to have the requested suffix
		static vector<TString> postpend( vector<TString> input, TString suffix );

		//	Add the LHCb label to plots
		//	use false at the end to simulation
		static TPaveText* addLHCbLabel(TString footer, bool DATA=true);

		//	Remove the suffix which is after the "deliminator" character
		//	For multiple occurances of the character the last on in the whitrong is chosen
		static TString remove_suffix( TString whole_string, TString deliminator );
		static string remove_suffix( string whole_string, string deliminator );

		//	      Is a TString empty ?
		static bool is_empty( TString input );

		//      Compare a string to a TString
		static int comp_TString( string test, TString input );

		//      Compare a TString to a string
		static int TString_comp( TString test, string input );

		//  Useful for appendding a filename
		static void filename_append( TString input, TString *output, TString ext );

		//      Wrapper
		static TString filename_append( TString input, TString ext );

		//	Append the filename with the deisred extention
		static string filename_append( string input, string ext );

		//      Pass this an array of strings and it will find all strings matching a substring
		static vector<TString> filter_names( vector<TString> all_names, TString substring );

		//	Remove the extention from an array of strings
		static void strip_strings( vector<string>* list, string ext );

		//	Create a new string composed of the 
		static vector<string> strip_all_strings( vector<string>* list, string ext );

		//	Copied functions from within RapidFit StringProcessing to try and guarantee that any changes I make don't break RapidFit
		static vector<string> SplitString( string, char );
		static int CharacterPosition( string, char );
		static vector<int> StringPositions( string, string );
		static void RemoveCharacter( string&, char );
		static string ReplaceString( string&, string, string );
		static void RemoveWhiteSpace( vector<string>& );
		static vector<string> CombineUniques( vector<string>, vector<string> );
		static int VectorContains( vector<TString> const&, TString const& );
		static int VectorContains( vector<string> const*, string const* );
		static int VectorContains( const vector<string>& inputVec, const string& inputStr );
		static TString CondenseStrings( vector<string>, int, int );

		//	Filter the vector TStrings and return a sublist which contains the second argument
		static vector<TString> GetStringContaining( vector<TString>, TString );
		static vector<TString> StripStrings( vector<TString>, TString );

		//	convert a full vector to the most convenient object to have at the time
		static vector<string> TString2string( vector<TString> input );
		static vector<TString> string2TString( vector<string> input );

		//	Remove the suffix from the input string and return the string
		static TString RemoveSuffix( TString input, TString suffix );

		static TString getFileName( TString input );

		static TString Clean( TString input );

		static TString prettyPrint(Double_t value);

		static TString prettyPrint(Double_t val, Double_t err);

		static vector<TString> filter_names( vector<TString> all_names, string substring );


		//	Return a flat string of the current time
		static string TimeString();
};

#endif
