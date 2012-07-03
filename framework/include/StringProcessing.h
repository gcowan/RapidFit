/**
        @class StringProcessing

        Collection of static functions for string processing

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef STRING_PROCESSING_H
#define STRING_PROCESSING_H

//	ROOT Headers
#include "TString.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class StringProcessing
{
	public:
		static string TimeString();
		static vector<string> SplitString( const string, const char );
		static int CharacterPosition( const string, const char );
		static vector<int> StringPositions( const string, const string );
		static void RemoveCharacter( string&, const char );
		static string ReplaceString( const string&, const string, const string );
		static void RemoveWhiteSpace( vector<string>& );

		static string LatexSafe( const string input );
		static string LatexSafe( const TString input );

		/*!
		 * @brief combine 2 vectors of strings into one unique one
		 *
		 * This creates a new vector of strings which is populated first by the unique strings from List1
		 * then appends the list with unique strings from List2
		 *
		 * @param List1  first list of input strings, not assumed to be full of unique strings
		 * @param List2  second list of input strings, not assumed to be full of unique strings
		 *
		 * @return Returns a vector of strings which are all unique
		 */
		static vector<string> CombineUniques( const vector<string> List1, const vector<string> List2 );

		static vector<string> RemoveCommon( const vector<string>, const vector<string> );
		static void RemoveElement( vector<string>& Input, string Element );

		static int VectorContains( vector<string> const*, string const* );
		static int VectorContains( const vector<string>&, const string& );
		static TString CondenseStrings( const vector<string>&, const int, const int );

		static vector<TString> GetStringContaining( const vector<TString>, const TString );
		static vector<TString> StripStrings( const vector<TString>, const TString );
		static vector<string> Convert( const vector<TString> );

                static int GetNumberOnLeft( const string& );
                static string AddNumberToLeft( const string&, const int& );
                static string RemoveFirstNumber( const string& );
                static vector<string> FillList( const int, const int=0 );

                static string AddNames( const string&, const string& );
                static string MultNames( const string&, const string& );

		//	ported from utils / StringOperations
		//	Remove characters not nice in filenames and multiple '_' occurances
		static TString Clean( const TString input );
		//      Is a TString empty ?
                static bool is_empty( const TString input );
};

#endif

