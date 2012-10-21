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
		//	Return the current time as a time-stamp string
		static string TimeString();

		//	Split a string into a vector of strings by the deliminator
		static vector<string> SplitString( const string, const char );

		//	Get the Position of a Character in a string (first occurance)
		static int CharacterPosition( const string, const char );

		//	Get all of the string Positions of a chatacter in a string
		static vector<int> StringPositions( const string, const string );

		//	Remove the first occurance? of a selected character
		static void RemoveCharacter( string&, const char );

		//	Replace an instance of a substring within a 
		static string ReplaceString( const string&, const string, const string );
		static void RemoveWhiteSpace( vector<string>& );

		//	Replace non-Latex safe characters with '_'
		static string LatexSafe( const string input );

		//	Replace non-Latex safe chatacters with '_'
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

		static vector<string> RemoveDuplicates( const vector<string>, vector<string>& duplicated );
		static vector<string> RemoveCommon( const vector<string>, const vector<string> );
		static void RemoveElement( vector<string>& Input, string Element );

		static int VectorContains( vector<string> const*, string const* );
		static int VectorContains( const vector<string>&, const string& );
		static TString CondenseStrings( const vector<string>&, const int, const int );

		static vector<TString> GetStringContaining( const vector<TString>, const TString );
		static vector<TString> StripStrings( const vector<TString>, const TString );
		static vector<string> Convert( const vector<TString> );

		//	Get the left most character in a string in Numerical Form
		static int GetNumberOnLeft( const string& );

		//	Add a number to the start of a string
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

		/*!
		 * Return the full Path of a file if it exists:
		 * The priority of the file existing is:
		 *		local			i.e. ./filename.ext
		 *		RAPIDFITROOT=$PWD	i.e. $PWD/pdfs/configdata/filename.ext
		 *		RAPDFITROOT is defined	i.e. $RAPIDFITROOT/pdfs/configdata/filename.ext
		 */
		static string FindFileName( const string fileName );
};

#endif

