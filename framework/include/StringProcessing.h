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
		/*!
		 * @brief Return the current time as a time-stamp flat string
		 *
		 * @return Returns the current time (down to seconds) as a flat string
		 */
		static string TimeString();

		/*!
		 * @brief Split a string into a vector of strings by the deliminator
		 *
		 * @Param input this is the input string to be split
		 *
		 * @Param delim this is the deliminiator char to split the string by
		 *
		 * @return This Returns a vector of sub-strings from the input string split by the delim, a vector of 1 string if the delim is not found
		 */
		static vector<string> SplitString( const string input, const char delim );

		/*!
		 * @brief Get the Position of a Character in a string (first occurance)
		 *
		 * @Param input this is the string to search for the character of interest
		 *
		 * @Param find this is the char that is searched for in the string
		 *
		 * @return Returns the absolute position of the find character in the input string, or -1 if it isn't found
		 */
		static int CharacterPosition( const string input, const char find );

		/*!
		 * @brief Get all of the string Positions of a character in a string
		 *
		 * @Param input1	This is the String thought to contain input2
		 *
		 * @Param input2	This is the string to be searched for in input1
		 *
		 * @return returns a vector of the positions the substring input2 appears in input1
		 */
		static vector<int> StringPositions( const string input1, const string input2 );

		/*!
		 * @brief Remove the Leading character from a string specified by the user
		 *
		 * @param Input			This is the string which is to have it's leading character removed if it matches SearchCharacter
		 *
		 * @param SearchCharacter	This is the character which is to be removed
		 *
		 * @return void
		 */
		static void RemoveLeadingCharacters( string & Input, const char SearchCharacter );

		/*!
		 * @brief Remove the first occurance? of a selected character
		 *
		 * @param inputStr		String to have characters removed from
		 *
		 * @param inputChar		Character that is to be removed in the inputStr
		 *
		 * @return void
		 */
		static void RemoveCharacter( string& inputStr, const char inputChar );

		/*!
		 * @brief Replace an instance of a substring within a 
		 *
		 * @param inputStr	This is the string which should have it's substrings replaced
		 *
		 * @param FindStr	This is the substring which should be replaced
		 *
		 * @param ReplaceStr	This is the string to replace the FindStr in inputStr with
		 *
		 * @return Returns a new string which has had every occurance of FindStr replaced with ReplaceStr
		 */
		static string ReplaceString( const string& inputStr, const string FindStr, const string ReplaceStr );

		/*!
		 * @brief Remove Leading Whitespace characters from a given vector of strings
		 *
		 * @param input		This is the vector of strings that any leading whitespace is to be removed from
		 *
		 * @return void
		 */
		static void RemoveLeadingWhiteSpace( vector<string>& input );

		/*!
		 * @brief Remove leading whitespace tab characters from an input vector of strings
		 *
		 * @param input		This is the vector of strings which will be modified to have the leading tabs removed
		 *
		 * @return void
		 */
		static void RemoveWhiteSpace( vector<string>& input );

		/*!
		 * @brief Replace non-Latex safe characters with '_'
		 *
		 * @param input		This is the string which is intended to be made 'latex safe;
		 *
		 * @return returns a 'Latex safe' copy of the input string
		 */
		static string LatexSafe( const string input );

		/*!
		 * @brief Replace non-Latex safe chatacters with '_'
		 *
		 * @param input		This is the TString which is intended to be made 'latex safe'
		 *
		 * @return returns a 'Latex safe' copy of the input string
		 */
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

		/*!
		 * @brief This will return a new vector of output strings constructed from the unique strings within the inputVec
		 *
		 * @param inputVec		This is the input vector of strings which may not be unique
		 *
		 * @param duplicated		This is the vector of strings which occur multiple times in the inputVec
		 *
		 * @return This is is a new vector of strings assembled from the unique entries in the inputVec 
		 */
		static vector<string> RemoveDuplicates( const vector<string> inputVec, vector<string>& duplicated = *((vector<string>*) NULL) );

		/*!
		 * @brief This routine performs a lookup for the stringRef object within the vectorRef and returns the location of the found string (-1 if not found)
		 *
		 * @param vectorRef		This is a reference of a vector which is to be searched for the stringRef
		 *
		 * @param stringRef		This is a reference to the string that is to be searched for within the vectorRef
		 *
		 * @return Returns an integer containing the location of the stringRef within the vectorRef or -1 if not found
		 */
		static int VectorContains( vector<string> const* vectorRef, string const* stringRef );

		/*!
		 * @brief This routine performs a lookup for the inputVec object within the inputStr and returns the location of the found string (-1 if not found)
		 *
		 * @param inputVec		This is the vector os strings which is to be searched for the inputStr
		 *
		 * @param inputStr		This is the string which is to be looked for within the inputVec
		 *
		 * @return Returns an integer containing the location of the inputStr within the inputVec or -1 if not found
		 */
		static int VectorContains( const vector<string>& inputVec, const string& inputStr );

		/*!
		 * @brief Return the a TString object which is composed of the selected strings within the vector
		 *
		 * @param stringVec		This is a vector of strings which are intended to be made into a single TString
		 *
		 * @param lower			This is the lowest index to start within the stringVec
		 *
		 * @param uppper		This is the highest index to start within the stringVec
		 *
		 * @return returns a new TString object which is composed of elements from within the inputVector within limits
		 */
		static TString CondenseStrings( const vector<string>& stringVec, const int lower, const int upper );

		/*!
		 * @brief return the TStrings which contain the search string as a substring
		 *
		 * @param inputVec		This is the vector of strings to be searched for the inputString
		 *
		 * @param inputString		This is the string which is to be searched for within the entries of the inputVec
		 *
		 * @return Returns a new vector of TStrings which each contain the inputString object as a substring
		 */
		static vector<TString> GetStringContaining( const vector<TString> inputVec, const TString inputString );

		/*!
		 * @brief This returns a new vector of TStrings which have had the chosen ext removed if found
		 *
		 * @param list		This is the list of strings which should have the ext removed
		 *
		 * @param ext		This is the ext which should be removed from the input elements of the list if found
		 *
		 * @return Returns a new vector of TStrings which have the ext removed from each element
		 */
		static vector<TString> StripStrings( const vector<TString> list , const TString ext );

		/*!
		 * @brief This converts a vector of TStrings to a vector of strings
		 *
		 * @param input		This is the vector of TStrings which should be converted into a vector of strings
		 *
		 * @return Returns a new vector of strings based on the input vector of TStrings
		 */
		static vector<string> Convert( const vector<TString> input );

		/*!
		 * @brief Get the left most character in a string in Numerical Form
		 *
		 * @param inputStr		This is the string which should have the leading character converted to an integer
		 *
		 * @return Returns the left most character converted to an integer using atoi
		 */
		static int GetNumberOnLeft( const string& inputStr );

		/*!
		 * @brief Add a number to the start of a string
		 *
		 * @param inputStr		This is the input string which should have a new character added to the start
		 *
		 * @param inputNum		This is the integer which should be converted to a character to be added to the start of the inputStr
		 *
		 * @return Returns a new string which has the inputNum converted to a character and the rest of the string matching the inputStr
		 */
		static string AddNumberToLeft( const string& inputStr, const int& inputNum );

		/*!
		 * @brief This returns a substring which corresponds to the input string minus the first character
		 *
		 * @param input			This is th input string which should have the first character removed
		 *
		 * @return returns a new string based on the input string minus the first character
		 */
		static string RemoveFirstNumber( const string& input );

		/*!
		 * @brief This creates a new list of strings which match the natural numbers between start and end
		 *
		 * @param start			This is the start of the series of numbers
		 *
		 * @param end			This is the end of the series of numbers
		 *
		 * @return Returns a vector of strings which are numbers between start and end
		 */
		static vector<string> FillList( const int start, const int end =0 );

		/*!
		 * @brief This returns a new string composed of "name1+name2" with empty strings replaced with unknown
		 *
		 * @param name1			This is the first input string
		 *
		 * @param name2			This is the second input string
		 *
		 * @return returns a new string which is composed of the input strings and "+"
		 */
		static string AddNames( const string& name1, const string& name2 );

		/*!
		 * @brief This returns a new string composed of "name1xname2" with empty strings replaced with unknown
		 *
		 * @param name1			This is the first input string
		 *
		 * @param name2			This is the second input string
		 *
		 * @return returns a new string which is composed of the input strings and "x"
		 */
		static string MultNames( const string& name1, const string& name2 );

		/*!
		 * @brief This removes 'unsafe' characters from a TString an replaces them with a '_'
		 *      ported from utils / StringOperations
		 *      Remove characters not nice in filenames and multiple '_' occurances
		 *
		 * @param input This is the FileName which may contain unsafe characters
		 *
		 * @return returns a new TString with unsafe characters removed from the filename
		 */
		static TString Clean( const TString input );

		/*!
		 * @brief Is a TString empty?
		 *
		 * @param input This is the TString that is being tested
		 *
		 * @return true if the TString is empty or just contains whitespace, false if not
		 */
		static bool is_empty( const TString input );

		/*!
		 * @brief Returns the Full path of a given file after looking in some 'sane' expected places for it
		 *
		 * Return the full Path of a file if it exists:
		 * The priority of the file existing is:
		 *		local			i.e. ./filename.ext
		 *		RAPIDFITROOT=$PWD	i.e. $PWD/pdfs/configdata/filename.ext
		 *		RAPDFITROOT is defined	i.e. $RAPIDFITROOT/pdfs/configdata/filename.ext
		 *
		 * @param fileName This is the fileName to start from to try and resolve to find a file on disk
		 *
		 * @return Returns the path to a file if found that allows the file to be opened with open or TFile constructor
		 */
		static string FindFileName( const string fileName, bool quiet=false );

		/*!
		 * @brief 
		 *
		 * Returns a vector populated with the contents of the original vector
		 * If the element is NOT in the input vector this returns the input
		 * else
		 * The vector that is returned has the element of choice moved to the start and the order of the elements otherwise is preserved
		 *
		 * @param input This is the vector of string which should be searched
		 *
		 * @param element This is the string element which should be moved to the [0]th element of the vector if it has been found
		 *
		 * @returns a vector of srings with the chosen element (if found) moved to the [0]th element of the vector
		 */
		static vector<string> MoveElementToStart( const vector<string> input, const string element );

	private:
		StringProcessing();
};

#endif

