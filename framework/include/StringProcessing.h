/**
        @class StringProcessing

        Collection of static functions for string processing

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef STRING_PROCESSING_H
#define STRING_PROCESSING_H

//	ROOT Headers
#include "TString.h"
//	System Headers
#include <string>
#include <vector>

using namespace std;

class StringProcessing
{
	public:
		static vector<string> SplitString( string, char );
		static int CharacterPosition( string, char );
		static vector<int> StringPositions( string, string );
		static void RemoveCharacter( string&, char );
		static string ReplaceString( string&, string, string );
		static void RemoveWhiteSpace( vector<string>& );
		static vector<string> CombineUniques( vector<string>, vector<string> );
		static int VectorContains( vector<string> const*, string const* );
		static TString CondenseStrings( vector<string>, int, int );
};

#endif
