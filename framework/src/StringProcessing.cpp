/**
  @class StringProcessing

  Collection of static functions for string processing

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	ROOT Headers
#include "TString.h"
//	RapidFit Headers
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <algorithm>

//Split a string every time you find a given character
vector<string> StringProcessing::SplitString( string Input, char SplitCharacter )
{
	vector<string> splitParts;

	while(true)
	{
		//Search for the character
		int position = CharacterPosition( Input, SplitCharacter );

		if ( position == -1 )
		{
			//Ignore empty strings
			if ( Input != "" )
			{
				splitParts.push_back( Input );
			}

			//If it's not found, you've reached the end
			break;
		}
		else
		{
			//Split the string at the character postion
			string tempString( Input, 0, unsigned(position) );
			string newInput( Input, unsigned(position + 1) );

			//Ignore empty strings
			if ( tempString != "" )
			{
				splitParts.push_back( tempString );
			}
			Input = newInput;
		}
	}

	return splitParts;
}

//Return the position of the first instance of a character in a string
int StringProcessing::CharacterPosition( string Input, char SearchCharacter )
{
	for (unsigned int characterIndex = 0; characterIndex < Input.size(); ++characterIndex)
	{
		//Look for the character
		if ( Input[characterIndex] == SearchCharacter )
		{
			return int(characterIndex);
		}

		//If you get to the end, character not found
		if ( characterIndex == Input.size() - 1 )
		{
			return -1;
		}
	}
	return -1;
}

//Return the position of all instances of a string in another string
vector<int> StringProcessing::StringPositions( string Input, string SearchString )
{
	vector<int> positions;
	int numberDiscarded = 0;

	while (true)
	{
		int firstCharacterPosition = CharacterPosition( Input, SearchString[0] );
		if ( firstCharacterPosition == -1 )
		{
			//If you can't find the first character of the string, you won't find the string!
			return positions;
		}
		else
		{
			bool found = true;
			for (unsigned int characterIndex = 0; characterIndex < SearchString.size(); ++characterIndex )
			{
				//Compare each subsequent character
				if ( Input[ characterIndex + unsigned(firstCharacterPosition) ] != SearchString[characterIndex] )
				{
					found = false;
					break;
				}
			}

			if (found)
			{
				positions.push_back( firstCharacterPosition + numberDiscarded );
			}

			//Search the rest of the string
			string tempString( Input, unsigned(firstCharacterPosition + 1) );
			numberDiscarded += unsigned(firstCharacterPosition + 1);
			Input = tempString;
		}
	}

	return positions;
}

//Remove any instances of a particular character in a string
void StringProcessing::RemoveCharacter( string & Input, char SearchCharacter )
{
	char* passedCharacters = new char[ Input.size() + 1 ];
	int addedCharacters = 0;

	for (unsigned int characterIndex = 0; characterIndex < Input.size(); ++characterIndex )
	{
		if ( Input[characterIndex] != SearchCharacter )
		{
			passedCharacters[addedCharacters] = Input[characterIndex];
			++addedCharacters;
		}
	}

	passedCharacters[addedCharacters] = '\0';
	Input = passedCharacters;
}

//Replace any instances of a particular character in a string
string StringProcessing::ReplaceString( string & Input, string FindString, string ReplaceString )
{
	string output;

	//Search the input character by character
	for (unsigned int characterIndex = 0; characterIndex < Input.size(); ++characterIndex )
	{
		//Find out if the full search string is present
		bool isInstance = true;
		for (unsigned int testIndex = 0; testIndex < FindString.size(); ++testIndex )
		{
			if ( Input[ characterIndex + testIndex ] != FindString[testIndex] )
			{
				isInstance = false;
				break;
			}
		}

		//Replace or push back
		if (isInstance)
		{
			for (unsigned int replaceIndex = 0; replaceIndex < ReplaceString.size(); ++replaceIndex )
			{
				output.push_back( ReplaceString[replaceIndex] );
			}

			characterIndex += unsigned(FindString.size() - 1);
		}
		else
		{
			output.push_back( Input[characterIndex] );
		}
	}

	return output;
}

//Remove white space from passed lines
void StringProcessing::RemoveWhiteSpace( vector<string> & newContent )
{
	vector<string> output;

	for (unsigned int lineIndex = 0; lineIndex < newContent.size(); ++lineIndex )
	{
		//Remove tabs
		RemoveCharacter( newContent[lineIndex], '\t' );

		//Remove empty lines
		if ( newContent[lineIndex] != "" )
		{
			output.push_back( newContent[lineIndex] );
		}
	}

	newContent = output;
}

//Return a vector containing all the unique strings from the two input vectors
vector<string> StringProcessing::CombineUniques( vector<string> VectorOne, vector<string> VectorTwo )
{
	vector<string> result;
	vector<string>::iterator stringIterator;

	//Don't assume VectorOne is unique
	for ( stringIterator = VectorOne.begin(); stringIterator != VectorOne.end(); ++stringIterator )
	{
		if ( VectorContains( &result, &(*stringIterator) ) == -1 )
		{
			result.push_back( *stringIterator );
		}
	}

	//Now add in VectorTwo
	for ( stringIterator = VectorTwo.begin(); stringIterator != VectorTwo.end(); ++stringIterator )
	{
		if ( VectorContains( &result, &(*stringIterator) ) == -1 )
		{
			result.push_back( *stringIterator );
		}
	}

	return result;
}

//Return the position of a search string within a vector of strings, or -1 if not found
int StringProcessing::VectorContains( vector<string> const* InputVector, string const* SearchString )
{
	vector<string>::const_iterator begin = InputVector->begin();
	vector<string>::const_iterator ending = InputVector->end();
	vector<string>::const_iterator result = find( begin, ending, *SearchString);
	int position=int( result - InputVector->begin() );
	if( position >= 0 && position < int( InputVector->size() ) ) return position;

	//If you've got this far, it wasn't found
	return -1;
}

//Return the a TString object which is composed of the selected strings within the vector
TString StringProcessing::CondenseStrings( vector<string> input_vec, int lolim=-1, int hilim=-1 )
{
	int temp_int = int( input_vec.size() );
	if( ( hilim > temp_int ) || (hilim < 0) ) hilim = temp_int;
	if( lolim < 0 ) lolim = 0;

	TString Returnable_String( "" );
	for( int iter = lolim; iter < hilim ; ++iter )
	{
		Returnable_String.Append( input_vec[iter] );
	}
	return Returnable_String;
}

vector<TString> StringProcessing::GetStringContaining( vector<TString> list, TString search_str )
{
	vector<TString> output_list;

	for( unsigned int i=0; i< list.size(); ++i )
	{
		string input = list[i].Data();
		size_t found = input.find( search_str.Data() );
		if( found != string::npos )
		{
			output_list.push_back( list[i] );
		}
	}

	return output_list;
}

vector<TString> StringProcessing::StripStrings( vector<TString> list, TString ext )
{
	vector<TString> output_list;
	string remove = ext.Data();
	for( unsigned int i=0; i< list.size(); ++i )
	{
		string object = list[i].Data();
		size_t found = object.find( remove );
		TString TObject_Str;
		if( found != string::npos )
		{
			TObject_Str = object.substr( 0, object.length()-remove.length() ).c_str();
		} else {
			TObject_Str = object.c_str();
		}
		output_list.push_back( TObject_Str );
	}
	return output_list;
}

vector<string> StringProcessing::Convert( vector<TString> input )
{
	vector<string> output;
	for( unsigned int i=0; i< input.size(); ++i )
	{
		output.push_back( input[i].Data() );
	}
	return output;
}

