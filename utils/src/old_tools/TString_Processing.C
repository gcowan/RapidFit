//	ROOT Headers
#include "TString.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TDirectory.h"
//	RapidFit Utils Header
#include "TString_Processing.h"
//	System Headers
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace::std;

vector<TString> postpend( vector<TString> input, TString suffix )
{
	vector<TString> output;
	for( vector<TString>::iterator i=input.begin(); i!= input.end(); ++i )
	{
		TString full_str = *i; full_str.Append( suffix );
		output.push_back( full_str );
	}
	return output;
}

//	Is a TString empty ?
bool is_empty( TString input )
{
	string temp("");
	if ( temp.compare( string(input.Data()) ) == 0 ) return true;
	return false;
}

//	Compare a string to a TString
int comp_TString( string test, TString input )
{
	return test.compare(input.Data());
}

//	Compare a TString to a string
int TString_comp( TString test, string input )
{
	return string(test.Data()).compare(input);
}

//  Useful for appendding a filename
void filename_append( TString input, TString *output, TString ext )
{
	size_t found = string::npos;
	string temp(input.Data());
	found=temp.find_last_of(".");
	if (found!=string::npos)      // Only designed for appending filenames...
	{
		temp.insert( found, ext );
	}
	output->Append(temp);
}

//	Wrapper
TString filename_append( TString input, TString ext )
{
	TString new_string;
	filename_append( input, &new_string, ext );
	return new_string;
}

string filename_append( string input, string ext )
{
	return string( filename_append( TString( input ), TString( ext ) ).Data() );
}

//	Originally Written by Conor Fitzpatrick
inline TString prettyPrint(Double_t value){
	char pretty[20];
	TString prettyString;
	sprintf (pretty, "%1.3g",value);
	prettyString = pretty;
	return prettyString;
}

//	Originally Written by Conor Fitzpatrick
inline TString prettyPrint(Double_t val, Double_t err){
	TString outstr = "$";
	char errstr[20];
	char valstr[20];
	Int_t n = sprintf (errstr, "%1.1g",err);
	sprintf(valstr,"%1.*f",n-2,val);
	outstr += valstr;
	outstr += "\\pm";
	outstr += errstr;
	outstr+= "$";
	return outstr;
}

//	Originally Written by Conor Fitzpatrick
TPaveText* addLHCbLabel(TString footer, bool DATA){
	//                              
	TPaveText * label = new TPaveText(0.18, 0.73, 0.18, 0.88,"BRNDC");
	label->SetFillStyle(0);         //Transparent i.e. Opacity of 0 :D
	label->SetBorderSize(0);
	label->SetTextAlign(11);
	label->SetTextSize(Float_t(0.04));
	TText * labeltext = 0;
	TString labeltstring( "LHC#font[12]{b} 2011" );
	if( DATA ) labeltstring.Append( " Data" );
	if( !DATA ) labeltstring.Append( " Simulation" );
	labeltext = label->AddText( labeltstring );
	labeltext = label->AddText("#sqrt{s} = 7TeV");
	labeltext = label->AddText(footer);
	(void) labeltext;
	return label;
}


//      Pass this an array of strings and it will find all strings matching a substring
vector<TString> filter_names( vector<TString> all_names, string substring )
{
	//      To be populated and returned to the user
	vector<TString> returnable_names;

	//      Address of the end of the string object
	size_t found=string::npos;

	//      Loop over all strings in the vector
	for( unsigned short int i=0; i<all_names.size(); ++i )
	{       //      Again using STL functions :)
		string temp_str = (all_names[i].Data());
		//      Attempt to find the coordinate of the substring within the string
		found = temp_str.find( substring );
		if( found!=string::npos )       //      If the substring is found
		{
			returnable_names.push_back(all_names[i]);
		}
	}
	//      Return all strings found containing the substring
	return returnable_names;
}

//      Pass this an array of strings and it will find all strings matching a substring
vector<TString> filter_names( vector<TString> all_names, TString substring )
{
	return filter_names( all_names, string(substring.Data()) );
}

//	WIP
//      Remove the extention from an array of strings
void strip_strings( vector<string>* list, string ext )
{
	vector<string> new_list;
	(void) list; (void) ext; (void) new_list;
}

vector<string> strip_all_strings( vector<string>* list, string ext )
{
	vector<string> new_list;

	for( vector<string>::iterator list_i = list->begin(); list_i != list->end(); ++list_i )
	{
		size_t found = string::npos;
		found=list_i->find_last_of( ext );

		//	If this extention was in the string, redefine the string to be everything upto the extention
		if ( found!=string::npos )
		{
			char output_obj;
			char* output=&output_obj;
			list_i->copy( output, (list_i->size()-found), 0 );
			new_list.push_back( string(output) );
		}
		else{
			new_list.push_back( *list_i );
		}
	}

	return new_list;
}

TString remove_suffix( TString whole_string, TString deliminator )
{
	//	To be done
	(void) whole_string; (void) deliminator;
	return TString();
}

string remove_suffix( string whole_string, string deliminator )
{
	//	To be done
	(void) whole_string; (void) deliminator;
	return string();
}

//Split a string every time you find a given character
vector<string> SplitString( string Input, char SplitCharacter )
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
int CharacterPosition( string Input, char SearchCharacter )
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
vector<int> StringPositions( string Input, string SearchString )
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
void RemoveCharacter( string & Input, char SearchCharacter )
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
string ReplaceString( string & Input, string FindString, string ReplaceString )
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
void RemoveWhiteSpace( vector<string> & newContent )
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
vector<string> CombineUniques( vector<string> VectorOne, vector<string> VectorTwo )
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

int VectorContains( vector<TString> const& InputVector, TString const& SearchString )
{
	vector<string> temp_vec=TString2string(InputVector);
	string temp_str( SearchString.Data() );
	return VectorContains( &(temp_vec), &temp_str );
}

//Return the position of a search string within a vector of strings, or -1 if not found
int VectorContains( vector<string> const* InputVector, string const* SearchString )
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
TString CondenseStrings( vector<string> input_vec, int lolim=-1, int hilim=-1 )
{
	int temp_int = int( input_vec.size() );
	if( ( hilim > temp_int ) || (hilim < 0) ) hilim = temp_int;
	if( lolim < 0 ) lolim = 0;

	TString Returnable_String( "" );
	for( int iter = lolim; iter < hilim ; ++iter )
	{
		Returnable_String.Append( input_vec[(unsigned)iter] );
	}
	return Returnable_String;
}

vector<TString> GetStringContaining( vector<TString> list, TString search_str )
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

vector<TString> StripStrings( vector<TString> list, TString ext )
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

vector<string> TString2string( vector<TString> input )
{
	vector<string> output;
	for( vector<TString>::iterator i=input.begin(); i!=input.end(); ++i )   output.push_back( string( i->Data() ) );
	return output;
}

vector<TString> string2TString( vector<string> input )
{
	vector<TString> output;
	for( vector<string>::iterator i=input.begin(); i!=input.end(); ++i )    output.push_back( string( i->c_str() ) );
	return output;
}

TString RemoveSuffix( TString input, TString suffix )
{
	size_t found;
	string suffix_str( suffix.Data() );
	string input_str( input.Data() );
	found = input_str.find( suffix_str );
	string output_str;
	for( unsigned int i=0; i<found; ++i )
	{
		output_str.push_back( input_str[i] );
	}
	TString output( output_str.c_str() );
	return output;
}

TString getFileName( TString input )
{
	string input_str( input.Data() );
	size_t found=input_str.find_last_of("/\\");
	string output_str;
	for( size_t i=found+1; i<input_str.size(); ++i )
	{
		output_str.push_back( input_str[i] );
	}
	return TString( output_str.c_str() );
}

TString Clean( TString input )
{
	string temp(input.Data());

	replace(temp.begin(), temp.end(), '.', '_');
	replace(temp.begin(), temp.end(), '/', '_');
	replace(temp.begin(), temp.end(), ' ', '_' );

	vector<string::iterator> remove_list;
	string::iterator ch=temp.begin();
	for( unsigned int i=0; i< temp.size(); ++i, ++ch )
	{
		if( i+1 != temp.size() )
		{
			if( temp[i] == temp[i+1] )
			{	
				if( temp[i] == '_' )	remove_list.push_back( ch );
			}
		}
	}
	if( temp[0] == '_' ) remove_list.push_back( temp.begin() );

	for( unsigned int i=0; i< remove_list.size(); ++i )
	{
		temp.erase( remove_list[i] );
	}

	return TString( temp.c_str() );
}

