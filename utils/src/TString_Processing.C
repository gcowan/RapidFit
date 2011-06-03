//	ROOT Headers
#include "TString.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TDirectory.h"
//	RapidFit Utils Header
#include "TString_Processing.h"
//	RapisFit Headers
#include "StringProcessing.h"
//	System Headers
#include <string>
#include <vector>
#include <iostream>

using namespace::std;

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

