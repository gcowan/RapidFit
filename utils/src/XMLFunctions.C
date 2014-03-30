
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"

#include "StringOperations.h"
#include "ROOT_File_Processing.h"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace::std;

void RestoreXML( vector<string> input_filenames, vector<string> other_params )
{
	string SaveAs="--SaveAs";
	string SaveAs2="-SaveAs";

	TString XMLName = "XMLFile_";

	for( unsigned int i=0; i< other_params.size(); ++i )
	{
		string thisString = other_params[i];
		size_t found = thisString.find(SaveAs);
		size_t found2 = thisString.find(SaveAs2);

		if( found2 != string::npos )
		{
			found = found2;
		}

		if( found != string::npos )
		{
			string newName = StringOperations::SplitString( thisString, ':' )[1];
			string xmlext=".xml";
			size_t foundext = newName.find(xmlext);
			if( foundext != string::npos ) newName = StringOperations::SplitString( thisString, '.' )[0];

			XMLName = newName.c_str();
			XMLName.Append("_");
		}
	}

	for( unsigned int i=0; i< input_filenames.size(); ++i )
	{
		TString thisXMLName = XMLName;

		thisXMLName+=i; thisXMLName.Append("_");

		thisXMLName.Append( StringOperations::TimeString() );

		thisXMLName.Append( ".xml" );

		cout << "Requesting FittingXML" << endl;

		TTree* runtimeXML = ROOT_File_Processing::GetTree( input_filenames[i], "FittingXML" );

		vector<string>* thisXML = new vector<string>();

		runtimeXML->SetBranchAddress( "FittingXML", &thisXML );

		runtimeXML->GetEntry(0);

		cout << "Saving output XML to file:\t" << thisXMLName << endl;
		cout << endl;

		stringstream full_xml;

		for( unsigned int j=0; j< thisXML->size(); ++j )
		{
			full_xml << (*thisXML)[j] << endl;
		}

		ofstream output_xmlFile;
		output_xmlFile.open( thisXMLName.Data() );

		output_xmlFile << full_xml.str();

		output_xmlFile.close();
	}
}

void GetToyXML( vector<string> input_filenames, vector<string> other_params )
{
	string SaveAs="--SaveAs";
	string SaveAs2="-SaveAs";

	TString XMLName = "ToyXMLFile_";

	for( unsigned int i=0; i< other_params.size(); ++i )
	{
		string thisString = other_params[i];
		size_t found = thisString.find(SaveAs);
		size_t found2 = thisString.find(SaveAs2);

		if( found2 != string::npos )
		{
			found = found2;
		}

		if( found != string::npos )
		{
			string newName = StringOperations::SplitString( thisString, ':' )[1];
			string xmlext=".xml";
			size_t foundext = newName.find(xmlext);
			if( foundext != string::npos ) newName = StringOperations::SplitString( thisString, '.' )[0];

			XMLName = newName.c_str();
			XMLName.Append("_");
		}
	}

	for( unsigned int i=0; i< input_filenames.size(); ++i )
	{
		TString thisXMLName = XMLName;

		thisXMLName+=i; thisXMLName.Append("_");

		thisXMLName.Append( StringOperations::TimeString() );

		thisXMLName.Append( ".xml" );

		TTree* runtimeXML = ROOT_File_Processing::GetTree( input_filenames[i], "XMLForToys" );

		string* thisXML = new string();

		runtimeXML->SetBranchAddress( "ToyXML", &thisXML );

		runtimeXML->GetEntry(0);

		cout << "Saving output XML to file:\t" << thisXMLName << endl;
		cout << endl;

		stringstream full_xml;

		full_xml << (*thisXML) << endl;

		ofstream output_xmlFile;
		output_xmlFile.open( thisXMLName.Data() );

		output_xmlFile << full_xml.str();

		output_xmlFile.close();
	}
}

void GetProjectionXML( vector<string> input_filenames, vector<string> other_params )
{
	string SaveAs="--SaveAs";
	string SaveAs2="-SaveAs";

	TString XMLName = "ProjectionXMLFile_";

	for( unsigned int i=0; i< other_params.size(); ++i )
	{
		string thisString = other_params[i];
		size_t found = thisString.find(SaveAs);
		size_t found2 = thisString.find(SaveAs2);

		if( found2 != string::npos )
		{
			found = found2;
		}

		if( found != string::npos )
		{
			string newName = StringOperations::SplitString( thisString, ':' )[1];
			string xmlext=".xml";
			size_t foundext = newName.find(xmlext);
			if( foundext != string::npos ) newName = StringOperations::SplitString( thisString, '.' )[0];

			XMLName = newName.c_str();
			XMLName.Append("_");
		}
	}

	for( unsigned int i=0; i< input_filenames.size(); ++i )
	{
		TString thisXMLName = XMLName;

		thisXMLName+=i; thisXMLName.Append("_");

		thisXMLName.Append( StringOperations::TimeString() );

		thisXMLName.Append( ".xml" );

		TTree* runtimeXML = ROOT_File_Processing::GetTree( input_filenames[i], "XMLForProjections" );

		string* thisXML = new string();

		runtimeXML->SetBranchAddress( "ProjectionXML", &thisXML );

		runtimeXML->GetEntry(0);

		cout << "Saving output XML to file:\t" << thisXMLName << endl;
		cout << endl;

		stringstream full_xml;

		full_xml << (*thisXML) << endl;

		ofstream output_xmlFile;
		output_xmlFile.open( thisXMLName.Data() );

		output_xmlFile << full_xml.str();

		output_xmlFile.close();
	}
}
