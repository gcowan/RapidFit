
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TObject.h"

#include "DebugClass.h"
#include "StringProcessing.h"

#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <complex>

using namespace::std;

bool DebugClass::DebugAllStatus = false;

vector<string> DebugClass::classes_to_debug = vector<string>();

void DebugClass::SetDebugAll( const bool input )
{
	DebugAllStatus = input;
}

void DebugClass::SetClassNames( const vector<string> input )
{
	DebugClass::classes_to_debug = input;
}

vector<string> DebugClass::GetClassNames()
{
	return DebugClass::classes_to_debug;
}

bool DebugClass::DebugThisClass( const string& name )
{
	if( name.empty() ) return false;
	if( DebugClass::DebugAllStatus )
	{
		return true;
	}
	else
	{
		if( DebugClass::classes_to_debug.empty() )
		{
			return false;
		}
		else
		{
			string thisName=name;
			//cout << thisName << " : " << classes_to_debug[0] << endl;
			int num = -1;
			num = StringProcessing::VectorContains( &DebugClass::classes_to_debug, &thisName );
			//cout << num << endl;
			if( num == -1 ) return false;
			else return true;
		}
	}
	return false;
}

void DebugClass::SegFault()
{
	int *p = NULL;
	*p = 1;
	throw "error";
}

template<class T> void DebugClass::Dump2File( const string fileName, const vector<T> objects )
{
	stringstream fileContent;

	for( unsigned int i=0; i< objects.size(); ++i )
	{
		fileContent << objects[i] << endl;
	}

	ofstream outputFile;
	outputFile.open( fileName.c_str() );

	outputFile << fileContent.str();

	outputFile.close();
}

// explicit instantiations
template void DebugClass::Dump2File( const string, const vector<int> );
template void DebugClass::Dump2File( const string, const vector<char> );
template void DebugClass::Dump2File( const string, const vector<string> );
template void DebugClass::Dump2File( const string, const vector<double> );
template void DebugClass::Dump2File( const string, const vector<unsigned int> );
template void DebugClass::Dump2File( const string, const vector<complex<int> > );
template void DebugClass::Dump2File( const string, const vector<complex<double> > );

template<class T> void DebugClass::AppendToFile( const string fileName, const vector<T> objects )
{
	stringstream fileAppend;

	for( unsigned int i=0; i< objects.size(); ++i )
	{
		fileAppend << objects[i] << "\t";
	}

	fileAppend << endl;

	ofstream outputFile;
	outputFile.open( fileName.c_str(), std::fstream::app );

	outputFile << fileAppend.str();

	outputFile.close();
}

template void DebugClass::AppendToFile( const string fileName, const vector<int> objects );
template void DebugClass::AppendToFile( const string fileName, const vector<char> objects );
template void DebugClass::AppendToFile( const string fileName, const vector<string> objects );
template void DebugClass::AppendToFile( const string fileName, const vector<double> objects );
template void DebugClass::AppendToFile( const string fileName, const vector<unsigned int> objects );
template void DebugClass::AppendToFile( const string fileName, const vector<complex<int> > objects );
template void DebugClass::AppendToFile( const string fileName, const vector<complex<double> > objects );

template<class T> void DebugClass::Dump2TTree( const string fileName, const vector<T> objects, const string ttreeName, const string branchName )
{
	TFile* thisFile = new TFile( fileName.c_str(), "UPDATE" );

	TString TTreeName;
	if( ttreeName.empty() ) TTreeName = "someTTree";
	else TTreeName = ttreeName;

	TTree* thisTTree = NULL;
	
	thisTTree = (TTree*) thisFile->Get( TTreeName );
	
	if( thisTTree == NULL ) new TTree( TTreeName, TTreeName );

	T someObject;

	TString output_ext = extension<T>();

	TString BranchName;

	if( branchName.empty() ) BranchName = "someData";
	else BranchName = branchName;

	TBranch* thisBranch = thisTTree->Branch( BranchName, &someObject, BranchName+output_ext );

	thisTTree->SetEntries( (unsigned int)objects.size() );

	for( unsigned int i=0; i< objects.size(); ++i )
	{
		someObject = objects[i];
		thisBranch->Fill();
	}

	thisTTree->Write("", TObject::kOverwrite);

	thisFile->Write("", TObject::kOverwrite);
	thisFile->Close();
}

template void DebugClass::Dump2TTree( const string, const vector<int>, string, string );
template void DebugClass::Dump2TTree( const string, const vector<double>, string, string );
template void DebugClass::Dump2TTree( const string, const vector<float>, string, string );

template<class T> void DebugClass::Dump2File( const string fileName, const vector< vector<T> > objects )
{
	stringstream fileContent;

	for( unsigned int i=0; i< objects.size(); ++i )
	{
		for( unsigned int j=0; j< objects[i].size(); ++j )
		{
			fileContent << objects[i][j] << "   ";
		}
		fileContent << endl;
	}

	ofstream outputFile;
	outputFile.open( fileName.c_str() );

	outputFile << fileContent.str();

	outputFile.close();
}

template void DebugClass::Dump2File( const string, const vector<vector<int> > );
template void DebugClass::Dump2File( const string, const vector<vector<char> > );
template void DebugClass::Dump2File( const string, const vector<vector<string> > );
template void DebugClass::Dump2File( const string, const vector<vector<double> > );
template void DebugClass::Dump2File( const string, const vector<vector<unsigned int> > );
template void DebugClass::Dump2File( const string, const vector<vector<complex<int> > > );
template void DebugClass::Dump2File( const string, const vector<vector<complex<double> > > );

template<class T> void DebugClass::Dump2TTree( const string fileName, const vector< vector<T> > objects, const string ttreeName, const vector<string> branchName )
{
	TFile* thisFile = new TFile( fileName.c_str(), "APPEND" );

	TString TTreeName;
	if( ttreeName.empty() ) TTreeName = "someTTree";
	else TTreeName = ttreeName;

	TTree* thisTTree = NULL;
	
	thisTTree = (TTree*)thisFile->Get( TTreeName );
	
	if( thisTTree == NULL ) thisTTree = new TTree( TTreeName, TTreeName );

	T someObject;

	TString output_ext = extension<T>();

	for( unsigned int i=0; i< objects.size(); ++i )
	{
		TString BranchName;

		if( branchName.size() != objects.size() )
		{
			BranchName = "someData_"; BranchName += i;
		}
		else
		{
			if( branchName[i].empty() )	{	BranchName = "someData_"; BranchName+= i;	}
			else BranchName = branchName[i];
		}

		TBranch* thisBranch = thisTTree->Branch( BranchName, &someObject, BranchName+output_ext );

		thisTTree->SetEntries( (unsigned int)objects.size() );

		for( unsigned int j=0; j< objects[i].size(); ++j )
		{
			someObject = objects[i][j];
			thisBranch->Fill();
		}
	}

	thisTTree->Write("", TObject::kOverwrite);

	thisFile->Write("", TObject::kOverwrite);
	thisFile->Close();
}

template void DebugClass::Dump2TTree( const string, const vector<vector<int> >, string, vector<string> );
template void DebugClass::Dump2TTree( const string, const vector<vector<double> >, string, vector<string> );
template void DebugClass::Dump2TTree( const string, const vector<vector<float> >, string, vector<string> );

unsigned int DebugClass::GetGlobalCounter()
{
	return DebugClass::GlobalCounter;
}

void DebugClass::IncrementGlobalCounter()
{
	++DebugClass::GlobalCounter;
}

unsigned int DebugClass::GlobalCounter = 0;

map<string, string> DebugClass::extensions = map<string,string>();

template<typename T> TString DebugClass::extension()
{
	if (extensions.empty()) {
		extensions.insert(make_pair("f", "/F"));
		extensions.insert(make_pair("d", "/D"));
		extensions.insert(make_pair("i", "/I"));
		extensions.insert(make_pair("j", "/i"));
		extensions.insert(make_pair("y", "/l"));
		extensions.insert(make_pair("b", "/O"));
	}

	map<string, string>::const_iterator it = extensions.find(typeid(T).name());
	if (it != extensions.end()) {
		return TString(it->second.c_str());
	} else {
		cout << "Warning, cannot find extension for type " << typeid(T).name() << endl;
		return "";
	}
}

TString DebugClass::TGetUniqueROOTFileName()
{
	TString FileName = "RapidFitDebugClass_";
	FileName += DebugClass::GetGlobalCounter();
	DebugClass::IncrementGlobalCounter();
	FileName.Append( ".root" );
	return FileName;
}

string DebugClass::GetUniqueROOTFileName()
{
	return string( DebugClass::TGetUniqueROOTFileName().Data() );
}

TString DebugClass::TGetUniqueFileName()
{
	TString FileName = "RapidFitDebugClass_";
	FileName += DebugClass::GetGlobalCounter();
	DebugClass::IncrementGlobalCounter();
	FileName.Append( ".txt" );
	return FileName;
}

string DebugClass::GetUniqueFileName()
{
	return string( DebugClass::TGetUniqueFileName().Data() );
}

