
#pragma once
#ifndef _DEBUG_RAPIDFIT_HEADER
#define _DEBUG_RAPIDFIT_HEADER

#include "TString.h"

#include <string>
#include <vector>
#include <map>

using namespace::std;

class DebugClass
{
	public:
		DebugClass( const bool=false );
		DebugClass( const DebugClass& );

		void SetDebugAll( const bool=true );

		void SetClassNames( const vector<string> input );

		vector<string> GetClassNames() const;

		bool DebugThisClass( const string name ) const;

		static void SegFault();

		template<class T> static void Dump2File( const string fileName, const vector<T> objects );

		template<class T> static void Dump2TTree( const string fileName, const vector<T> objects, const string ttreeName="", const string branchName="" );

		template<class T> static void Dump2File( const string fileName, const vector< vector<T> > objects );

		template<class T> static void Dump2TTree( const string fileName, const vector< vector<T> > objects, const string ttreeName="", const vector<string> branchName=vector<string>() );

		static unsigned int GetGlobalCounter();

		static void IncrementGlobalCounter();

		template<typename T> static TString extension();

		static TString TGetUniqueROOTFileName();

		static string GetUniqueROOTFileName();

		static TString TGetUniqueFileName();

		static string GetUniqueFileName();

		template<class T> static void AppendToFile( const string fileName, const vector<T> objects );

	private:

		bool perform_debugging;
		vector<string> classes_to_debug;

		static unsigned int GlobalCounter;

		static map<string, string> extensions;
};

#endif

