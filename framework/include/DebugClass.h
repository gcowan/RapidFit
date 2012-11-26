
#pragma once
#ifndef _DEBUG_RAPIDFIT_HEADER
#define _DEBUG_RAPIDFIT_HEADER

#include <string>
#include <vector>

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
	private:

		bool perform_debugging;
		vector<string> classes_to_debug;

};

#endif

