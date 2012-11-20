
#pragma once
#ifndef _DEBUG_RAPIDFIT_HEADER
#define _DEBUG_RAPIDFIT_HEADER

#include <string>
#include <vector>

using namespace::std;

class DebugClass
{
	public:
		DebugClass( bool=false );
		DebugClass( const DebugClass& );
		~DebugClass(){};

		void SetDebugAll( bool=true );

		void SetClassNames( vector<string> input );

		vector<string> GetClassNames() const;

		bool DebugThisClass( const string name ) const;

	private:

		bool perform_debugging;
		vector<string> classes_to_debug;

};

#endif

