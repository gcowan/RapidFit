
#pragma once
#ifndef _DEBUG_RAPIDFIT_HEADER
#define _DEBUG_RAPIDFIT_HEADER

#include <string>
#include <vector>

using namespace::std;

class DebugClass
{
	public:
		DebugClass( bool=true );
		DebugClass( const DebugClass& );
		~DebugClass(){};
		
		bool GetStatus();
		void SetStatus( bool input ); 

		void SetClassNames( vector<string> input );

		vector<string> GetClassNames();

		bool DebugThisClass( string name );

	private:

		bool perform_debugging;
		vector<string> classes_to_debug;

};

#endif

