
#include "DebugClass.h"
#include "StringProcessing.h"

#include <string>
#include <vector>
#include <iostream>

using namespace::std;

DebugClass::DebugClass( bool input ) : perform_debugging( input ), classes_to_debug()
{
	if( !input ) classes_to_debug.push_back( "default" );
}

DebugClass::DebugClass( const DebugClass& input ) :
	perform_debugging( input.perform_debugging ), classes_to_debug( input.classes_to_debug )
{
}

bool DebugClass::GetStatus() const
{
	return perform_debugging;
}

void DebugClass::SetStatus( bool input )
{
	perform_debugging = input;
}

void DebugClass::SetClassNames( vector<string> input )
{
	classes_to_debug = input;
}

vector<string> DebugClass::GetClassNames() const
{
	return classes_to_debug;
}

bool DebugClass::DebugThisClass( const string name )
{
	if( !perform_debugging ) return false;
	else
	{
		if( classes_to_debug.empty() )
		{
			return true;
		}
		else
		{
			string thisName=name;
			//cout << thisName << " : " << classes_to_debug[0] << endl;
			int num = StringProcessing::VectorContains( &classes_to_debug, &thisName );
			//cout << num << endl;
			if( num == -1 ) return false;
			else return true;
		}
	}
}

