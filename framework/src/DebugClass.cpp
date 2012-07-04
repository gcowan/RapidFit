
#include "DebugClass.h"
#include "StringProcessing.h"

#include <string>
#include <vector>

using namespace::std;

DebugClass::DebugClass( bool input ) : perform_debugging( input ), classes_to_debug()
{
	if( !input ) classes_to_debug.push_back( "default" );
}

DebugClass::DebugClass( const DebugClass& input ) : perform_debugging( true ), classes_to_debug( input.classes_to_debug )
{
}

bool DebugClass::GetStatus()
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

vector<string> DebugClass::GetClassNames()
{
	return classes_to_debug;
}

bool DebugClass::DebugThisClass( string name )
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
			int num = StringProcessing::VectorContains( classes_to_debug, name );
			if( num == -1 ) return false;
			else return true;
		}
	}
}

