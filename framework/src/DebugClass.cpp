
#include "DebugClass.h"
#include "StringProcessing.h"

#include <string>
#include <vector>
#include <iostream>

using namespace::std;

DebugClass::DebugClass( const bool input ) : perform_debugging( input ), classes_to_debug()
{
}

DebugClass::DebugClass( const DebugClass& input ) :
	perform_debugging(false), classes_to_debug()
{
	if( &input != NULL )
	{
		perform_debugging = input.perform_debugging;
		classes_to_debug = input.classes_to_debug;
	}
}

void DebugClass::SetDebugAll( const bool input )
{
	perform_debugging = input;
}

void DebugClass::SetClassNames( const vector<string> input )
{
	classes_to_debug = input;
}

vector<string> DebugClass::GetClassNames() const
{
	return classes_to_debug;
}

bool DebugClass::DebugThisClass( const string name ) const
{
	if( name.empty() ) return false;
	if( perform_debugging )
	{
		return true;
	}
	else
	{
		if( classes_to_debug.empty() )
		{
			return false;
		}
		else
		{
			string thisName=name;
			//cout << thisName << " : " << classes_to_debug[0] << endl;
			int num = -1;
			num = StringProcessing::VectorContains( &classes_to_debug, &thisName );
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
	return;
}

