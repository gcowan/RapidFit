
#include "ComponentRef.h"

#include <string>

using namespace::std;

ComponentRef::ComponentRef( string input ) : thisSubComponent(NULL), thisIndex(-1), thisName( input )
{
}


ComponentRef::~ComponentRef()
{
	if( thisSubComponent != NULL ) delete thisSubComponent;
}

ComponentRef::ComponentRef( const ComponentRef& input ) :
	thisSubComponent( (input.thisSubComponent)==NULL?NULL:(new ComponentRef( *(input.thisSubComponent) )) ), thisIndex( input.thisIndex ), thisName( input.thisName )
{
}

void ComponentRef::addSubComponent( string input )
{
	if( thisSubComponent != NULL ) delete thisSubComponent;
	thisSubComponent = new ComponentRef( input );
}

ComponentRef* ComponentRef::getSubComponent()
{
	return thisSubComponent;
}

void ComponentRef::setComponentNumber( int input )
{
	thisIndex = input;
}

int ComponentRef::getComponentNumber()
{
	return thisIndex;
}

string ComponentRef::getComponentName()
{
	return thisName;
}

