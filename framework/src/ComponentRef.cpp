
#include "ComponentRef.h"

#include <string>

using namespace::std;

ComponentRef::ComponentRef( const string input, const string obs ) :
	thisSubComponent(NULL), thisIndex(-1), thisName( input ), obsName( obs )
{
}

ComponentRef::~ComponentRef()
{
	if( thisSubComponent != NULL ) delete thisSubComponent;
}

ComponentRef::ComponentRef( const ComponentRef& input ) :
	thisSubComponent( NULL ), thisIndex( input.thisIndex ),
	thisName( input.thisName ), obsName( input.obsName )
{
	if(  input.thisSubComponent != NULL )
	{
		thisSubComponent = new ComponentRef( *(input.thisSubComponent) );
	}
}

void ComponentRef::addSubComponent( const string input )
{
	//int stored_index = -1;
	//string stored_string = "";
	if( thisSubComponent != NULL )
	{
		//stored_index = thisSubComponent->getComponentNumber();
		//stored_string = thisSubComponent->getComponentName();
		delete thisSubComponent;
	}
	thisSubComponent = new ComponentRef( input, obsName );
	//if( input.compare( stored_string ) == 0 ) thisSubComponent->setComponentNumber( stored_index );
}

ComponentRef* ComponentRef::getSubComponent() const
{
	return thisSubComponent;
}

void ComponentRef::setComponentNumber( const int input )
{
	thisIndex = input;
}

int ComponentRef::getComponentNumber() const
{
	return thisIndex;
}

string ComponentRef::getComponentName() const
{
	return thisName;
}

string ComponentRef::getObservableName() const
{
	return obsName;
}

void ComponentRef::setObservableName( const string input )
{
	obsName = input;
}

