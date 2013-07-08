/**
  @class XMLTag

  A class used to parse an xml config file. Represents a single xml tag and its contents.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

#include "TF1.h"
//	RapidFit Headers
#include "XMLTag.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <stdlib.h>
#include <string>

using namespace::std;

//Default constructor w Override tags
XMLTag::XMLTag( vector<pair<string,string> >* Override_Tags ): children(), value(), name("RapidFit"), parent(NULL), path(""), forbidden(Override_Tags)
{
}

//Constructor with correct arguments
XMLTag::XMLTag( string TagName, vector<string> TagContent, XMLTag* Parent ) : children(), value(), name(TagName), parent(Parent), path(Parent->GetPath()), forbidden(Parent->GetForbidden())
{
	path.Append( "/" + name );
	//Find the child tags
	children = FindTagsInContent( TagContent, value );
}

//Destructor
XMLTag::~XMLTag()
{
	//Delete the children!
	vector< XMLTag* >::iterator childIterator;
	for ( childIterator = children.begin(); childIterator != children.end(); ++childIterator )
	{
		delete *childIterator;
	}
}

vector<pair<string,string> >* XMLTag::GetForbidden() const
{
	return forbidden;
}

//Return the complete Path of the Object within the XML Structure
string XMLTag::GetPath() const
{
	return path.Data();
}

//Return the tag name
string XMLTag::GetName() const
{
	return name;
}

//Return the children
vector< XMLTag* > XMLTag::GetChildren() const
{
	string NameStr="Name";
	vector<string> children_names;
	int NamePos=0;
	for( vector<XMLTag*>::const_iterator child_i = children.begin(); child_i != children.end(); ++child_i )
	{
		children_names.push_back( (*child_i)->GetName() );
	}
	NamePos = StringProcessing::VectorContains( &children_names, &NameStr );
	if( NamePos != -1  )
	{
		for( unsigned int i=0; i< children.size(); ++i )
		{
			children[i]->PrePendPath( "/" + children[(unsigned)NamePos]->GetValue()[0] );
		}
	}

	return children;
}

void XMLTag::RemoveChild( int num )
{
	if( num<0 || num>=(int)children.size() ) return;

	vector<XMLTag*>::iterator toBeRemoved;

	int counter=0;
	for( vector<XMLTag*>::iterator child_i = children.begin(); child_i != children.end(); ++child_i, ++counter )
	{
		if( counter == num )
		{
			toBeRemoved = child_i;

			break;
		}
	}

	if( children[(unsigned)num] != NULL ) delete children[(unsigned)num];
	children.erase( toBeRemoved );
}

vector<string> XMLTag::GetRAWValue() const
{
	return this->GetValue();
}

//Return the value
vector<string> XMLTag::GetValue() const
{
	vector<string> new_value = value;
	if( value.empty() ) return vector<string>();
	//cout << "MyPath " << path << endl;
	//path = parent->GetPath();
	//cout << "Parent Path " << path << endl;
	//path.Append( "/" + name );
	//cout << "ModPar Path " << path << endl;
	for( unsigned int i=0; i< children.size(); ++i )
	{
		children[i]->RegeneratePath();
	}
	//cout << path << endl;
	string pathStr = path.Data();
	vector<string> forbidden_paths;
	if( forbidden != NULL )
	{
		for( vector<pair<string,string> >::const_iterator fpath_i = forbidden->begin(); fpath_i != forbidden->end(); ++fpath_i )
		{
			forbidden_paths.push_back( fpath_i->first );
		}
		int TagPos = StringProcessing::VectorContains( &forbidden_paths, &pathStr );
		if( TagPos != -1 )
		{
			value[0] = (*forbidden)[(unsigned)TagPos].second;
			forbidden->erase(forbidden->begin()+TagPos);
			new_value[0] = value[0];
		}
	}
	if( forbidden != NULL ) if( !forbidden->empty() ) cout << "XMLTag- USING: " << path << "\t\t" << name << "\t\t" << new_value[0] << endl;
	if( value.empty() )
	{
	//	cerr << "Requested value of tag " << name << ", but the value is empty" << endl;
	//	throw(-765);
		return vector<string>();
	}
	else
	{
		return new_value;
	}
}

//Find first level tags in the content
vector< XMLTag* > XMLTag::FindTagsInContent( vector<string> Content, vector<string> & Value )
{
	vector< XMLTag* > childTags;

	while (true)
	{
		//Find the first tag
		int lineNumber, linePosition;
		string tagName = FindNextTagOpen( Content, lineNumber, linePosition );

		//Create the value
		if ( lineNumber == -1 && linePosition == -1 )
		{
			//If no tag found, all content goes into value
			StringProcessing::RemoveWhiteSpace(Content);
			for ( unsigned int contentIndex = 0; contentIndex < Content.size(); ++contentIndex )
			{
				Value.push_back( Content[contentIndex] );
			}
			break;
		}
		else
		{
			//Everything before the tag goes into value
			vector<string> newValue = this->SplitContent( Content, lineNumber, linePosition, int(tagName.size() + 2) );
			StringProcessing::RemoveWhiteSpace(newValue);
			for ( unsigned int contentIndex = 0; contentIndex < newValue.size(); ++contentIndex )
			{
				Value.push_back( newValue[contentIndex] );
			}

			//Create the child tag
			vector<string> childContent = this->FindTagContent( tagName, Content );
			childTags.push_back( new XMLTag( tagName, childContent, this ) );
		}
	}

	return childTags;
}

//Find the next tag opening in content, and return the name
string XMLTag::FindNextTagOpen( vector<string> Content, int & LineNumber, int & LinePosition )
{
	//Loop over all lines
	for ( unsigned int lineIndex = 0; lineIndex < Content.size(); ++lineIndex )
	{
		//	Ignore lines starting with '#'
		if( Content[lineIndex].find('#')!=0 )
		{
			int openPosition = StringProcessing::CharacterPosition( Content[lineIndex], '<' );
			if ( openPosition > -1 )
			{
				int closePosition = StringProcessing::CharacterPosition( Content[lineIndex], '>' );
				if ( closePosition > openPosition )
				{
					//Tag is complete
					string this_name( Content[lineIndex], unsigned(openPosition + 1), unsigned(closePosition - openPosition - 1) );

					//Error check
					if ( this_name.size() == 0 )
					{
						cerr << "Found tag with no name in line: \"" << Content[lineIndex] << "\"" << endl;
						exit(1);
					}

					if ( this_name[0] == '/' )
					{
						cerr << "Found a closing tag <" << this_name << "> when opening tag expected" << endl;
						exit(1);
					}

					LineNumber = int(lineIndex);
					LinePosition = openPosition;
					return this_name;
				}
				else
				{
					//Incomplete tag
					cerr << "Incomplete XML tag in line: \"" << Content[lineIndex] << "\"" << endl;
					exit(1);
				}
			}
		}
	}

	//If you get to the end, no tags found
	LineNumber = -1;
	LinePosition = -1;
	return "";
}

//Find all instances of tags with this name opening in the content
void XMLTag::FindTagOpens( string TagName, vector<string> Content, vector<int> & LineNumbers, vector<int> & LinePositions )
{
	LineNumbers.clear();
	LinePositions.clear();
	string openTag = "<" + TagName + ">";

	//Loop over all lines
	for ( unsigned int lineIndex = 0; lineIndex < Content.size(); ++lineIndex )
	{
		vector<int> positions = StringProcessing::StringPositions( Content[lineIndex], openTag );
		for ( unsigned int positionIndex = 0; positionIndex < positions.size(); ++positionIndex )
		{
			LineNumbers.push_back(int(lineIndex));
			LinePositions.push_back( positions[positionIndex] );
		}
	}
}

//Find all instances of tag with this name closing in the content
void XMLTag::FindTagCloses( string TagName, vector<string> Content, vector<int> & LineNumbers, vector<int> & LinePositions )
{
	LineNumbers.clear();
	LinePositions.clear();
	string closeTag = "</" + TagName + ">";

	//Loop over all lines
	for ( unsigned int lineIndex = 0; lineIndex < Content.size(); ++lineIndex )
	{
		vector<int> positions = StringProcessing::StringPositions( Content[lineIndex], closeTag );
		for ( unsigned int positionIndex = 0; positionIndex < positions.size(); ++positionIndex )
		{
			LineNumbers.push_back(int(lineIndex));
			LinePositions.push_back( positions[positionIndex] );
		}
	}

	if ( LineNumbers.size() < 1 )
	{
		cerr << "Tag " << TagName << " is not closed" << endl;
		exit(1);
	}
}

//Find the content of a tag. Note that this method seems unnecessarily complicated because it must allow tags to have the same name as their parent tag
vector<string> XMLTag::FindTagContent( string TagName, vector<string> & Content )
{
	vector<int> openLines, openPositions, closeLines, closePositions;
	FindTagOpens( TagName, Content, openLines, openPositions );
	FindTagCloses( TagName, Content, closeLines, closePositions );
	int splitIndex=0;;

	//Check to make sure tag order is preserved
	for ( unsigned int closeIndex = 0; closeIndex < closeLines.size(); ++closeIndex )
	{
		unsigned int opens = 0;
		for ( unsigned int openIndex = 0; openIndex < openLines.size(); ++openIndex )
		{
			if ( openLines[openIndex] < closeLines[closeIndex] )
			{
				//New tag of same name is opened before this close
				++opens;
			}
			else if ( openLines[openIndex] == closeLines[closeIndex] && openPositions[openIndex] < closePositions[closeIndex] )
			{
				//New tag of same name is opened before this close on same line
				++opens;
			}
			else
			{
				//New tag is not opened before this close
				break;
			}
		}

		if ( opens == closeIndex )
		{
			//Tag is complete
			splitIndex = int(closeIndex);
			break;
		}
	}

	//Split off the tag content
	return this->SplitContent( Content, closeLines[unsigned(splitIndex)], closePositions[unsigned(splitIndex)], int(TagName.size() + 3) );
}

//Split a vector string about the specified place
vector<string> XMLTag::SplitContent( vector<string> & Content, int LineNumber, int LinePosition, int GapLength )
{
	vector<string> head, tail;

	for ( int lineIndex = 0; lineIndex < int(Content.size()); ++lineIndex )
	{
		if ( lineIndex > LineNumber )
		{
			//All lines before the split line go in "head"
			tail.push_back( Content[unsigned(lineIndex)] );
		}
		else if ( lineIndex < LineNumber )
		{
			//All lines after the split line go in "tail"
			head.push_back( Content[unsigned(lineIndex)] );
		}
		else
		{
			//The first part of the split line
			string begin( Content[unsigned(lineIndex)], 0, unsigned(LinePosition) );
			if ( begin != "" )
			{
				head.push_back(begin);
			}

			//The second part of the split line
			string end( Content[unsigned(lineIndex)], unsigned(LinePosition + GapLength) );
			if ( end != "" )
			{
				tail.push_back(end);
			}
		}
	}

	Content = tail;
	return head;
}


bool XMLTag::GetBooleanValue( const XMLTag* input )
{
	string value = input->GetValue()[0].c_str();
	return strcasecmp( value.c_str() , "TRUE" ) == 0 ;
}

int XMLTag::GetIntegerValue( const XMLTag* input )
{
	return atoi( input->GetValue()[0].c_str() );
}

string XMLTag::GetStringValue( const XMLTag* input )
{
	return input->GetValue()[0];
}

double XMLTag::GetDoubleValue( const XMLTag* input )
{
	return strtod( input->GetValue()[0].c_str(), NULL );
}

double XMLTag::GetTF1Eval( const XMLTag* input )
{
	TF1* thisFunc = new TF1( "someFunc", input->GetValue()[0].c_str() );
	double returnable = thisFunc->Eval( 0 );
	delete thisFunc;
	return returnable;
}

void XMLTag::PrePendPath( const string input )
{
	path = parent->GetPath();
	path.Append( input.c_str() );
	path.Append( "/" + name );
	for( unsigned int i=0; i< children.size(); ++i )
	{
		children[i]->RegeneratePath();
	}
}

void XMLTag::AppendPath( const string input )
{
	path.Append( input.c_str() );
	for( unsigned int i=0; i< children.size(); ++i )
	{
		children[i]->RegeneratePath();
	}
}

void XMLTag::RegeneratePath()
{
	TString parent_path;
	if( parent != NULL ) parent_path = parent->GetPath();
	path = parent_path;
	path.Append( "/" + name );
	for( unsigned int i=0; i< children.size(); ++i )
	{
		children[i]->RegeneratePath();
	}
}

