/**
        @class XMLTag

        A class used to parse an xml config file. Represents a single xml tag and its contents.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "XMLTag.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <stdlib.h>

//Default Constructor
XMLTag::XMLTag() : children(), value(), name("RapidFit"), parent(NULL), path(""), forbidden()
{
}

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

vector<pair<string,string> >* XMLTag::GetForbidden()
{
	return forbidden;
}

//Return the complete Path of the Object within the XML Structure
string XMLTag::GetPath()
{
	return string( path.Data() );
}

//Return the tag name
string XMLTag::GetName()
{
	string new_name = name;
	return new_name;
}

//Return the children
vector< XMLTag* > XMLTag::GetChildren()
{
	string NameStr="Name";
	vector<string> children_names;
	int NamePos=0;
	for( vector<XMLTag*>::iterator child_i = children.begin(); child_i != children.end(); ++child_i )
	{
		children_names.push_back( (*child_i)->GetName() );
	}
	NamePos = StringProcessing::VectorContains( &children_names, &NameStr );
	if( NamePos != -1  )
	{
		size_t found = string::npos;
		found = string( path.Data() ).find( children[NamePos]->GetValue()[0] );
		if( found == string::npos )
		{
			path.Append( "/" + children[NamePos]->GetValue()[0] );
		}
	}
	vector<XMLTag*> new_children = children;
	//cout << path << endl;
	return new_children;
}

//Return the value
vector<string> XMLTag::GetValue()
{
	vector<string> new_value = value;
	path = parent->GetPath();
	path.Append( "/" + name );
	//cout << path << endl;
	string pathStr = path.Data();
	vector<string> forbidden_paths;
	for( vector<pair<string,string> >::iterator fpath_i = forbidden->begin(); fpath_i != forbidden->end(); ++fpath_i )
	{
		forbidden_paths.push_back( fpath_i->first );
	}
	int TagPos = StringProcessing::VectorContains( &forbidden_paths, &pathStr );
	if( TagPos != -1 )
	{
		new_value[0] = (*forbidden)[TagPos].second;
		forbidden->erase(forbidden->begin()+TagPos);
	}
	//cout << path << "\t" << name << "\t" << "\t" << new_value[0] << endl;
	if ( value.size() == 0 )
	{
		cerr << "Requested value of tag " << name << ", but the value is empty" << endl;
		exit(1);
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
			vector<string> newValue = SplitContent( Content, lineNumber, linePosition, int(tagName.size() + 2) );
			StringProcessing::RemoveWhiteSpace(newValue);
			for ( unsigned int contentIndex = 0; contentIndex < newValue.size(); ++contentIndex )
			{
				Value.push_back( newValue[contentIndex] );
			}

			//Create the child tag
			vector<string> childContent = FindTagContent( tagName, Content );
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
					string name( Content[lineIndex], unsigned(openPosition + 1), unsigned(closePosition - openPosition - 1) );
	
					//Error check
						if ( name.size() == 0 )
					{
						cerr << "Found tag with no name in line: \"" << Content[lineIndex] << "\"" << endl;
						exit(1);
					}
		
					if ( name[0] == '/' )
					{
						cerr << "Found a closing tag <" << name << "> when opening tag expected" << endl;
						exit(1);
					}

					LineNumber = int(lineIndex);
					LinePosition = openPosition;
					return name;
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
	return SplitContent( Content, closeLines[unsigned(splitIndex)], closePositions[unsigned(splitIndex)], int(TagName.size() + 3) );
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
