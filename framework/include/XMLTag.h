/*!
 * @class XMLTag
 *
 * @brief A class used to parse an xml config file. Represents a single xml tag and its contents.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef XML_TAG_H
#define XML_TAG_H

//	ROOT Headers
#include "TString.h"
//	System Headers
#include <vector>
#include <string>

using namespace::std;

class XMLTag
{
	public:
		/*!
		 * @brief Constructor Which accepts (Optional) Ovveride XML arguments from the command line
		 *
		 * @param Input   This contains a vector of pairs with each pair containing an absolute path of an XML tag and a new value to take
		 */
		XMLTag( vector<pair<string,string> >* Input =NULL );

		/*!
		 * @brief Constructor used internally to Generate a new XMLTag which has a knowledge of it's parent and can have children of it's own
		 *
		 * @param Name      Name of the new XMLTag
		 * @param Children  Either children or contents of Tag
		 * @param Parent    Pointer to the parent of this tag
		 */
		XMLTag( string Name, vector<string> Children, XMLTag* Parent = NULL );

		/*!
		 * @brief Standard Destructor
		 */
		~XMLTag();

		/*!
		 * @brief Returns a vector of Tags found within the content
		 *
		 * @param Content   The raw text to be passed to the interprator
		 * @param Value     Populated if no new tag is found or everything before the next open tag
		 *
		 * @return          A vector of all children tags within Content (empty if none found)
		 */
		vector<XMLTag*> FindTagsInContent( const vector<string> Content, vector<string>& Value );

		/*!
		 * @brief Get the path of this XML tag relative to the top level of the XML File
		 *
		 * @return the path of the XML Tag in the form '/top/middle/bottom/this'
		 */
		string GetPath() const;

		/*!
		 * @brief Modify the Internal path of this XMLTag with a string in front of the tag name
		 *
		 * @input String that is to be pre-pended to the internal path
		 */
		void PrePendPath( const string input );

		/*!
		 * @brief Modify the Internal path of this XMLTag with a string after the name of this tag
		 *
		 * @input String that is to be appended to the internal path
		 */
		void AppendPath( const string input );

		/*!
		 * @brief Return the Name of this Tag
		 *
		 * @return returns the name of this XMLTag as a string
		 */
		string GetName() const;

		/*!
		 * @brief Return the Override Tags and the override values
		 *
		 * @return a pointer to the list of override tags and the values that they should take
		 */
		vector<pair<string,string> >* GetForbidden() const;

		/*!
		 * @brief Return a vector of pointers to the children Tags of this Tag
		 *
		 * @return  Vector containing pointers to the children, empty if none
		 */
		vector< XMLTag* > GetChildren() const;

		/*!
		 * @brief Remove a chosen child from an XMLTag
		 *
		 * @param num This is the number of the XMLTag to be removed from this tag
		 *
		 * @return void
		 */
		void RemoveChild( int num );

		/*!
		 * @brief This method will evaluate the first string content as a boolean and return a c++ bool
		 *
		 * @brief boolean calculated from case insensitive check for the string true in a string
		 */
		static bool GetBooleanValue( const XMLTag* input );

		/*!
		 * @brief This method will evaluate the first string content as an integer value and return a c++ integer
		 *
		 * @return a signed integer based on the first string content in XMLTag
		 */
		static int GetIntegerValue( const XMLTag* input );

		/*!
		 * @brief This method will return the first string content of this XMLTag as the value of this XMLTag
		 *
		 * @return a string which is the first string in the XMLTag
		 */
		static string GetStringValue( const XMLTag* input );

		/*!
		 * @brief This method will evaluate the first string content as a double value and return a c++ double
		 *
		 * @return a double based on the first string content in XMLTag
		 */
		static double GetDoubleValue( const XMLTag* input );

		/*!
		 * @brief This method will evaluate the first string content as a TF1 expression value and return a c++ double
		 *
		 * @return a double based on the first string content in the XMLTag evaluated as a TF1
		 */
		static double GetTF1Eval( const XMLTag* input );

		vector<string> GetRAWValue() const;
	private:
		/*!
		 * Don't Copy the class this way!
		 */
		XMLTag ( const XMLTag& );

		/*!
		 * Don't Copy the class this way!
		 */
		XMLTag& operator = ( const XMLTag& );

		/*!
		 * @brief Finds the next open tag that is not precedded by a '#' character
		 */
		string FindNextTagOpen( vector<string>, int&, int& );

		/*!
		 * @brief Find all instances of tags with this name opening in the content
		 */
		void FindTagOpens( string, vector<string>, vector<int>&, vector<int>& );

		/*!
		 * @brief Find all instances of tag with this name closing in the content
		 */
		void FindTagCloses( string, vector<string>, vector<int>&, vector<int>& );

		/*!
		 * @brief Internal method to trigger the re-generation of the path inside this XMLTag
		 */
		void RegeneratePath();

		/*!
		 * @brief Find the content of a tag.
		 *
		 * Note that this method seems unnecessarily complicated because it must allow tags to have the same name as their parent tag
		 */
		vector<string> FindTagContent( string, vector<string>& );

		/*!
		 * @brief Split a vector string about the specified place
		 */
		vector<string> SplitContent( vector<string>&, int, int, int );

		vector<XMLTag*> children;			/*!	Children within this XML Tag		*/
		mutable vector<string> value;			/*!	Value of this XML Tag			*/
		string name;					/*!	Name of this XML Tag			*/
		XMLTag* parent;					/*!	Address of the parent of this tag	*/
		mutable TString path;				/*!	Path of this XML Tag			*/
		vector<pair<string, string> >* forbidden;	/*!	override xml path and values		*/

		/*!
		 * @brief Get the Value of this Tag
		 *
		 * This will be multi valued for tags which have multiple values on multiple lines
		 *
		 * @return Return a list of the vector of Values
		 */
		vector<string> GetValue() const;
};

#endif

