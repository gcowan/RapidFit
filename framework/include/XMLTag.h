/**
        @class XMLTag

        A class used to parse an xml config file. Represents a single xml tag and its contents.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef XML_TAG_H
#define XML_TAG_H

//	System Headers
#include <vector>
#include <string>

using namespace std;

class XMLTag
{
	public:
		XMLTag();
		XMLTag( string, vector<string> );
		~XMLTag();

		static vector< XMLTag* > FindTagsInContent( vector<string>, vector<string>& );
		string GetName();
		vector< XMLTag* > GetChildren();
		vector<string> GetValue();

	private:
		static string FindNextTagOpen( vector<string>, int&, int& );
		static void FindTagOpens( string, vector<string>, vector<int>&, vector<int>& );
		static void FindTagCloses( string, vector<string>, vector<int>&, vector<int>& );
		static vector<string> FindTagContent( string, vector<string>& );
		static vector<string> SplitContent( vector<string>&, int, int, int );

		vector< XMLTag* > children;
		vector<string> value;
		string name;
};

#endif
