/**
        @class XMLTag

        A class used to parse an xml config file. Represents a single xml tag and its contents.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef XML_TAG_H
#define XML_TAG_H

//	ROOT Headers
#include "TString.h"
//	System Headers
#include <vector>
#include <string>

using namespace std;

class XMLTag
{
	public:
		XMLTag();
		XMLTag( vector<pair<string,string> >* );
		XMLTag( string, vector<string>, XMLTag* );
		~XMLTag();

		vector< XMLTag* > FindTagsInContent( vector<string>, vector<string>& );
		string GetPath();
		void AppendPath( string );
		string GetName();
		vector<pair<string,string> >* GetForbidden();
		vector< XMLTag* > GetChildren();
		vector<string> GetValue();
		bool NotAppended();

	private:
                //      Uncopyable!
                XMLTag ( const XMLTag& );
                XMLTag& operator = ( const XMLTag& );

		string FindNextTagOpen( vector<string>, int&, int& );
		void FindTagOpens( string, vector<string>, vector<int>&, vector<int>& );
		void FindTagCloses( string, vector<string>, vector<int>&, vector<int>& );
		vector<string> FindTagContent( string, vector<string>& );
		vector<string> SplitContent( vector<string>&, int, int, int );

		vector< XMLTag* > children;
		vector<string> value;
		string name;
		XMLTag* parent;
		TString path;
		vector<pair<string, string> >* forbidden;
};

#endif
