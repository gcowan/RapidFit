
#ifndef __XMLFUNCTIONS_UTILS_H
#define __XMLFUNCTIONS_UTILS_H

#include "TObject.h"

#include <string>
#include <vector>

using namespace::std;

class XMLUtilFunctions : public TObject
{

	public:

		static void RestoreXML( vector<string> input_filenames, vector<string> other_params );

		static void GetToyXML( vector<string> input_filenames, vector<string> other_params );

		static void GetProjectionXML( vector<string> input_filenames, vector<string> other_params );

		static void Print();

		XMLUtilFunctions() {};
		virtual ~XMLUtilFunctions() {};
		ClassDef( XMLUtilFunctions, 1 );

};

#endif

