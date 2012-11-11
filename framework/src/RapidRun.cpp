//	ROOT Headers
#include "TList.h"
#include "TObject.h"
#include "TObjString.h"
#include "TString.h"
//	RapidFit Headers
#include "main.h"
#include "RapidRun.h"
#include "XMLConfigReader.h"
//	System Headers
#include <vector>
#include <string>
#include <iostream>

using std::cout;
using std::endl;

ClassImp( RapidRun );//needed for Cint

RapidRun::RapidRun() : args()
{
}

RapidRun::RapidRun( TList* _args ) : args( _args )
{
}

int RapidRun::run()
{
	//convert the TList args into a main compatible list
	//const Int_t argc = args->GetSize();
	//memory is dynamically allocated
	//char** argv = new char*[size_t(argc)];

	vector<string> argv_str;

	TObjLink *lnk = args->FirstLink();
	Int_t count = 0;
	while (lnk) {
		string this_str = string( ((TObjString*)lnk->GetObject())->String() );
		cout << "String from List " << this_str << endl;
		argv_str.push_back( this_str );
		cout << "String Storeds in argv: " << argv_str.back() << endl;
		lnk = lnk->Next();
		++count;
	}
	cout << "RapidRun Launching RapidFit!" <<endl;


	vector<char*> argv;
	for( unsigned int i=0; i< argv_str.size(); ++i )
	{
		const char* this_char = argv_str[i].c_str();
		argv.push_back( const_cast<char*>(this_char) );
	}

	Int_t result = RapidFit( (int)argv.size(), &(argv[0]) );

	//	Cleanup
	//for( int i=0; i<argc; ++i )
	//{
	//	delete argv[i];
	//}
	//delete argv;

	return result;
}

