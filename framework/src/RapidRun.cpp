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
#include <cstring>
#include <iostream>

using std::cout;
using std::endl;

ClassImp( RapidRun )//needed for Cint

RapidRun::RapidRun( TList* _args ) : args( _args )
{
}

int RapidRun::run(){
	//convert the TList args into a main compatible list
	const Int_t argc = args->GetSize();
	//memory is dynamically allocated
	char** argv = new char*[size_t(argc)];

	TObjLink *lnk = args->FirstLink();
	Int_t count = 0;
	while (lnk) {
		const char* str = ((TObjString*)lnk->GetObject())->String();
		argv[ count ] = new char[ strlen(str) ];
		strcpy( argv[ count ], str);
		lnk = lnk->Next();
		count++;
	}
	cout << "RapidRun Launching RapidFit!" <<endl;

	Int_t result = RapidFit(count, argv);

	//	Cleanup
	for( int i=0; i<argc; ++i )
	{
		delete argv[i];
	}
	delete argv;

	return result;
}

