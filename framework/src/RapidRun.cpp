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
#include <sstream>

using std::cout;
using std::endl;

ClassImp( RapidRun );//needed for Cint

RapidRun::RapidRun() : args()
{
}

RapidRun::RapidRun( TList* _args ) : args( _args )
{
}

bool RapidRun::isGridified()
{
	return runningOnGrid;
}

void RapidRun::setGridification( bool input )
{
	runningOnGrid = input;
}

bool RapidRun::runningOnGrid = false;

int RapidRun::run()
{
	cout << "Running Start" << endl;

	RapidRun::setGridification( true );

	vector<string> argv_str;
	TObjLink *lnk = args->FirstLink();
	while (lnk)
	{
		TString thisTString = ((TObjString*)lnk->GetObject())->GetString();
		cout << "String from List " << thisTString.Data() << endl;
		stringstream thisStream;
		thisStream << thisTString.Data();
		argv_str.push_back( thisStream.str() );
		cout << "String Stored in argv: " << argv_str.back() << endl;
		lnk = lnk->Next();
	}
	cout << "RapidRun Launching RapidFit!" <<endl;

	Int_t result = RapidFit( argv_str );

	return result;
}

