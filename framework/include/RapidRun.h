#ifndef RAPIDRUN_H
#define RAPIDRUN_H


//	ROOT Headers
#include "TList.h"
#include "TObject.h"
#include "TObjString.h"
#include "TString.h"
//	RapidFit Headers
#include "main.h"
//	System Headers
#include <memory>


class RapidRun : public TObject
{
	public:
		RapidRun() : args(){}//needed by Root IO
		RapidRun( TList* );
		int run();

	private:
		std::auto_ptr<TList> args;
		ClassDef( RapidRun, 1 )//Needed for Cint
};

#endif
