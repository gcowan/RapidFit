
#pragma once
#ifndef RAPIDRUN_H
#define RAPIDRUN_H


//	ROOT Headers
#include "Rtypes.h"
#include "TList.h"
#include "TObject.h"
#include "TObjString.h"
#include "TString.h"
//	RapidFit Headers
#include "main.h"
#include "Mathematics.h"
//	System Headers
#include <memory>

#define RAPIDFIT_LIBRARY_NAME "libRapidRun"

class RapidRun : public TObject
{
	public:
		RapidRun(); //needed by Root IO
		RapidRun( TList* );
		int run();

		static void setGridification( bool input );
		static bool isGridified();

	private:
		std::auto_ptr<TList> args;

		static bool runningOnGrid;

	public:
		ClassDef( RapidRun, 2 ); //Needed for Cint
};

#endif

