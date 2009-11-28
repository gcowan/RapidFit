// $Id: FumiliWrapper.h,v 1.1 2009/11/10 10:35:44 gcowan Exp $
/**
        @class FumiliWrapper

        A wrapper for Fumili, implementing IMinimiser

	I think this has defeated me, there is no easy way to get
	this to work. With Fumili you need to calculate the
	derivatives (there is a Gradient() method that needs to be
	defined in FumiliFCNBase()) which is in general difficult.
	Here is an example of how this minimisation technique can 
	be used:

	http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/

	This is using it without ROOT, just as we are using Minuit2.
	There are some advanced ways to use fitters in ROOT (see 
	TVirtualFitter and related classes).

	http://seal.web.cern.ch/seal/MathLibs/Minuit2/html/annotated.html
	
        @author Greig A Cowan greig.cowan@cern.ch
	@date 2009-10-09
*/

#ifndef FUMILI_WRAPPER_H
#define FUMILI_WRAPPER_H

#include "IMinimiser.h"
#include "Minuit2/MnMigrad.h"
#include "FumiliFunction.h"
#include "FitResult.h"

using namespace ROOT::Minuit2;

class FumiliWrapper : public IMinimiser
{
	public:
		FumiliWrapper();
		~FumiliWrapper();

		//Interface functions
		virtual void Minimise( FitFunction* );
		virtual FitResult * GetFitResult();
		virtual void ContourPlots( vector< pair< string, string > > );

	private:
		//MnMigrad minuit;
		FumiliFunction * function;
		FitResult * fitResult;
		vector< pair< string, string > > contours;
};

#endif
