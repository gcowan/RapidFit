/**
        @class RooJPsiPhiWrapper

        A wrapper for a RooFit format PDF describing Bs decay to J/Psi Phi

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "RooJPsiPhiWrapper.h"
#include "RooPdf_Bs2JPsiPhi.h"

//Default constructor
RooJPsiPhiWrapper::RooJPsiPhiWrapper()
{
	prototypeDataPoint.push_back( "time" );
	prototypeDataPoint.push_back( "ctheta_tr" );
	prototypeDataPoint.push_back( "phi_tr" );
	prototypeDataPoint.push_back( "ctheta_1" );

	prototypeParameterSet.push_back( "gamma" );
	prototypeParameterSet.push_back( "dgam" );
	prototypeParameterSet.push_back( "Rt" );
	prototypeParameterSet.push_back( "Rp" );
	prototypeParameterSet.push_back( "delta1" );
	prototypeParameterSet.push_back( "delta2" );
	prototypeParameterSet.push_back( "delta_ms" );
	prototypeParameterSet.push_back( "phi_s" );
	prototypeParameterSet.push_back( "tagFraction" );
	prototypeParameterSet.push_back( "resolution" );

	if ( SetUpWrapper( prototypeDataPoint, 4, prototypeParameterSet, 10 ) )
	{
		wrappedPDF = new RooPdf_Bs2JPsiPhi("JPsiPhi", "Bs to J psi phi PDF", 1,//Dummy for "btype". Doesn't seem to be used
				*observables[0], *observables[1], *observables[2], *observables[3],
				*variables[0], *variables[1], *variables[2], *variables[3],
				*variables[4], *variables[5], *variables[6], *variables[7],
				*variables[8], *variables[9]);
		valid = true;
	}
	else
	{
		valid = false;
	}
}

//Constructor with correct arguments
RooJPsiPhiWrapper::RooJPsiPhiWrapper( vector<string> PDFObservables, vector<string> PDFParameters)
{
	prototypeDataPoint = PDFObservables;
	prototypeParameterSet = PDFParameters;

	if ( SetUpWrapper( PDFObservables, 4, PDFParameters, 10 ) )
	{
		wrappedPDF = new RooPdf_Bs2JPsiPhi("JPsiPhi", "Bs to J psi phi PDF", 1,//Dummy for "btype". Doesn't seem to be used
				*observables[0], *observables[1], *observables[2], *observables[3],
				*variables[0], *variables[1], *variables[2], *variables[3],
				*variables[4], *variables[5], *variables[6], *variables[7],
				*variables[8], *variables[9]);
		valid = true;
	}
	else
	{
		valid = false;
	}
}

//Destructor
RooJPsiPhiWrapper::~RooJPsiPhiWrapper()
{
}
