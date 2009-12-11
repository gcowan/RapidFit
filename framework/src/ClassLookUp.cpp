/**
        @class ClassLookUp

        Central place to hold the methods for returning an instance of a class with a given name.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "ClassLookUp.h"
#include "RaPDF_Bs2DsPi.h"
#include "RaPDF_Bs2JpsiPhi.h"
#include "RaPDF_Bs2JpsiPhiNew.h"
#include "RaPDF_Bs2JpsiPhi_withTimeRes.h"
#include "RaPDF_Bs2JpsiPhi_sWave.h"
#include "RaPDF_Bs2JpsiPhiMassSignal.h"
#include "RaPDF_Bs2JpsiPhiLongLivedBkg.h"
#include "Bs2JpsiPhiLongLivedBkg_withTimeRes.h"
#include "RaPDF_Bs2JpsiPhiPromptBkg.h"
#include "Bs2JpsiPhiPromptBkg_withTimeRes.h"
#include "RaPDF_Bs2JpsiPhiMassBkg.h"
#include "RaPDF_Bs2PhiPhi.h"
#include "MinuitWrapper.h"
#include "Minuit2Wrapper.h"
#include "FumiliWrapper.h"
#include "NegativeLogLikelihood.h"
#include "Foam.h"
#include "AcceptReject.h"
#include "JPsiPhiDataGenerator.h"
#include <stdlib.h>

//Look up the name of a PDF, return an appropriate instance of IPDF
IPDF * ClassLookUp::LookUpPDFName( string Name, vector<string> PDFObservables, vector<string> PDFParameters )
{
        if ( Name == "RaPDF_Bs2JpsiPhiNew" )
        {
	        //Default JPsiPhi
	        return new RaPDF_Bs2JpsiPhiNew();
        }
        else if ( Name == "RaPDF_Bs2PhiPhi" )
        {
	        //Default PhiPhi
                return new RaPDF_Bs2PhiPhi();
        }
	else if ( Name == "RaPDF_Bs2JpsiPhi_withTimeRes" )
        {
                // Bs2JPsiPhi with analytic time resolution
                return new RaPDF_Bs2JpsiPhi_withTimeRes();
        }
        else if ( Name == "RaPDF_Bs2DsPi" )
        {
                // DsPi
                return new RaPDF_Bs2DsPi();
        }
	else if ( Name == "RaPDF_Bs2JpsiPhiMassSignal" )
        {
                //Default JPsiPhi signal mass PDF
                return new RaPDF_Bs2JpsiPhiMassSignal();
        }
        else if ( Name == "RaPDF_Bs2JpsiPhiLongLivedBkg" )
        {
	        //Long lived background for JPsiPhi
                return new RaPDF_Bs2JpsiPhiLongLivedBkg();
        }
        else if ( Name == "Bs2JpsiPhiLongLivedBkg_withTimeRes" )
        {
                //Long lived background for JPsiPhi with time resolution (convolved gaussian)
                return new Bs2JpsiPhiLongLivedBkg_withTimeRes();
        }
        else if ( Name == "RaPDF_Bs2JpsiPhiPromptBkg" )
        {
        	// This one does not work at the moment. Use time resolution.
		cerr << "This one does not work at the moment. Use time resolution" << endl;
	        return new RaPDF_Bs2JpsiPhiPromptBkg();
        }

        else if ( Name == "Bs2JpsiPhiPromptBkg_withTimeRes" )
        {
	        //Prompt background for JPsiPhi, with time resolution (convolved gaussian)
		return new Bs2JpsiPhiPromptBkg_withTimeRes();
        }

	else if ( Name == "RaPDF_Bs2JpsiPhi_sWave" )
        {
	        //JPsiPhi signal PDF including the s-wave contribution
                return new RaPDF_Bs2JpsiPhi_sWave();
        }
        else if ( Name == "RaPDF_Bs2JpsiPhiMassBkg" )
        {
                //Default JPsiPhi prompt bkg mass signal
                return new RaPDF_Bs2JpsiPhiMassBkg();
	}
	else
	{
		cerr << "Unrecognised PDF name: " << Name << endl;
		exit(1);
	}
}

//Look up the name of a fit function, and return an appropriate instance
FitFunction * ClassLookUp::LookUpFitFunctionName( string Name )
{
	if ( Name == "NegativeLogLikelihood" )
	{
		return new NegativeLogLikelihood();
	}
	else
	{
		cerr << "Unrecognised function to minimise: " << Name << endl;
		exit(1);
	}
}

//Look up the name of a minimiser, and return an appropriate instance
IMinimiser * ClassLookUp::LookUpMinimiserName( string Name, int NumberParameters )
{
	if ( Name == "Minuit" )
	{
		return new MinuitWrapper(NumberParameters);
	}
	else if ( Name == "Minuit2" )
	{
		return new Minuit2Wrapper();
	}
	else if ( Name == "Fumili" )
        {
                return new FumiliWrapper();
        }
	else
	{
		cerr << "Unrecognised minimiser name: " << Name << endl;
		exit(1);
	}
}

//Look up the name of a data generator, and return an appropriate instance
IDataGenerator * ClassLookUp::LookUpDataGenerator( string Name, PhaseSpaceBoundary * Boundary, IPDF * Function )
{
	if ( Name == "Foam" )
	{
		return new Foam( Boundary, Function );
	}
	else if ( Name == "AcceptReject" )
	{
		return new AcceptReject( Boundary, Function );
	}
	else if ( Name == "JPsiPhiDataGenerator" )
	{
		return new JPsiPhiDataGenerator( Boundary, Function );
	}
	else
	{
		cerr << "Unrecognised data generator name: " << Name << endl;
		exit(1);
	}
}
