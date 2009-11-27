/**
        @class ClassLookUp

        Central place to hold the methods for returning an instance of a class with a given name.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "ClassLookUp.h"
//#include "RooGaussianWrapper.h"
//#include "RooJPsiPhiWrapper.h"
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
#include "InvalidObject.h"
#include "MinuitWrapper.h"
#include "Minuit2Wrapper.h"
#include "FumiliWrapper.h"
#include "NegativeLogLikelihood.h"
#include "Foam.h"
#include "AcceptReject.h"
#include "JPsiPhiDataGenerator.h"

//Look up the name of a PDF, return an appropriate instance of IPDF
IPDF * ClassLookUp::LookUpPDFName( string Name, vector<string> PDFObservables, vector<string> PDFParameters )
{
	/*
	if ( Name == "RooGaussianWrapper" )
	{
		if ( PDFObservables.size() == 0 && PDFParameters.size() == 0 )
		{
			//Default gaussian
			return new RooGaussianWrapper();
		}
		else
		{
			return new RooGaussianWrapper( PDFObservables, PDFParameters );
		}
	}
	else if ( Name == "RooJPsiPhiWrapper" )
	{
		if ( PDFObservables.size() == 0 && PDFParameters.size() == 0 )
		{
			//Default JPsiPhi
			return new RooJPsiPhiWrapper();
		}
		else
		{
			return new RooJPsiPhiWrapper( PDFObservables, PDFParameters );
		}
	}
	*/
        if ( Name == "RaPDF_Bs2JpsiPhiNew" )
        {
                //if ( PDFObservables.size() == 0 && PDFParameters.size() == 0 )
                //{
                        //Default JPsiPhi
                        return new RaPDF_Bs2JpsiPhiNew();
                //}
        }
	if ( Name == "RaPDF_Bs2JpsiPhi_withTimeRes" )
        {
                // Bs2JPsiPhi with analytic time resolution
                return new RaPDF_Bs2JpsiPhi_withTimeRes();
        }
        if ( Name == "RaPDF_Bs2DsPi" )
        {
                // DsPi
                return new RaPDF_Bs2DsPi();
        }
        

	else if ( Name == "RaPDF_Bs2JpsiPhiMassSignal" )
        {
                //if ( PDFObservables.size() == 0 && PDFParameters.size() == 0 )
                //{
                        //Default JPsiPhi signal mass PDF
                        return new RaPDF_Bs2JpsiPhiMassSignal();
                //}
        }
        else if ( Name == "RaPDF_Bs2JpsiPhiLongLivedBkg" )
        {
                        return new RaPDF_Bs2JpsiPhiLongLivedBkg();
        }
        else if ( Name == "Bs2JpsiPhiLongLivedBkg_withTimeRes" )
        {
                        return new Bs2JpsiPhiLongLivedBkg_withTimeRes();
        }
        else if ( Name == "RaPDF_Bs2JpsiPhiPromptBkg" )
        {
        	// This one does not work at the moment. Use time resolution.
	        return new RaPDF_Bs2JpsiPhiPromptBkg();
        }

        else if ( Name == "Bs2JpsiPhiPromptBkg_withTimeRes" )
        {
		return new Bs2JpsiPhiPromptBkg_withTimeRes();
        }

	else if ( Name == "RaPDF_Bs2JpsiPhi_sWave" )
        {
                        return new RaPDF_Bs2JpsiPhi_sWave();
        }
        else if ( Name == "RaPDF_Bs2JpsiPhiMassBkg" )
        {
                //if ( PDFObservables.size() == 0 && PDFParameters.size() == 0 )
                //{
                        //Default JPsiPhi prompt bkg
                        return new RaPDF_Bs2JpsiPhiMassBkg();
                //}
        }
	else
	{
		cerr << "Unrecognised PDF name: " << Name << endl;
		return new InvalidObject( "Unrecognised PDF name: " + Name );
	}
}

//Look up the name of a fit function, and return an appropriate instance
FitFunction * ClassLookUp::LookUpFitFunctionName( string Name, string Weight )
{
	if ( Name == "NegativeLogLikelihood" )
	{
		if ( Weight == "" )
		{
			return new NegativeLogLikelihood();
		}
		else
		{
			return new NegativeLogLikelihood(Weight);
		}
	}
	else
	{
		cerr << "Unrecognised function to minimise: " << Name << endl;
		return new InvalidObject( "Unrecognised function name: " + Name );
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
		return new InvalidObject( "Unrecognised minimiser name: " + Name );
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
		return new InvalidObject( "Unrecognised data generator name: " + Name );
	}
}
