/**
        @class ClassLookUp

        Central place to hold the methods for returning an instance of a class with a given name.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "ClassLookUp.h"
#include "Bs2JpsiPhi.h"
#include "Bs2JpsiPhi_mistagObservable.h"
#include "Bs2JpsiPhi_mistagObservable_withTimeRes.h"
#include "Bs2JpsiPhi_mistagObservable_withAngAcc.h"
#include "Bs2JpsiPhi_mistagObservable_withAverageAngAcc.h"
#include "Bs2JpsiPhi_sWave.h"
#include "Bs2JpsiPhi_mistagParameter.h"
#include "Bs2JpsiPhi_mistagParameter_withTimeRes.h"
#include "Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.h"
#include "Bs2JpsiPhi_mistagParameter_alt.h"
#include "Bs2JpsiPhiMassSignal.h"

#include "Bd2JpsiKstar_withTimeRes_withAverageAngAcc.h"
#include "Bd2JpsiKstar_sWave.h"
#include "Bs2DsPi.h"
#include "Bs2DsPi_mistagParameter.h"
#include "Bs2DsPi_acc.h"
#include "Bs2PhiPhi.h"
#include "Bs2DsPiMassSignal.h"
#include "Bs2DsPiBkg_withTimeRes.h"

#include "Bs2JpsiPhiLongLivedBkg.h"
#include "Bs2JpsiPhiLongLivedBkg_withTimeRes.h"
#include "Bs2JpsiPhiPromptBkg.h"
#include "Bs2JpsiPhiPromptBkg_withTimeRes.h"
#include "Bs2JpsiPhiPromptBkg_withTimeResDouble.h"
#include "Bs2JpsiPhiMassBkg.h"

#include "MinuitWrapper.h"
#include "Minuit2Wrapper.h"
#include "FumiliWrapper.h"
#include "NegativeLogLikelihood.h"
#include "Foam.h"
#include "AcceptReject.h"
#include "JPsiPhiDataGenerator.h"
#include "SWeightPrecalculator.h"
#include <stdlib.h>

//Look up the name of a PDF, return an appropriate instance of IPDF
IPDF * ClassLookUp::LookUpPDFName( string Name, vector<string> PDFObservables, vector<string> PDFParameters )
{
        if ( Name == "Bs2JpsiPhi_mistagObservable" )
        {
	        //Default JPsiPhi
	        return new Bs2JpsiPhi_mistagObservable();
        }
	else if ( Name == "Bs2JpsiPhi_mistagObservable_withTimeRes" )
        {
                // Bs2JPsiPhi with analytic double gaussian time resolution
                return new Bs2JpsiPhi_mistagObservable_withTimeRes();
        }
        else if ( Name == "Bs2JpsiPhi_mistagObservable_withAngAcc" )
        {
                //JpsiPhi with angular acceptance fed in as "observables"
                return new Bs2JpsiPhi_mistagObservable_withAngAcc();
        }
        else if ( Name == "Bs2JpsiPhi_mistagObservable_withAverageAngAcc" )
        {
                //JpsiPhi with angular acceptance fed in as fixed physics parameters
                return new Bs2JpsiPhi_mistagObservable_withAverageAngAcc();
        }
	else if ( Name == "Bs2JpsiPhi_mistagParameter" )
        {
	        //Default JPsiPhi with mistag as a physics parameter
	        return new Bs2JpsiPhi_mistagParameter();
        }
	else if ( Name == "Bs2JpsiPhi_mistagParameter_withTimeRes" )
        {
	        //Default JPsiPhi with mistag as physics parameter and time res on
	        return new Bs2JpsiPhi_mistagParameter_withTimeRes();
        }
	else if ( Name == "Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc" )
        {
	        //Default JPsiPhi with mistag as physics parameter, average ang acc and time res on
	        return new Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc();
        }
	else if ( Name == "Bs2JpsiPhi_mistagParameter_alt" )
        {
	        //JPsiPhi from Pete with mistag as a physics paramter and single gaussian time res
	        return new Bs2JpsiPhi_mistagParameter_alt();
        }
        else if ( Name == "Bs2PhiPhi" )
        {
	        //Default PhiPhi
                return new Bs2PhiPhi();
        }
        else if ( Name == "Bs2DsPi" )
        {
                // DsPi
                return new Bs2DsPi();
        }
        else if ( Name == "Bs2DsPi_mistagParameter" )
        {
			// DsPi
			return new Bs2DsPi_mistagParameter();
        }
        else if ( Name == "Bs2DsPi_acc" )
        {
			// DsPi w. acceptance
			return new Bs2DsPi_acc();
        }
        else if ( Name == "Bs2DsPiMassSignal" )
        {
                //Default DsPi signal mass PDF
                return new Bs2DsPiMassSignal();
        }
		else if ( Name == "Bs2JpsiPhiMassSignal" )
        {
                //Default JPsiPhi signal mass PDF
                return new Bs2JpsiPhiMassSignal();
        }
        else if ( Name == "Bs2DsPiBkg_withTimeRes" )
        {
		// DsPi
		return new Bs2DsPiBkg_withTimeRes();
        }
		else if ( Name == "Bs2JpsiPhiLongLivedBkg" )
        {
	        //Long lived background for JPsiPhi
                return new Bs2JpsiPhiLongLivedBkg();
        }
        else if ( Name == "Bs2JpsiPhiLongLivedBkg_withTimeRes" )
        {
                //Long lived background for JPsiPhi with time resolution (convolved gaussian)
                return new Bs2JpsiPhiLongLivedBkg_withTimeRes();
        }
        else if ( Name == "Bs2JpsiPhiPromptBkg" )
        {
        	// This one does not work at the moment. Use time resolution.
		cerr << "This one does not work at the moment. Use time resolution" << endl;
		exit(1);
	        //return new Bs2JpsiPhiPromptBkg();
        }
        else if ( Name == "Bs2JpsiPhiPromptBkg_withTimeRes" )
        {
	        //Prompt background for JPsiPhi, with time resolution (convolved gaussian)
		return new Bs2JpsiPhiPromptBkg_withTimeRes();
        }
        else if ( Name == "Bs2JpsiPhiPromptBkg_withTimeResDouble" )
        {
	        //Prompt background for JPsiPhi, with time resolution (double convolved gaussian)
			return new Bs2JpsiPhiPromptBkg_withTimeResDouble();
        }
		else if ( Name == "Bs2JpsiPhi_sWave" )
        {
	        //JPsiPhi signal PDF including the s-wave contribution
                return new Bs2JpsiPhi_sWave();
        }
        else if ( Name == "Bs2JpsiPhiMassBkg" )
        {
                //Default JPsiPhi prompt bkg mass signal
                return new Bs2JpsiPhiMassBkg();
	}
		else if ( Name == "Bd2JpsiKstar_withTimeRes_withAverageAngAcc" )
        {
			// Bd2JPsiKstar with analytic double gaussian time resolution
			return new Bd2JpsiKstar_withTimeRes_withAverageAngAcc();
        }

                 else if ( Name == "Bd2JpsiKstar_sWave" )
         {
                         // Bd2JPsiKstar with analytic double gaussian time resolution and sWave
                         return new Bd2JpsiKstar_sWave();
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

//Look up the name of a precalculator, and return an appropriate instance
IPrecalculator * ClassLookUp::LookUpPrecalculator( string Name, IPDF * FirstPDF, IPDF * SecondPDF, ParameterSet * FitParameters, string WeightName )
{
	if ( Name == "SWeightPrecalculator" )
	{
		return new SWeightPrecalculator( FirstPDF, SecondPDF, FitParameters, WeightName );
	}
	else
	{
		cerr << "Unrecognised precalculator name: " << Name << endl;
		exit(1);
	}
}
