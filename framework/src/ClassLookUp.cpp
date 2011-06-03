/**
        @class ClassLookUp

        Central place to hold the methods for returning an instance of a class with a given name.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "ClassLookUp.h"
// Signal PDFs set 1  (originally Greig and Conor mostly)
#include "Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.h"
#include "Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave.h"
#include "Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms.h"

#include "Bs2JpsiPhi_SignalAlt_MP_v1.h"
#include "Bs2JpsiPhi_SignalAlt_MP_dev.h"
#include "Bs2JpsiPhi_SignalAlt_MO_v1.h"
#include "Bs2JpsiPhi_SignalAlt_MO_dev.h"


#include "Bs2JpsiPhiMassSignal.h"

#include "Bd2JpsiKstar_withTimeRes_withAverageAngAcc.h"
#include "Bd2JpsiKstar_sWave.h"
#include "Bd2JpsiKstar_sWave_rTerms.h"
#include "Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms.h"

#include "Bs2DsPi.h"
#include "Bs2DsPi_mistagParameter.h"
#include "Bs2DsPi_acc.h"
#include "Bs2PhiPhi.h"
#include "Bs2DsPiMassSignal.h"
#include "Bs2DsPiBkg_withTimeRes.h"

#include "Bs2JpsiPhiLongLivedBkg.h"
#include "Bs2JpsiPhiLongLivedBkg_II.h"
#include "Bs2JpsiPhiLongLivedBkg_withTimeRes.h"
#include "Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist.h"

#include "Bs2JpsiPhiPromptBkg_withTimeRes.h"
#include "Bs2JpsiPhiPromptBkg_withTimeResDouble.h"
#include "Bs2JpsiPhiPromptBkg_tripleGaussian.h"
#include "Bs2JpsiPhiMassBkg.h"
#include "Bs2JpsiPhiMassBkgLL.h"
#include "LongLivedBkg_3Dangular.h"

#include "Exponential.h"

#include "MinuitWrapper.h"
#include "Minuit2Wrapper.h"
#include "FumiliWrapper.h"
#include "NegativeLogLikelihood.h"
#include "Foam.h"
#include "AcceptReject.h"
#include "JPsiPhiDataGenerator.h"
#include "SWeightPrecalculator.h"

//	System Headers
#include <stdlib.h>

//Look up the name of a PDF, return an appropriate instance of IPDF
IPDF * ClassLookUp::LookUpPDFName( string Name, vector<string> PDFObservables, vector<string> PDFParameters, PDFConfigurator configurator )
{

	IPDF* returnable_PDF = NULL;

	vector<string> null_vec = PDFObservables;
	vector<string> null_vec2 = PDFParameters;
	while( !null_vec.empty() ) null_vec.pop_back();   // what are these lines for ???
	while( !null_vec2.empty() ) null_vec2.pop_back();

		if ( Name == "Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc" )
        {
	        //Default JPsiPhi
		returnable_PDF = new Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc();
        }

		else if ( Name == "Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave" )
        {
	        //Default JPsiPhi with s wave
	        returnable_PDF = new Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave();
        }
		else if ( Name == "Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms" )
        {
	        //Default JPsiPhi with s wave
	        returnable_PDF = new Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms();
        }
		else if ( Name == "Bs2JpsiPhi_SignalAlt_MP_v1" )
        {
	        //JPsiPhi from Pete with sWave
		returnable_PDF = new  Bs2JpsiPhi_SignalAlt_MP_v1();
        }
		else if ( Name == "Bs2JpsiPhi_SignalAlt_MP_dev" )
        {
	        //JPsiPhi from Pete with sWave
	        returnable_PDF = new Bs2JpsiPhi_SignalAlt_MP_dev(configurator);
        }
		else if ( Name == "Bs2JpsiPhi_SignalAlt_MO_v1" )
        {
	        //JPsiPhi from Pete with sWave
	        returnable_PDF = new Bs2JpsiPhi_SignalAlt_MO_v1();
        }
		else if ( Name == "Bs2JpsiPhi_SignalAlt_MO_dev" )
        {
	        //JPsiPhi from Pete with sWave
	        returnable_PDF = new Bs2JpsiPhi_SignalAlt_MO_dev(configurator);
        }
        else if ( Name == "Bs2PhiPhi" )
        {
	        //Default PhiPhi
                returnable_PDF = new Bs2PhiPhi();
        }
        else if ( Name == "Bs2DsPi" )
        {
                // DsPi
                returnable_PDF = new Bs2DsPi();
        }
        else if ( Name == "Bs2DsPi_mistagParameter" )
        {
			// DsPi
			returnable_PDF = new Bs2DsPi_mistagParameter();
        }
        else if ( Name == "Bs2DsPi_acc" )
        {
			// DsPi w. acceptance
			returnable_PDF = new Bs2DsPi_acc();
        }
        else if ( Name == "Bs2DsPiMassSignal" )
        {
                //Default DsPi signal mass PDF
                returnable_PDF = new Bs2DsPiMassSignal();
        }
		else if ( Name == "Bs2JpsiPhiMassSignal" )
        {
                //Default JPsiPhi signal mass PDF
                returnable_PDF = new Bs2JpsiPhiMassSignal();
        }
        else if ( Name == "Bs2DsPiBkg_withTimeRes" )
        {
		// DsPi
		returnable_PDF = new Bs2DsPiBkg_withTimeRes();
        }
		else if ( Name == "Bs2JpsiPhiLongLivedBkg" )
        {
	        //Long lived background for JPsiPhi
                returnable_PDF = new Bs2JpsiPhiLongLivedBkg();
        }
		else if ( Name == "Bs2JpsiPhiLongLivedBkg_II" )
        {
	        //Long lived background for JPsiPhi with different parameters
			returnable_PDF = new Bs2JpsiPhiLongLivedBkg_II();
        }
        else if ( Name == "Bs2JpsiPhiLongLivedBkg_withTimeRes" )
        {
                //Long lived background for JPsiPhi with time resolution (convolved gaussian)
                returnable_PDF = new Bs2JpsiPhiLongLivedBkg_withTimeRes();
        }
        else if ( Name == "Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist" )
        {
			//Long lived background for JPsiPhi with time resolution (convolved gaussian)
			returnable_PDF = new Bs2JpsiPhiLongLivedBkg_withTimeRes_withAngDist();
        }
        else if ( Name == "Bs2JpsiPhiPromptBkg_withTimeRes" )
        {
	        //Prompt background for JPsiPhi, with time resolution (convolved gaussian)
		returnable_PDF = new Bs2JpsiPhiPromptBkg_withTimeRes();
        }
        else if ( Name == "Bs2JpsiPhiPromptBkg_withTimeResDouble" )
        {
	        //Prompt background for JPsiPhi, with time resolution (double convolved gaussian)
			returnable_PDF = new Bs2JpsiPhiPromptBkg_withTimeResDouble();
        }
        else if ( Name == "Bs2JpsiPhiPromptBkg_tripleGaussian" )
        {
	        //Prompt background for JPsiPhi, with time resolution (double convolved gaussian)
			returnable_PDF = new Bs2JpsiPhiPromptBkg_tripleGaussian();
        }
        else if ( Name == "Bs2JpsiPhiMassBkg" )
        {
                //Default JPsiPhi prompt bkg mass signal
                returnable_PDF = new Bs2JpsiPhiMassBkg();
	}
        else if ( Name == "Exponential" )
        {
            	//Default JPsiPhi prompt bkg mass signal
                returnable_PDF = new Exponential();
	}
        else if ( Name == "Bs2JpsiPhiMassBkgLL" )
        {
			//Default JPsiPhi prompt bkg mass signal
			returnable_PDF = new Bs2JpsiPhiMassBkgLL();
		}

	else if ( Name == "Bd2JpsiKstar_withTimeRes_withAverageAngAcc" )
        {
			// Bd2JPsiKstar with analytic double gaussian time resolution
			returnable_PDF = new Bd2JpsiKstar_withTimeRes_withAverageAngAcc();
        }

        else if ( Name == "Bd2JpsiKstar_sWave" )
        {
                         // Bd2JPsiKstar with analytic double gaussian time resolution and sWave
                         returnable_PDF = new Bd2JpsiKstar_sWave();
        }

	else if ( Name == "Bd2JpsiKstar_sWave_rTerms" )
        {
                         // Bd2JPsiKstar with analytic double gaussian time resolution and sWave
                         returnable_PDF = new Bd2JpsiKstar_sWave_rTerms();
        }


	else if ( Name == "Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms" )
        {
                         // Bd2JPsiKstar with analytic double gaussian time resolution and sWave
                         returnable_PDF = new Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms();
        }

	else if ( Name == "LongLivedBkg_3Dangular" )
        {
                         returnable_PDF = new LongLivedBkg_3Dangular(configurator);
        }

	else
	{
		cerr << "Unrecognised PDF name: " << Name << endl;
		exit(1);
	}
	return returnable_PDF;
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
IPrecalculator * ClassLookUp::LookUpPrecalculator( string Name, IPDF * FirstPDF, IPDF * SecondPDF, vector<ParameterSet*> FitParameters, string WeightName )
{
	if ( Name == "SWeightPrecalculator" )
	{
		return new SWeightPrecalculator( FirstPDF, SecondPDF, FitParameters.back(), WeightName );
	}
	else
	{
		cerr << "Unrecognised precalculator name: " << Name << endl;
		exit(1);
	}
}
