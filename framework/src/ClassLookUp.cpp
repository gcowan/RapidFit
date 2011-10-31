/**
  @class ClassLookUp

  Central place to hold the methods for returning an instance of a class with a given name.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	ROOT Headers
#include "TSystem.h"

//	RapidFit Headers
#include "ClassLookUp.h"
#include "BasePDF.h"
#include "NormalisedSumPDF.h"
#include "ProdPDF.h"
#include "SumPDF.h"
#include "MinuitWrapper.h"
#include "Minuit2Wrapper.h"
#include "FumiliWrapper.h"
#include "NegativeLogLikelihood.h"
#include "NegativeLogLikelihoodThreaded.h"
#include "Foam.h"
#include "AcceptReject.h"
#include "JPsiPhiDataGenerator.h"
#include "SWeightPrecalculator.h"
#include "RapidRun.h"

//	System Headers
#include <stdlib.h>
#include <dlfcn.h>
#include <fstream>
#include <iostream>
#ifdef _WIN32
	//#include "dummy.h"	//	This needs correcting
#elif __APPLE__
	#include <mach-o/dyld.h>
	#include <syslimits.h>
#else
	#include <linux/limits.h>
	#include <unistd.h>
	#include <stdint.h>
#endif

//	Get the path of the running executable on this system
//	Arch independent, and an example of the windows code has not been dreamt up
//	(and I don't care to unles someone throws their toys out of the pram!)
char* ClassLookUp::getSelfPath()
{
	//void* handle = NULL;
	static char filename[PATH_MAX];
	uint32_t pathsize = PATH_MAX;
	(void) pathsize;

	#ifdef _WIN32
		cout << "Windows hasn't been implemented yet... any takers?" << endl;
		exit(-42);
	#elif __APPLE__
		_NSGetExecutablePath(filename, &pathsize);
	#else
		ssize_t len = readlink("/proc/self/exe", filename, sizeof(filename)-1);
		if (len != -1) {
			filename[len] = '\0';
		} else {
			filename[0] = 'F' ;
			filename[1] = 'A' ;
			filename[2] = 'I' ;
			filename[3] = 'L' ;
			filename[4] = '\0';
		}
	#endif
	
	return filename;
}

//	Resolve the object 'Name' within the binary object containing RapidFit
void* ClassLookUp::getObject( string Name )
{
	void* returnable_object;

	string pathName = ClassLookUp::getSelfPath();
	string root_exe = "root.exe";
	string failed_str = "FAIL";


	//	If RapidFit is compiled/run as a CINT library
	//	DO NOT INTEROGATE THE OUTSIDE WORLD (you probably won't like what you find on the grid...)
	//	This needs to be protected with a macro due to ROOT needlessley SCREAMING at you with unwanted errors
	if( pathName.find( root_exe ) != string::npos )
	{
		returnable_object = (void*) gSystem->DynFindSymbol( RAPIDFIT_LIBRARY_NAME, Name.c_str() );
	} else {
		void *handle = dlopen( ClassLookUp::getSelfPath(), RTLD_LAZY | RTLD_GLOBAL );
		returnable_object = dlsym( handle, Name.c_str() );
	}

	if( returnable_object == NULL )
	{
		cerr << "Error:\t\tCan't find symbol List or Object:\t" << Name << endl << endl;
		exit(-45);
	}

	return returnable_object;
}

IPDF* ClassLookUp::LookUpPDFName( string Name, vector<string> PDFObservables, vector<string> PDFParameters, PDFConfigurator* configurator )
{
	(void) PDFObservables;
	(void) PDFParameters;

	string pdf_creator_Name = "CreatePDF_"+Name;

	CreatePDF_t* pdf_creator = (CreatePDF_t*) ClassLookUp::getObject( pdf_creator_Name );

	IPDF* returnable_PDF = pdf_creator( configurator );

	returnable_PDF->SetName( Name ); 

	return returnable_PDF;
}

IPDF* ClassLookUp::CopyPDF( IPDF* inputPDF )
{
	IPDF* returnable_PDF=NULL;
	string Name = inputPDF->GetName();

	if( Name == "NormalisedSum" )
	{
		returnable_PDF = (IPDF*) new NormalisedSumPDF( *(const NormalisedSumPDF*) inputPDF );
	}
	else if( Name == "Prod" )
	{
		returnable_PDF = (IPDF*) new ProdPDF( *(const ProdPDF*) inputPDF );
	}
	else if( Name == "Sum" )
	{
		returnable_PDF = (IPDF*) new SumPDF( *(const SumPDF*) inputPDF );
	}
	else
	{
		string pdf_copy_Name = "CopyPDF_"+Name;

		CopyPDF_t* pdf_copy = (CopyPDF_t*) ClassLookUp::getObject( pdf_copy_Name );

		returnable_PDF = pdf_copy( *inputPDF );
	}

	returnable_PDF->SetName( Name );
	returnable_PDF->Remove_Cache( true );

	if( !inputPDF->GetActualParameterSet()->GetAllNames().empty() )
	{
		returnable_PDF->SetPhysicsParameters( inputPDF->GetActualParameterSet() );
		returnable_PDF->UpdateIntegralCache();
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
	if ( Name == "NegativeLogLikelihoodThreaded" )
	{
		return new NegativeLogLikelihoodThreaded();
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
