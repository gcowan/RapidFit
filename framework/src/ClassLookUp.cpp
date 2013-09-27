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
#include "IPDF.h"
#include "NormalisedSumPDF.h"
#include "ProdPDF.h"
#include "SumPDF.h"
#include "MinuitWrapper.h"
#include "Minuit2Wrapper.h"
#include "FumiliWrapper.h"
#include "NegativeLogLikelihood.h"
#include "NegativeLogLikelihoodThreaded.h"
#include "NegativeLogLikelihoodThreadedNew.h"
#include "NegativeLogLikelihoodNumerical.h"
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
	//#include "dummy.h"		//	This needs correcting
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
		//	ROOT equivalent to dlopen
		returnable_object = (void*) gSystem->DynFindSymbol( RAPIDFIT_LIBRARY_NAME, Name.c_str() );
	}
	else
	{
		//	POSIX for *nix based systems
		//	There should probably be a wrapper here for windows based builds but we don't anticipate this
		void *handle = dlopen( ClassLookUp::getSelfPath(), RTLD_LAZY | RTLD_NODELETE );
		returnable_object = dlsym( handle, Name.c_str() );
	}

	if( returnable_object == NULL )
	{
		cerr << "Error:\t\tCan't find symbol List or Object:\t" << Name << endl << endl;
	}

	return returnable_object;
}

//	Given a PDF object name we know that the PDF constructor is wrapped in an object with a derrived name
IPDF* ClassLookUp::LookUpPDFName( string Name, PDFConfigurator* configurator )
{
	string pdf_creator_Name = "CreatePDF_"+Name;

	CreatePDF_t* pdf_creator = (CreatePDF_t*) ClassLookUp::getObject( pdf_creator_Name );

	if( pdf_creator == NULL )
	{
		cerr << "Cannot Find PDF Named: " << Name << " You probably provided the wrong name in your XML." << endl << endl;
		exit(-15513);
	}

	IPDF* returnable_PDF = (IPDF*) pdf_creator( configurator );

	returnable_PDF->SetConfigurator( configurator );

	string thisLabel = configurator->GetPDFLabel();

	if( !thisLabel.empty() ) returnable_PDF->SetLabel( thisLabel );

	return returnable_PDF;
}

IPDF* ClassLookUp::CopyPDF( const IPDF* inputPDF )
{
	IPDF* returnable_PDF=NULL;

	string Name = inputPDF->GetName();

	if( inputPDF->IsCopyConstructorSafe() )
	{

		//	The PDF already knows where it's template is, no need to look for it again
		if( inputPDF->GetCopyConstructor() != NULL )
		{
			CopyPDF_t* pdf_copy = inputPDF->GetCopyConstructor();
			returnable_PDF = (IPDF*) pdf_copy( *inputPDF );
		}
		else
		{
			//	These special case PDFs explicitly need to be declared here as they're special case objects heavily integrated into the framework
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
				//	Each PDF has a C wrapper function with an unmangled name of CopyPDF_SomePDF
				string pdf_copy_Name = "CopyPDF_"+Name;

				//	Find this object in the object which has been loaded as a library
				CopyPDF_t* pdf_copy = (CopyPDF_t*) ClassLookUp::getObject( pdf_copy_Name );

				if( pdf_copy == NULL )
				{
					cerr << "Cannot Find Copy Constructor for: " << Name << " This is highly unsual and was thought to be impossible... I cannot continue" << endl;
					cerr << "You may have used SetName and incorrectly named the PDF, Use SetLabel if this is the case." << endl << endl;
					exit(21841);
				}

				//	Give the PDF object explicit knowledge of the path of it's constructor template
				returnable_PDF = (IPDF*) pdf_copy( *inputPDF );
				inputPDF->SetCopyConstructor( pdf_copy );
			}
		}
	}
	else
	{
		PDFConfigurator* thisConfig = new PDFConfigurator( *(inputPDF->GetConfigurator()) );
		if( !thisConfig->hasConfigurationValue( "RAPIDFIT_SAYS_THIS_IS_A_COPY", "True" ) )
		{
			if( thisConfig->hasConfigurationValue( "RAPIDFIT_SAYS_THIS_IS_A_COPY", "False" ) )
			{
				// OK Someone is intentionally messing around with the framework's behaviour!
			}
			else
			{
				thisConfig->addConfigurationParameter( "RAPIDFIT_SAYS_THIS_IS_A_COPY:True" );
			}
		}
		returnable_PDF = ClassLookUp::LookUpPDFName( Name, thisConfig ); 
		inputPDF->SetCopyConstructor( NULL );
		delete thisConfig;
	}

	returnable_PDF->Can_Remove_Cache( false );

	return returnable_PDF;
}

//Look up the name of a fit function, and return an appropriate instance
IFitFunction * ClassLookUp::LookUpFitFunctionName( string Name )
{
	if ( Name == "NegativeLogLikelihood" )
	{
		IFitFunction* returnable = (IFitFunction*) new NegativeLogLikelihood();
		returnable->SetThreads(0);
		return returnable;
	}
	else if ( Name == "NegativeLogLikelihoodThreaded" )
	{
		IFitFunction* returnable = (IFitFunction*) new NegativeLogLikelihoodThreaded();
		returnable->SetThreads(-1);
		return returnable;
	}
	else if ( Name == "NegativeLogLikelihoodThreadedNew" )
	{
		IFitFunction* returnable = (IFitFunction*) new NegativeLogLikelihoodThreadedNew();
		returnable->SetThreads(-1);
		return returnable;
	}
	else if ( Name == "NegativeLogLikelihoodNumerical" )
	{
		IFitFunction* returnable = (IFitFunction*) new NegativeLogLikelihoodNumerical();
		returnable->SetThreads(-1);
		return returnable;
	}
	else
	{
		cerr << "Unrecognised function to minimise: " << Name << endl << endl;
		exit(-2355);
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
		cerr << "Unrecognised minimiser name: " << Name << endl << endl;
		exit(23478);
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
		cerr << "Unrecognised data generator name: " << Name << endl << endl;
		exit(366);
	}
}

//Look up the name of a precalculator, and return an appropriate instance
IPrecalculator * ClassLookUp::LookUpPrecalculator( string Name, string WeightName, FitResult* inputResult, unsigned int config )
{
	if ( Name == "SWeightPrecalculator" )
	{
		return new SWeightPrecalculator( inputResult, WeightName, config );
	}
	else
	{
		cerr << "Unrecognised precalculator name: " << Name << endl << endl;
		exit(-2138);
	}
}

