/*!
 * This class is intended as a simple interface for studies that people will use RapidFit for
 * It provides a template of 'essential' functions that every study should implement without being too strict
 *
 * This is more of a template of how to do something systematically within RapidFit and interface it to the rest of the framework
 *
 * This should be strictly split into BaseStudy and IStudy, but I am lazy and didn't want to generate a BaseStudy for a few methods, but that should be easy if required
 */

#pragma once
#ifndef IStudy_H
#define IStudy_H

///	RapidFit Headers
#include "XMLConfigReader.h"
#include "FitResultVector.h"
#include "ParameterSet.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
#include "ConstraintFunction.h"
#include "PDFWithData.h"
///	System Headers
#include <vector>
#include <string>

using namespace::std;

class IStudy
{
	public:
		/*!
		 * @brief Can't have virtual public Constructors, and it doesn't make sense to either
		 *        By default assume we're only being passed references to objects, and as such do not need to look after them internally
		 */
		virtual ~IStudy()
		{
			if( delete_objects == true )
			{
				while( !pdfsAndData.empty() ) { if( pdfsAndData.back() != NULL ) { delete pdfsAndData.back(); } pdfsAndData.pop_back(); }
				if( studyParameters != NULL ) { delete studyParameters; }
				if( theMinimiser != NULL ) delete theMinimiser;
				if( allResults != NULL ) delete allResults;	//Potentially causes SERIOUS PROBLEMS with unchecked code in RapidFit!!!
				while( !allConstraints.empty() ) { if( allConstraints.back() != NULL ) { delete allConstraints.back(); } allConstraints.pop_back(); }
			}
			if( debug != NULL ) delete debug;
		};

		/*!
		 * @brief Interface Function to call the Study to now actually start
		 *
		 * @param Verbosity (optional) This provides a hook to change how verbose we want the study to be, useful for debugging something
		 *
		 * @return Void
		 */
		virtual void DoWholeStudy( int Verbosity = -999 ) = 0;

		/*!
		 * @brief Interface Function to call the Study to get the results of the Study
		 *
		 * @return Returns a FitResultVector instance which contains all of the Output from the individual fits in the Study
		 */
		virtual FitResultVector* GetStudyResult() = 0;

		/*!
		 * @brief Interface Function to call the Study to change the Number of Repeats that should be performed
		 */
		virtual void SetNumRepeats( int ) = 0;

		/*!
		 * @brief This Allows for CommandLine defined Parameters to be introduced into the study here
		 *
		 * @param Input   These are a vector of the raw CommandLine arguments to define a Parameter that is important to the Study
		 *
		 * NOT essential, but all Studies in theory can implement this to meet their own ends
		 *
		 * @return Void
		 */
		virtual void SetCommandLineParams( vector<string> Input ) = 0;

		/*!
		 * @brief Interface (and implementation) to Print some more verbose Debugging Information on this study
		 *
		 * @return Void
		 */
		virtual void Print()
		{
			cout << "Study Performed Using: " << endl << endl;
			cout << pdfsAndData.size() << " pdfsAndData: " << endl;
			for( unsigned int i=0; i< pdfsAndData.size(); ++i ) { pdfsAndData[i]->Print(); cout << endl; }
			cout << endl << "studyParameters: " << endl;
			studyParameters->Print(); cout << endl;
			cout << "Minimser: " << endl; theMinimiser->Print(); cout << endl;
			cout << "Results: " << endl; allResults->Print(); cout << endl;
			cout << allConstraints.size() << " ConstraintFunction: " << endl;
			for( unsigned int i=0; i< allConstraints.size(); ++i ) { allConstraints[i]->Print(); cout << endl; }
			cout << numberStudies << " numberStudies " << endl;
		};

		virtual void SetDebug( DebugClass* input_debug )
		{
			if( debug != NULL ) delete debug;
			debug = new DebugClass( *input_debug );
		}

	private:
		/*!
		 * Don't Copy the class this way!
		 */
		IStudy ( const IStudy& );
		/*!
		 * Don't Copy the class this way!
		 */
		IStudy& operator= ( const IStudy& );

	protected:
		/*!
		 *	All Studies should provide a minimum of these objects
		 *	Yes, this shouldn't technically be in an interface class, however I'm trying to create a common code structure with this
		 *
		 *      Having these objects here allows for a much easier to maintain codebase if all future studies can use this naming scheme
		 */

		/*!
		 * This is a vector of Pointers to the DataSet(s) and the PDFs we wish to evaluate them with
		 */
		vector<PDFWithData*> pdfsAndData;

		/*!
		 * This a pointer to the ParameterSet being used within the Study
		 */
		ParameterSet* studyParameters;

		/*!
		 * This is the configuration class for the Minimiser which is to be passed to FitAssembler for minimisations
		 */
		MinimiserConfiguration* theMinimiser;

		/*!
		 * This is the configuration class for the FitFunction which is to be passed to the FitAssembler for minimisations
		 */
		FitFunctionConfiguration* theFunction;

		/*!
		 * This is the object which should contain all of the FitResults during the fit
		 */
		FitResultVector* allResults;

		/*!
		 * This is a vector of External Constraints which are to be passed to the FitAssembler for minimisations
		 */
		vector<ConstraintFunction*> allConstraints;

		/*!
		 * This is the XML file which was used to initialize the whole study. It's a bit bad to rely on this to generate objects but it does make life easier
		 */
		XMLConfigReader* xmlConfig;

		/*!
		 * This is the stored number of Repeats which were requested for this study
		 */
		int numberStudies;

		/*!
		 * Should the IStudy baseclass clean up after you?
		 * Default of true makes the IStudy destructor cleanup all of the internally looked after objects, false means you wish to write your own proper destructor
		 */
		bool delete_objects;

		DebugClass* debug;

		/*!
		 * Provide a default constructor to initialize the objects to NULL or empty,
		 *
		 * again shouldn't be in an interface, but less work than a BaseStudy simply for (con/de)structors
		 */
		IStudy() :
			pdfsAndData(), studyParameters(), theMinimiser(NULL), theFunction(NULL), allResults(NULL),
			allConstraints(), numberStudies(-1), delete_objects(false), xmlConfig(NULL), debug(new DebugClass(false) )
	{};
};

#endif

