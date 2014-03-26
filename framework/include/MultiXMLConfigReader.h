/*!
 * @class XMLConfigReader
 *
 * @brief Opens an xml config file and uses it to create RapidFit data objects
 *
 * This can be thought of as the initialization engine of RapidFit
 *
 * It is capable of calling the correct constructors and assembling the correct objects for passing back into the main code
 *
 *
 * @ warning XMLConfigurationReader DOES NOT UNDER ANY CIRCUMSTANCES LOOK AFTER THE OBJECTS IT PASSES A POINTER FOR
 *                                  clean up after yourself!
 *
 *
 * This does rely on XMLTag but in principle it can be exchanged for a proper XML parser
 *
 * Benefits?	Unknown, none of us are XML guru's
 * 
 * Cons?        Can't comment out code as quickly/flexibly
 *              Would be difficult to engineer --OverrideXML into a new parser it barely works here ;)
 *              Would have to be installed on all platforms and the authors have a distrust of ROOT XML
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef Multi_XML_CONFIG_READER_H
#define Multi_XML_CONFIG_READER_H

//	RapidFit Headers
#include "XMLTag.h"
#include "ParameterSet.h"
#include "PDFWithData.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
#include "OutputConfiguration.h"
#include "DataSetConfiguration.h"
#include "ComponentPlotter.h"
#include "ConstraintFunction.h"
#include "ScanParam.h"
#include "PrecalculatorConfig.h"
#include "DebugClass.h"
#include "XMLConfigReader.h"
//	System Headers
#include <vector>
#include <string>

using namespace::std;

class MultiXMLConfigReader : public I_XMLConfigReader
{
	public:

		MultiXMLConfigReader( vector<string> fileNames, vector<pair<string, string> >* OverrideXML = NULL );

		/*!
		 * Default Destructor
		 */
		~MultiXMLConfigReader();

		/*!
		 * @brief Returns the result of a boolean decision to determine if the XML has passed some basic sanity checks for having a complete configuration
		 *
		 * @return True if the file can be used by RapidFit, False if it cannot
		 */
		bool IsValid() const;

		/*!
		 * @brief Return the Input XML from the file
		 *
		 * @return This has removed all whitespaces from the XML so output should be treated as a template
		 */
		vector<string> GetXML() const;

		/*!
		 * @brief Returns the ParameterSet for the whole Fit.
		 *
		 * When an XML File contains multiple ParameterSets it will attempt to merge them into one large ParameterSet
		 * with one of each unique Parameter
		 *
		 * When Parameters are Provided at the Command Line they will Override the existing Parameters in the XML or will Add to it
		 *
		 * @param CommandLineParam  This parameter contains strings from Runtime to override or add a Physics Parameter to the ParameterSet
		 *
		 * @return Returns a Pointer to the ParameterSet object and gives up any control over it.
		 */
		ParameterSet* GetFitParameters( vector<string> CommandLineParam = vector<string>() );

		/*!
		 * @brief Returns the MinimiserConfiguration class as defined in the XML
		 *
		 * @return Returns a Pointer to the MinimiserConfiguration class which is capable of Creating the Minimiser required for the Fit
		 */
		MinimiserConfiguration* GetMinimiserConfiguration();

		/*!
		 * @brief Returns the FitFunctionConfiguration class as defined in the XML
		 *
		 * @return Returns a Pointer to the FitFunctionConfiguration class which is capable of Creating the FitFunction required for the Fit
		 */
		FitFunctionConfiguration* GetFitFunctionConfiguration();

		/*!
		 * @brief Returns the OutputConfiguration as defined in the XML
		 *
		 * @return Returns a Pointer to the OutputConfiguration which produces more detaild output around the fit minima
		 */
		OutputConfiguration* GetOutputConfiguration();

		/*!
		 * @brief Returns a vector of PDFWithData objects, one for each ToFit containing a PDF in the XML
		 *
		 * @param StartingValues   This contains the starting values passed from the command line to control where in the file the fits read in data
		 *
		 * @return Returns pointers to the different PDFWithData objects present in the XML
		 */
		vector<PDFWithData*> GetPDFsAndData( vector<int> StartingValues = vector<int>() );

		/*!
		 * @brief Returns the list of External Constrains as defined in the XML
		 *
		 * @returns a vector of separate ConstraintFunctions as defined in the XML
		 */
		vector<ConstraintFunction* > GetConstraints();

		/*!
		 * @brief Returns a vector of the PhaseSpaceBoundaries associated with the PDFs in the XML
		 */
		vector<PhaseSpaceBoundary*> GetPhaseSpaceBoundaries();

		/*!
		 * @brief Construct the PrecalculatorConfig Object in order to Weight a DataSet using an algorithm
		 */
		PrecalculatorConfig* GetPrecalculatorConfig();

		/*!
		 * @brief Get the size of each of the DataSets
		 */
		vector<int> GetAllDataSetSizes();

		/*!
		 * @brief Get the Start Entries of Each DataSet in the XML
		 */
		vector<int> GetAllStartEntries();

		/*!
		 * @brief Get the Number of Repeats as defined in the top level of the XML
		 */
		int GetNumberRepeats();

		/*!
		 * @brief Get the Random Number Generator Seed as defined in the top level of the XML
		 */
		unsigned int GetSeed();

		/*!
		 * @brief Set the Random Number Generator Seed externally to override any number in the XML
		 */
		void SetSeed( unsigned int new_seed );


		unsigned int GetOriginalSeed() const;

	private:
		/*!
		 * Don't Copy the class this way!
		 */
		MultiXMLConfigReader ( const MultiXMLConfigReader& );

		/*!
		 * Don't Copy the class this way!
		 */
		MultiXMLConfigReader& operator = ( const MultiXMLConfigReader& );

		vector<XMLConfigReader*> XMLReaders;
		int storedSeed;
		int storedRepeats;

		vector<string> _fileNames;
};

#endif


