/*!
 * @ingroup Configurators  This Generator class constricts DataSets which can be cast as IDataSets
 * @class DataSetConfiguration
 *
 * @brief A class for holding the data to create a data set
 *
 * @author Benjamin Wynne
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef DATA_SET_CONFIGURATION_H
#define DATA_SET_CONFIGURATION_H

///	ROOT Headers
#include "TNtuple.h"
///	RapidFit Headers
#include "IPDF.h"
#include "IDataSet.h"
#include "FitResultVector.h"
///	System Headers
#include <string>
#include <vector>

using namespace::std;

class DataSetConfiguration
{
	public:
		/*!
		 * @brief Constructor for a file based DataSet
		 *
		 * @param DataSource    This is the name of the DataSource we wish to use, generally 'File' unless something more complex has been coded up
		 *
		 * @param DataNumber    This is the maximum number of DataPoints to request, (ROOT files can NOT deliver more datapoints than they have)
		 *
		 * @param cut           This is a cut to be applied to ROOT DataSets which allows for final cuts to be make as part of the fitting process
		 *
		 * @param DataArguments    This is a list of arguments which have been passed to the DataSetConfiguraton class i.e. somefile.root
		 *
		 * @param DataArgumentNames  This is a list of the corresponding names to each of the arguments which are passed to the constructor i.e. Filename
		 *
		 * @param starting_entry   This is the starting entry that you should start reading a file from. BEWARE this is applied BEFORE THE CUT
		 *
		 * @param InputBoundary  This is an optional PhaseSpaceBoundary defined at Construction
		 */
		DataSetConfiguration( string DataSource, long DataNumber, string cut, vector<string> DataArguments, vector<string> DataArgumentNames, int starting_entry = 0, PhaseSpaceBoundary* InputBoundary = NULL );

		/*!
		 * @brief Constructor for a toy based DataSet
		 *
		 * @param DataSource    This is the name of the DataSource for this Toy, generally Foam
		 *
		 * @param DataNumber    This is the Number of DataPoints we wish to create from this source
		 *
		 * @param cut           This is a cut to be applied on the DataSet, should be applied after it has been constructed.
		 *                      ATM this does NOT serve any function, should this be removed from the interface?
		 *
		 * @param DataArguments    This is a list of arguments which have been passed to the DataSetConfiguraton class i.e. somefile.root                   
                 *  
                 * @param DataArgumentNames  This is a list of the corresponding names to each of the arguments which are passed to the constructor i.e. Filename
                 *
                 * @param DataPDF       This is the PDF to be used for generating the new DataSet
                 *
                 * @param InputBoundary  This is an optional PhaseSpaceBoundary defined at Construction
		 */
		DataSetConfiguration( string DataSource, long DataNumber, string cut, vector<string> DataArguments, vector<string> DataArgumentNames, IPDF * DataPDF, PhaseSpaceBoundary* InputBoundary = NULL );

		DataSetConfiguration ( const DataSetConfiguration& );

		/*!
		 * @brief Destructor
		 */
		~DataSetConfiguration();

		/*!
		 * @brief Set the Physics Parameters for use in construction of a toy DataSet
		 *
		 * @param Input   This is a Pointer to the ParameterSet we wish to use as truth for generating a Toy DataSet
		 *
		 * @return  true if successful, false is not
		 */
		bool SetPhysicsParameters( ParameterSet* Input );

		/*!
		 * @brief Set the source of this DataSet
		 *
		 * @param Input   New DataSource we wish to use, i.e. now we have done a fit to data, switch to Foam for a toy Study
		 *
		 * @return
		 */
		bool SetSource( string Input );

		/*!
		 * @Get the source of this DataSet
		 *
		 * This is the name of the source of the DataSet i.e, Foam/File/...
		 */
		string GetSource() const;

		/*!
		 * @brief External Interface to get a DataSet, either File or Toy Based
		 *
		 * This will actually perform the step of aquiring the DataSet from this Configuraton class
		 *
		 * @param InputBoundary   This is the PhaseSpace that is populated by a DataSet
		 *
		 * @param InputPDF        This is the PDF used in generating a Toy DataSet
		 *
		 * @param number_events   This (if defined) will override the requested DataSet size
		 *
		 * @return Returns a pointer to the Newly created DataSet   THIS CLASS DOES NOT LOOK AFTER CLEANING UP THE DATASET
		 */
		IDataSet* MakeDataSet( PhaseSpaceBoundary* InputBoundary, IPDF* InputPDF, int number_events =-1 );

		/*!
		 * @brief Get a pointer to the PDF used in Generation a Toy DataSet
		 *
		 * @return This simply returns a pointer to the PDF used in the Generation of the DataSet
		 */
		IPDF* GetGenerationPDF() const;

		/*!
		 * @brief Method for turning the debugging output on/off
		 *
		 * @param Input   true will turn debugging on, false will leave it disabled
		 *
		 * @return Void
		 */
		void SetDebug( bool Input );

		/*!
		 * @brief Return the XML required to reconstruct this class
		 *
		 * @return Returns the XML in a 'flat' string
		 */
		string XML() const;

		void SetDebug( DebugClass* input_debug );

		PhaseSpaceBoundary* GetPhaseSpace() const;

		int GetDataSetSize() const;

	private:
		/*!
		 * @brief Don't Copy the class this way!
		 */
		DataSetConfiguration& operator = ( const DataSetConfiguration& );

		/*!
		 * @brief This is the Function which Constructs a Toy DataSet (generally Foam) then stores the output in a .root file
		 *        Then passes it through the File Interpretor to make use of more advanced ROOT features
		 *
		 * @param arguments     These are the arguments passed to the DataSetGenerator from the XML
		 *
		 * @param ArgumentNames These are the names of the arguments passed to the DataSetGenerator
		 *
		 * @param internalBoundary  This is the PhaseSpaceBoundary which should be occupied with the required Boundary
		 *
		 * @param numberEvents  This is the number of events requested to be constructed within the PhaseSpaceBoundary
		 *
		 * @param FitPDF        This is the PDF that should be used to construct the DataPoints within the PhaseSpaceBoundary
		 *
		 * @return Pointer to the new DataSet which has just been constructed. The Generator does not handle destroying of the DataSet
		 */
		IDataSet* FoamFile( vector<string> arguments, vector<string> ArgumentNames, PhaseSpaceBoundary* internalBoundary, int numberEvents, IPDF* FitPDF );


		/*!
		 * @brief This is the Function which Loads a Everything required to construct a Toy DataSet and returns the Constructed DataSet
		 *
		 * @param source           This is the Name of the DataSetGenerator that should be loaded
		 *
		 * @param internalBoundary 
		 *
		 * @param numberEvents     This is the number of events requested to be constructed within the PhaseSpaceBoundary
		 *
		 * @param FitPDF           This is the PDF that should be used to construct the DataPoints within the PhaseSpaceBoundary
		 *
		 * @return Returns the DataSet which has just been constructed
		 */
		IDataSet* LoadGeneratorDataset( string source, PhaseSpaceBoundary* internalBoundary, int numberEvents, IPDF* FitPDF );

		/*!
		 * Private functions for reading objects and assembling DataSet
		 */

		IDataSet * LoadDataFile( vector<string>, vector<string>, PhaseSpaceBoundary*, long );	/*! @brief Undocumented	*/
		IDataSet * LoadAsciiFileIntoMemory( string, long, PhaseSpaceBoundary* );		/*! @brief Undocumented	*/
		IDataSet * LoadRootFileIntoMemory( string, string, long, PhaseSpaceBoundary* );		/*! @brief Undocumented	*/

		/*!
		 * @brief Private method for polling a ROOT file for the ntuple path
		 *
		 * This will return the path of the first found Ntuple in a file
		 *
		 * @param Input   This will be where we look first for a Ntuple in a file, if one is NOT found, then checks are made to find the first sensible NTuple
		 *
		 * @return This will return a sensible path to the first Ntuple Found if the Input is incorrect, returns Input by default
		 */
		string getNtuplePath( string );

		/*!
		 * @brief This is awesome :)
		 *
		 * This function takes your current path in the ROOT universe, searches recursively within this 'directory' to find all classes which inherit explicitly from inherit_type
		 * Then it populates the found_names object with relative paths and names of objects it discovered
		 *
		 * Basically this means that the NTuple path in RapidFit is not an optional requirement for files containing more than one dataset :D
		 *
		 * @param current_path   This is the path of where you currently are. Could be extracted internally to thif function from gDirectory, but this creates... issues
		 *                       Simply provide the internal path from gDirectory when you make the call
		 *
		 * @param found_names    This is a pointer to a vector of pairs which will be populated with the names and relative paths of the discovered objects
		 *
		 * @param inherit_type   This is the class you wish for all objects you care about to inherit from
		 *
		 * @param relative_path  This is a pointer to a TString of the current path. This is explicitly changed depending on where in the search you are
		 *
		 * @return Void  This returns all the information by manipulating the input objects!
		 */
		void get_object_list( TString current_path, vector<pair<string,string> > *found_names, TClass* inherit_type, TString *relative_path );

		string source;			/*!	Name of the source type	*/
		string cutString;		/*!	The CutString used when reading in from a ROOT file	*/
		long numberEvents;		/*!	Number of events to read in from an input file		*/

		vector<string> arguments, argumentNames;

		IPDF * generatePDF;		/*!	PDF used in generation of a toy DataSet			*/

		bool separateGeneratePDF, parametersAreSet;
		int Start_Entry;		/*!	The start entry in the file that the dataset should start looking	*/

		bool DEBUG_DATA;		/*!	Useful flag for turning on/off some information when debugging the dataset	*/

		PhaseSpaceBoundary* internalBoundary;	/*!	Internal pointer to the PhaseSpaceBoundary that corresponds to the DataSet last created or first one to be if one does not already exist */

		IDataSet* internalRef;		/*!	This is the internal reference to the DataSet that has just been created. It is NOT to be destroyed here	*/

		DebugClass* debug;

		string fileName;
};

#endif

