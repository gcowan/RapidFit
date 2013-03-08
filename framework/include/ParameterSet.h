/**
 * @class ParameterSet
 *
 * A collection of physics parameters
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef PARAMETER_SET_H
#define PARAMETER_SET_H

///	ROOT Headers
#include "TString.h"
///	RapidFit Headers
#include "ObservableRef.h"
#include "PhysicsParameter.h"
///	System Headers
#include <vector>
#include <string>

using namespace::std;

class ParameterSet
{
	public:
		/*!
		 * @brief This is the constructor which effectively merges all given ParameterSet objects
		 *
		 * @warning This will take the First occurance of each PhysicsParameter it finds!
		 *
		 * @param Input this is the list of ParameterSets which you intend to merge
		 */
		ParameterSet( vector<ParameterSet*> Input );

		/*!
		 * @brief This is the constructor which creates a set of pointers to PhysicsParameters
		 *
		 * @warning This does NOT check for duplicate PhysicsParameters!!!
		 *
		 * @param Input the list of names you wish to create Parameters for
		 */
		ParameterSet( vector<string> Input );

		/*!
		 * @brief Copy Constructor
		 */
		ParameterSet( const ParameterSet& Input );

		/*!
		 * @brief Assignment Operator
		 */
		ParameterSet& operator = ( const ParameterSet& );

		/*!
		 * @brief Comparison Operator
		 */
		bool operator== ( const ParameterSet& );

		/*!
		 * @brief Destructor
		 */
		~ParameterSet();

		/*!
		 * @brief This gives a list of all known Physics Parameters in this instance
		 *
		 * @return Returns the list of all of the Physics Parameter names in a vector of strings
		 */
		vector<string> GetAllNames() const;

		/*!
		 * @brief This gives a list of all known Free Physics Parameters in this instance
		 *
		 * @return Returns the list of all of the Free Physics Parameter names in a vector of strings
		 */
		vector<string> GetAllFloatNames() const;

		/*!
		 * @brief This gives a list of all known Fixed Physics Parameters in this instance
		 *
		 * @warning This is technically everything NOT FREE (but not free is fixed isn't it?)
		 *
		 * @return Returns the list of all of the Fixed Physics Parameter names in a vector of strings
		 */
		vector<string> GetAllFixedNames() const;

		/*!
		 * @brief This returns the Parameter Requested by the Name
		 *
		 * @param Name   This contains the Name of the Parameter being Requested
		 *
		 * @return returns the Parameter with the requested Name
		 */
		PhysicsParameter * GetPhysicsParameter( const string ) const;

		/*!
		 * @brief This returns the Parameter Requested by the Name and stores the location
		 *
		 * @param Input    This contains the Name of the Requested Parameter and will store the position of the Parameter once it's been found for future reference
		 *
		 * @return returns the Parameter with the requested Name, or at the given location if it's already been found once before and it's location stored
		 */
		PhysicsParameter * GetPhysicsParameter( const ObservableRef& ) const;

		/*!
		 * @brief This changes the requested PhysicsParameter to be the same as the given Parameter
		 *
		 * @param Name    This is the name of the PhysicsParameter wished to be changed
		 *
		 * @param Input   This is the input Parameter which should be used to change the existing Parameter
		 *
		 * @return true, parameter exists and was changed, false this doesn't exist and/or wasn't changed
		 *         This is not checked could this be made void?
		 */
		bool SetPhysicsParameter( string, PhysicsParameter* );

		/*!
		 * @brief This changes the requested PhysicsParameter to be the same as the given arguments
		 *
		 * @param Name      This is the name of the PhysicsParameter wished to be changed
		 *
		 * @param Value     This is the value the PhyscisParameter should now take
		 *
		 * @param Minimum   This is the new Minimum of the PhysicsParameter
		 *
		 * @param Maximum   This is the new Maximum of the PhysicsParameter
		 *
		 * @param Step      This is the new Step Size of the PhysicsParameter
		 *
		 * @param NewName   This is the new Name that the PhysicsParameter should take   ONLY CHANGED IN THE PhysicsParameter OBJECT!!!
		 *
		 * @param NewUnit   This is the new Unit that the PhysicsParameter should take
		 *
		 * @return true, parameter exists and was changed, false this doesn't exist and/or wasn't changed
		 *         This is not checked could this be made void?
		 */
		bool SetPhysicsParameter( string Name, double Value, double Minimum, double Maximum, double Step, string NewName, string NewUnit );

		/*!
		 * @brief This sets the PhysicsParameters to be the same as the provided ParameterSet
		 *
		 * @param Input    This is a pointer to the ParameterSet to be taken as the new Values, This is not altered by this function
		 *
		 * @return true = everything went as expected , false an error occurred (This is 'never' tested should this just be a void routine?)
		 */
		bool SetPhysicsParameters( const ParameterSet* Input );

		/*!
		 * @brief This is used from the various FitFunctions to translate new input values into PhysicsParameter Values
		 *
		 * @warning Not very nice in OO programming terms, and unsafe. Much faster though
		 *
		 * @param Input    This is a pointer to an array of values which is implicitly assumed to be the same size as the number of variables in this ParameterSet
		 *
		 * @return true = everything went as expected , false an error occurred (This is 'never' tested should this just be a void routine?)
		 */
		void UpdatePhysicsParameters( double* Input, int npar =-1 );

		/*!
		 * @brief This returns the Values of all of the Physics Parameters in a flat 1D double array
		 *
		 * @return This returns a pointer to a 1D array which contains the values of all of the PhysicsParameters in the
		 */
		double* GetPhysicsParameters() const;

		/*!
		 * @brief This is used from the various FitFunctions to translate new input values into PhysicsParameter Values
		 *
		 * @param Input This is the vector of new values for the Parameters
		 *
		 * @return true = everything went as expected , false an error occurred (This is 'never' tested should this just be a void routine?)
		 */
		void SetPhysicsParameters( vector<double> Input );

		/*!
		 * @brief This Adds a New PhysicsParameter which is the same as the PhysicsParameter Provided
		 *
		 * @warning This creates a new or modifies an exiting internal object and does not look after the input Parameter
		 *
		 * @param InputParam    This is the PhysicsParameter wished to be added to this ParameterSet
		 *
		 * @param replace       Should the internal values be modified? true = yes, false = no   default(no)
		 *
		 * @return true = everything went as expected , false an error occurred (This is 'never' tested should this just be a void routine?)
		 */
		bool AddPhysicsParameter( const PhysicsParameter* InputParam, bool replace =false );

		/*!
		 * @brief This Adds the given ParameterSet to this set
		 *
		 * @warning This creates a new or modifies an exiting internal object and does not look after the input ParameterSet
		 *
		 * @param InputParam    This is the ParameterSet wished to be added to this ParameterSet
		 *
		 * @param replace       Should the internal values be modified? true = yes, false = no   default(no)
		 *
		 * @return true = everything went as expected , false an error occurred (This is 'never' tested should this just be a void routine?)
		 */
		bool AddPhysicsParameters( const ParameterSet* InputParam, bool replace =false );

		/*!
		 * @brief Output some debugging info
		 *
		 * @return Void
		 */
		void Print() const;

		/*!
		 * @brief Return the required XML for this ParameterSet
		 *
		 * @return Returns the Full ParameterSet XML section of the XML
		 */
		string XML() const;

		/*!
		 * @brief 
		 *
		 * @return
		 */
		size_t GetUniqueID() const;

		/*!
		 * @brief
		 *
		 * @return void
		 */
		void FloatedFirst();

	private:

		void SetUniqueID( size_t );

		mutable vector<PhysicsParameter*> allParameters;	/*!	vector of pointers to all of the Physics Parameters managed by this ParameterSet			*/
		vector<string> allNames;			/*!	vector of strings of the names of all of the Physics Parameters						*/
		mutable size_t uniqueID;

		mutable vector<ObservableRef> allInternalNames;
		mutable vector<ObservableRef> allForeignNames;
};

#endif

