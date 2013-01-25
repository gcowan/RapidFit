/**
  @class ResultParameterSet

  A set of physics parameters after fitting

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

#pragma once
#ifndef RESULT_PARAMETER_SET_H
#define RESULT_PARAMETER_SET_H

//	RapidFit Headers
#include "ResultParameter.h"
#include "ParameterSet.h"
#include "RapidFitMatrix.h"
//	System Headers
#include <vector>
#include <string>

using namespace::std;

class ResultParameterSet
{
	public:
		/*!
		 * @brief Constructor
		 */
		ResultParameterSet( vector<string> );

		/*!
		 * @brief Copy Constructor
		 */
		ResultParameterSet( const ResultParameterSet& );

		/*!
		 * @brief Destructor
		 */
		~ResultParameterSet();

		/*!
		 * @brief Get the list of all ResultParameter names in this ResultParameterSet
		 *
		 * @return returns a vector of strings containing all ResultParameter names
		 */
		vector<string> GetAllNames() const;

		/*!
		 * @brief Get list of Floated/Free Parameters in this ResultParameterSet
		 *
		 * @return returns a vector of strings containing all floated ResultParameter names
		 */
		vector<string> GetAllFloatNames() const;

		/*!
		 * @brief Get list of Fixed Parameters in this ResultParameterSet
		 *
		 * @return returns a vector of strings containing all fixed ResultParameter names
		 */
		vector<string> GetAllFixedNames() const;

		/*!
		 * @brief Get the result parameter according to this index
		 *
		 * @param index This is the absolute index of the parameter wanted
		 *
		 * @return returns a pointer to the ResultParameter object inside this ResultParameterSet
		 */
		ResultParameter * GetResultParameter( int index ) const;

                /*!
                 * @brief Get the result parameter according to this name
                 *
                 * @warning This method always performs an expensive lookup without attempting to cache results
		 *
		 * @param name This is the name of the parameter wanted
		 *
                 * @return returns a pointer to the ResultParameter object inside this ResultParameterSet
                 */
		ResultParameter * GetResultParameter( string name ) const;

                /*!
                 * @brief Get the result parameter according to this ObservableRef
                 *
                 * @param name This is the ObservableRef containing the name of the object wanted
                 *
                 * @return returns a pointer to the ResultParameter object inside this ResultParameterSet
                 */
		ResultParameter * GetResultParameter( const ObservableRef& ) const;

		/*!
		 * @brief Set the result parameter to be equal to input
		 *
		 * @param name This is the name of the parameter requested to be set
		 *
		 * @param templateParam This is a ResultParameter that should be used to change the internal object to match
		 *
		 * @return true if successful, false if not
		 */
		bool SetResultParameter( string name, ResultParameter* templateParam );

		/*!
		 * @brief Set the result parameter to be equal to input
		 *
		 * @param name
		 *
		 * @param origvalue
		 *
		 * @param value
		 *
		 * @prarm error
		 *
		 * @param min
		 *
		 * @param max
		 *
		 * @param type
		 *
		 * @param unit
		 *
		 * @return true if successful, false if not
		 */
		bool SetResultParameter( string name, double value, double origvalue, double error, double min, double max, string type, string unit );

		/*!
		 * @brief Forcibly change or add a new parameter to this object
		 *
		 * @param templateParam
		 *
		 * @return true if successful, false if not
		 */
		bool ForceNewResultParameter( ResultParameter* templateParam );

		/*!
		 * @brief Forcibly change or add a new parameter to this input
		 *
		 * @return true if successful, false if not
		 */
		bool ForceNewResultParameter( string name, double value, double origvalue, double error, double min, double max, string type, string unit );

		/*!
		 * @brief Returns a ParameterSet with the same values and errors from this ResultParameterSet with correct Parameters fixed/free
		 *
		 * @return Returns a new ParameterSet with same properties as this ResultParameterSet
		 *
		 * @warning The returned object is not maintained by this class (clean up after yourself)
		 */
		ParameterSet* GetDummyParameterSet() const;

		/*!
		 * @brief Return the XML required to reconstruct this class
		 *
		 * @return Returns the XML in a 'flat' string
		 */
		void Print() const;

		/*!
		 * @brief Return the XML required to run some fit output/post-processing (i.e. all Parameters fixed)
		 *
		 * @return Returns the XML in a 'flat' string
		 */
		string FitXML() const;

		/*!
		 * @brief Return the XML required to run a set of toys with this ParameterSet (i.e. fixed/free)
		 *
		 * @return Returns the XML in a 'flat' string
		 */
		string ToyXML() const;

                /*!
                 * @brief Apply a covariance matrix to this ResultParameterSet to replace the stored errors
                 *
                 * @return void
                 */
		void ApplyCovarianceMatrix( RapidFitMatrix* Input );
	private:
		ResultParameterSet operator=( const ResultParameterSet& input );

		vector<ResultParameter*> allParameters;
		vector<string> allNames;

		/*!
		 * @brief Internal method for returning the XML required to regenerate this object
		 *
		 * @param input Should Free parameters be Free?
		 *
		 * @return Returns the XML in a 'flat' string
		 */
		string XML( const bool input=true ) const;
};

#endif

