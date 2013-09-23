/*!
 * @class PhysicsBottle
 *
 *  A collection of PDF-DataSet pairs to be fitted simultaneously with a given parameter set.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef PHYSICS_BOTTLE_H
#define PHYSICS_BOTTLE_H

//	RapidFit Headers
#include "IPDF.h"
#include "IDataSet.h"
#include "ParameterSet.h"
#include "ConstraintFunction.h"
//	System Headers
#include <vector>

class PhysicsBottle
{
	public:
		/*!
		 * @brief Constructor to construct the PhysicsBottle
		 *
		 * @param InputSet    This is the ParameterSet which is used by the PhysicsBottle
		 */
		PhysicsBottle( const ParameterSet* );

		/*!
		 * @brief Destructor
		 */
		~PhysicsBottle();

		/*!
		 * @brief Copy Constructor
		 */
		PhysicsBottle( const PhysicsBottle& );

		vector<IPDF*> GetAllPDFs() const;

		vector<IDataSet*> GetAllDataSets() const;

		void AddResult( const IPDF*, IDataSet* );

		void AddConstraint( const ConstraintFunction* );

		int NumberResults() const;

		IPDF* GetResultPDF( const int ) const;

		IDataSet* GetResultDataSet( const int ) const;

		vector< ConstraintFunction* > GetConstraints() const;

		ParameterSet * GetParameterSet() const;

		void SetParameterSet( const ParameterSet* );

		void Print() const;

	private:
		/*!
		 * Don't Copy the class this way!
		 */
		PhysicsBottle& operator = ( const PhysicsBottle& );

		vector< IPDF* > allPDFs;
		vector< IDataSet* > allDataSets;
		vector< ConstraintFunction* > allConstraints;
		ParameterSet * bottleParameters;
};

#endif

