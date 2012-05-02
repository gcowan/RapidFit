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
		PhysicsBottle( ParameterSet* );

		/*!
		 * @brief Destructor
		 */
		~PhysicsBottle();

		/*!
		 * @brief Copy Constructor
		 */
		PhysicsBottle(const PhysicsBottle&);


		void AddResult( IPDF*, IDataSet* );

		void AddConstraint( ConstraintFunction* );

		int NumberResults();

		IPDF* GetResultPDF(int);

		IDataSet* GetResultDataSet(int);

		vector< ConstraintFunction* > GetConstraints();

		ParameterSet * GetParameterSet();

		bool SetParameterSet( ParameterSet* );

	private:
		/*!
		 * Don't Copy the class this way!
		 */
		PhysicsBottle& operator = ( const PhysicsBottle& );

		vector< IPDF* > allPDFs;
		vector< IDataSet* > allDataSets;
		vector< ConstraintFunction* > allConstraints;
		ParameterSet * bottleParameters;
		bool finalised;
};

#endif
