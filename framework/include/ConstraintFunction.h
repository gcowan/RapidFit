/*!
 *  @class ConstraintFunction
 *
 *  @brief Where external, experimental constraints on PhysicsParameters are calculated
 * 
 *  @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef CONSTRAINT_FUNCTION_H
#define CONSTRAINT_FUNCTION_H

///	RapidFit Headers
#include "ParameterSet.h"
#include "ExternalConstraint.h"
///	System Headers
#include <vector>

using namespace::std;

class ConstraintFunction
{
	public:
		/*!
		 * @brief Constructor which loops over all ExternalConstraints in the XML
		 *
		 * @param Input   This is a list of the ExternalConstraints of parameters in the ParameterSet around their minima
		 */
		ConstraintFunction( const vector<ExternalConstraint*> Input );

		/*!
		 * @brief Copy Constructor
		 */
		ConstraintFunction( const ConstraintFunction& );

		/*!
		 * @brief Detructor
		 */
		~ConstraintFunction();

		/*!
		 * @brief Evaluate this Constraint Function(s) for this ParameterSet
		 *
		 * @oaram Input   This is a pointer to the ParameterSet being tested
		 *
		 * @return This returns a correction number which is +ve and is the effect of where parameters fall in these constraints
		 */
		double Evaluate( const ParameterSet* Input );

		/*!
		 * @brief Output some debugging info
		 *
		 * @return Void
		 */
		void Print() const;

		/*!
		 * @brief Generate the XML required to generate this object
		 *
		 * @return Return the XML for this in string Format
		 */
		string XML() const;

		/*!
		 * @brief
		 *
		 * @return
		 */
		vector<string> ConstrainedParameter() const;

	private:

		/*!
		 * Don't Copy the class this way!
		 */
		ConstraintFunction& operator= ( const ConstraintFunction& );


		vector< int > Found_Position;				/*	! Undocumented!	*/
		vector< ExternalConstraint* > allConstraints;		/*	Internal list of all the External Constraints requested	*/
};

#endif

