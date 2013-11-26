/*!
 * @class ExternalConstMatrix
 *
 * @brief A class that holds experimentally dervied constraints on fit parameters
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef EXTERNAL_CONSTMATRIX_H
#define EXTERNAL_CONSTMATRIX_H

#include "IConstraintFunction.h"
//	System Headers
#include <string>

using namespace::std;

class ExternalConstMatrix : public IConstraintFunction
{
	public:
		/*!
		 * @brief Constructor containing the external constraint name, value and error
		 *
		 * @param Name  Name of the Physics Parameter we're applying the External Constrint to
		 *
		 * @param Value This is CV of the external constraint being applied
		 *
		 * @param Error This is the Error on the Constraint 'width of the Gaussian'
		 */
		ExternalConstMatrix( string Name, string Value, string Error, string Correlations );

		ExternalConstMatrix( const ExternalConstMatrix& input );

		/*!
		 * Destructor
		 */
		~ExternalConstMatrix();

		/*!
		 * @brief Get the Name of the ExternalConstMatrix
		 *
		 * @return Returns the name of the PhysicsParameter Being Constrained
		 */
		string GetName() const;

		void SetPhysicsParameters( const ParameterSet* input );

		bool CanApply( const ParameterSet* input ) const;

		double GetChi2() const;

		void Print() const;

		string XML() const;

		bool isExternalConst() const;

		bool isExternalMatrix() const;
	private:
		string names;		/*!	External Constraint Name		*/
		string values, errors;	/*!	External Constraint Value and Error	 */
		string correlations;

		vector<double> names_val, values_val, errors_val;
		vector<vector<double> > corr_matrix;

		ParameterSet* internalParameterSet;
		vector<string> wantedParameters;
};

#endif

