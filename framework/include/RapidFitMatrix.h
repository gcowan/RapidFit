
#pragma once

#ifndef RAPIDFIT_MATRIX_H
#define RAPIDFIT_MATRIX_H

///     ROOT Headers
#include "TMatrixD.h"
#include "TMatrixDSym.h"
///	System Headers
#include <vector>
#include <string>

using namespace::std;

class RapidFitMatrix
{
	public:
		/*!
		 * @brief Default Constructor
		 */
		RapidFitMatrix() : thisMatrix(NULL), theseParameters()
		{}

		/*!
		 * @brief Destructor
		 */
		~RapidFitMatrix()
		{
			if( thisMatrix != NULL ) delete thisMatrix;
		}

		/*!
		 * @brief Copy Constructor
		 */
		RapidFitMatrix ( const RapidFitMatrix& input ) :
			thisMatrix(NULL), theseParameters(input.theseParameters)
		{
			if( input.thisMatrix != NULL )
			{
				//	This has not been checked and may be broken for non-square matricies depending on how you count columns/rows
				thisMatrix = new TMatrixDSym( input.thisMatrix->GetNcols() );
				for( int i=0; i< input.thisMatrix->GetNcols(); ++i )
				{
					for( int j=0; j< input.thisMatrix->GetNrows(); ++j )
					{
						(*thisMatrix)(i,j) = (*(input.thisMatrix))(i,j);
					}
				}
			}
		}

		/*!
		 * @brief a Pointer to the Matrix this class is in control of
		 */
		TMatrixDSym* thisMatrix;

		/*!
		 * @brief a list of Parameter Names which are in the same order as the axis of this matrix (assumes matrix to be square)
		 */
		vector<string> theseParameters;

	private:
		/*!
		 * @brief Don't copy the Class this Way!
		 */
		RapidFitMatrix& operator= ( const RapidFitMatrix& );
};

#endif

