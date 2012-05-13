
#ifndef CORRECT_COVA_H
#define CORRECT_COVA_H

///	ROOT Headers
#include "TMatrixD.h"
#include "TMatrixDSym.h"
///	RapidFit Headers
#include "IMinimiser.h"
#include "RapidFitMatrix.h"

using namespace::std;

class CorrectedCovariance
{
	public:
		/*!
		 * @brief Method for printing the contents of a Matrix object
		 *
		 * @warning Assumed square matrix atm
		 *
		 * @param Input the TMatrixD object to be printed
		 *
		 * @return Void
		 */
		static void DumpMatrix( TMatrixD* Input );

		/*!
		 * @brief Method for printing the contents of a Matrix object
		 *
		 * @warning Assumes square matrix atm
		 *
		 * @param Input the TMatrixDSym object to be printed
		 *
		 * @return Void
		 */
		static void DumpMatrix( TMatrixDSym* Input );

		/*!
		 * @brief Method to reduce a Covariant Matrix to only contain the relations between truely free fit Parameters
		 *
		 * This method checks the list of free Pararmeters in the fit and removes any parameters which have been externally
		 * constrained and returns the reduced correlation matrix between the remaining parameters
		 *
		 * Fixed Parameters introduce 0 rows/columns and make the Matrix singular
		 *
		 * Externally Constrained parameters give badly defined correlations (several orders of magnitudes out!!) and must be removed to keep results sane
		 *
		 * @param thisMininiser  This is a pointet to the Minimiser which contains information on the Parameters in the fit
		 *
		 * @param Raw_Covariance_Matrix  This is the Matrix which contains parameters which you require to be removed
		 *
		 * @return Returns a RapidFitMatrix which contains the Matrix and Names of the Parameters being related
		 */
		static RapidFitMatrix* GetReducedMatrix( IMinimiser* thisMinimiser, RapidFitMatrix* Raw_Covariance_Matrix );

		/*!
		 * @brief Method which will take the Minmiser once it has minimised and correct the correlation and parameter erros if Weights were used in the fit
		 *
		 * This method uses the same technique employed by RooFit:
		 * http://root.cern.ch/root/html/src/RooAbsPdf.cxx.html#fhMiHE  (11/05/2012)
		 *
		 * 
		 *
		 * @return Returns a RapidFitMatrix which contains the Matrix and Names of the Parameters being related
		 */
		static RapidFitMatrix* GetCorrectedCovarianceMatrix( IMinimiser* thisMinimiser );
};

#endif

