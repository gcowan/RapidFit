/**
        @class ToyStudy

        A class that automates a whole toy study

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef TOY_STUDY_H
#define TOY_STUDY_H

//	RapidFit Headers
#include "IStudy.h"
#include "PDFWithData.h"
#include "ParameterSet.h"
#include "FitResultVector.h"
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class ToyStudy	:	public IStudy
{
	public:
		/*!
		 * @brief Constructor
		 *
		 * @param inputMinimiser This is the InputMinimiser to be used to perform the minimisations
		 *
		 * @param inputFunction  This is the FitFunctionConfiguration to be used to construct the fits
		 *
		 * @param inputPDFWData  This is the vector of input PDFs and Datasets
		 *
		 * @param input_Const  This is a vector of external constraints to the fit
		 *
		 * @param num  Number of Repeats to perform in the study
		 */
		ToyStudy( MinimiserConfiguration* inputMinimiser, FitFunctionConfiguration* inputFunction, ParameterSet* inputParam,
				vector< PDFWithData* > inputPDFWData, vector< ConstraintFunction* > input_Const, int num );

		/*!
		 * @brief Destructor
		 */
		~ToyStudy();

		/*!
		 * @brief Function to call the Toy Study to start
		 *
		 * @param output This is the output verbosity level
		 */
		void DoWholeStudy( int output = -999 );

		/*!
		 * @brief Function to return the output from the whole Study
		 *
		 * @return This returns a pointer to the internal FitResultVector
		 *
		 * @warning The returned object is controlled and destoryed by this class
		 */
		FitResultVector* GetStudyResult();

		/*!
		 * @brief Function to Change the number of repeats that the Toy Study will perform
		 */
		void SetNumRepeats( int );

		/*!
		 * @brief Function to pass additional Command Line Parameters to the Fit
		 */
		void SetCommandLineParams( vector<string> );

		/*!
		 * @brief Calling this will stop the Toy Study attempting to obtain the required good number of fit results
		 *
		 * This is useful if you suspect the Study may fall over and perform infinite repeats
		 */
		void SetFixedNumberToys();

		/*!
		 * @brief Calling this will alow the output from all (including failed) toys to be saved
		 */
		void setSaveAllToys();

	private:
		//	Uncopyable!
		ToyStudy ( const ToyStudy& );
		ToyStudy& operator = ( const ToyStudy& );

		FitResult * GenerateAndMinimise();

		bool fixedNumToys;
		bool saveAllToys;
};

#endif

