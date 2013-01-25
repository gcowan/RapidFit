
#pragma once
#ifndef VectoredFeldmanCousins_H
#define VectoredFeldmanCousins_H

#include "IStudy.h"
#include "PDFWithData.h"
#include "FitResultVector.h"

#include <string>
#include <vector>

using namespace::std;

class VectoredFeldmanCousins : public IStudy
{
	public:

		/*!
		 * @brief Constructor Class
		 *
		 * @param input_GlobalResult  This is a 
		 *
		 * @param ResultsForFC  This is the vector of fit results to be used to run over to perform a VectoredFC
		 *
		 * @param inputNuisenceModel  This governs the treatment of the nuisence parameters in the study
		 *
		 * @param new_makeOutput  This is the OutputConfiguration which contains information on the control parameters for scans
		 *
		 * @param newMinimiser  This is the MinimiserConfiguration to be used for performing fits
		 *
		 * @param newFunction  This is the FitFunctionConfiguration to be used to construct fits
		 *
		 * @param new_xmlFile  This is the XML used as a template to construct the fits
		 *
		 * @param new_pdfsAndData  This is the vector of PDFWithData containing the PDFs and Datasets
		 *
		 */
		VectoredFeldmanCousins( FitResultVector* input_GlobalResult, FitResultVector* ResultsForFC, unsigned int inputNuisenceModel,
				OutputConfiguration* new_makeOutput, MinimiserConfiguration* newMinimiser, FitFunctionConfiguration* newFunction,
				XMLConfigReader* new_xmlFile, vector< PDFWithData* > new_pdfsAndData );

		/*!
		 * @brief Destructor Method
		 */
		~VectoredFeldmanCousins();

		/*!
		 * @brief Function for changing the level of Verbosity in the Study
		 */
		void SetOutput( int OutputLevel );

		/*!
		 * @brief Function which launches the FC study itself
		 */
		void DoWholeStudy( int OutputLevel =-1 );

		/*!
		 * @brief Function which returns the Result of the whole study
		 */
		FitResultVector* GetStudyResult();

		/*!
		 * @brief Function to control the number of toys we wish to run at each grid point
		 */
		void SetNumRepeats( int );

		/*!
		 * @brief NULL function to satisfy the IStudy method
		 *
		 * Can be used to pass additional options to the study in a standardised way i.e. to keep/throw bad fits and such-like
		 */
		void SetCommandLineParams( vector<string> );

	private:

		//	Can't be copied
		VectoredFeldmanCousins ( const VectoredFeldmanCousins& );
		VectoredFeldmanCousins& operator=( const VectoredFeldmanCousins& );

		/*!
		 * @brief Initialize internal PDFWithData objects which are used for generating the toys in an FC study
		 */
		void InitalizeNewPDFWithData();

		/*!
		 * @brief used for initializing the free/fixed parameter sets to be used in the current fit
		 */
		void InitializePhysicsParameters( ParameterSet* inputParameters );

		/*!
		 * @brief Used to get the datasets corresponding to the fit for a given toy at a given grid point
		 */
		vector<IDataSet*> GetNewDataSets( ParameterSet* input_params );

		/*!
		 * @brief Used for controlling the Output Level of the study during the fit to ensure the minimal information is output to screen
		 */
		void ResetOutput();

		ParameterSet* getParameterSet( ParameterSet* inputSet, ResultParameterSet* inputResult );

		/*!
		 * Internal Objects specific to the FC study
		 */

		FitResultVector* GlobalFitResult;
		ParameterSet* GlobalFitPhysicsParameters;
		FitResultVector* FitAtGridPoints;
		streambuf *cout_bak, *cerr_bak, *clog_bak;
		vector<PhaseSpaceBoundary*> allPhaseSpaces;
		vector<PDFWithData*> input_pdfsAndData;
		vector<string> controlled_parameters;
		unsigned int nuisenceModel;
		vector<IPDF*> stored_pdfs;
		vector<DataSetConfiguration*> stored_dataconfigs;
		bool sWeighted_study;
		vector<double> sweight_error;
		vector<double> generate_n_events;
		ParameterSet* ParameterSetWithFreeParameters;
		ParameterSet* ParameterSetWithFixedParameters;;
};

#endif

