
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

		VectoredFeldmanCousins( FitResultVector* input_GlobalResult, FitResultVector* ResultsForFC, unsigned int inputNuisenceModel, OutputConfiguration* new_makeOutput, MinimiserConfiguration* newMinimiser, FitFunctionConfiguration* newFunction, XMLConfigReader* new_xmlFile, vector< PDFWithData* > new_pdfsAndData );

		~VectoredFeldmanCousins();

		void SetOutput( int OutputLevel );

		void DoWholeStudy( int OutputLevel =-1 );

		FitResultVector* GetStudyResult();

		void SetNumRepeats( int );

		void SetCommandLineParams( vector<string> );

	private:
		//	Can't be copied
		VectoredFeldmanCousins ( const VectoredFeldmanCousins& );
		VectoredFeldmanCousins& operator=( const VectoredFeldmanCousins& );

		void InitalizeNewPDFWithData();

		void InitializePhysicsParameters( ParameterSet* inputParameters );

		vector<IDataSet*> GetNewDataSets( ParameterSet* input_params );

		void ResetOutput();

		ParameterSet* getParameterSet( ParameterSet* inputSet, ResultParameterSet* inputResult );

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

