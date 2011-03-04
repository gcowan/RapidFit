/**
	@class PDFWithData

	A class for creating/storing a PDF and its associated data set

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-5
*/

#ifndef PDF_WITH_DATA_H
#define PDF_WITH_DATA_H

#include "IPDF.h"
#include "PhaseSpaceBoundary.h"
#include "IPrecalculator.h"
#include "DataSetConfiguration.h"
#include <string>
#include <vector>

class PDFWithData
{
	public:
		PDFWithData();
		PDFWithData( IPDF*, PhaseSpaceBoundary*, vector< DataSetConfiguration* >, vector< IPrecalculator* > );
		~PDFWithData();

		void AddCachedData( vector<IDataSet*> );
		bool SetPhysicsParameters( ParameterSet* );
		IPDF * GetPDF();
		IDataSet * GetDataSet();
		DataSetConfiguration* GetDataSetConfig();

	private:
		IPDF * fitPDF;
		PhaseSpaceBoundary * inputBoundary;
		bool parametersAreSet;
		vector< IPrecalculator* > dataProcessors;
		vector< DataSetConfiguration* > dataSetMakers;
		vector<IDataSet*> cached_data;
};

#endif
