/**
	@class PDFWithData

	A class for creating/storing a PDF and its associated data set

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-5
*/

#ifndef PDF_WITH_DATA_H
#define PDF_WITH_DATA_H

//	RapidFit Headers
#include "IPDF.h"
#include "PhaseSpaceBoundary.h"
#include "IPrecalculator.h"
#include "DataSetConfiguration.h"
//	System Headers
#include <string>
#include <vector>

class PDFWithData
{
	public:
		PDFWithData();
		PDFWithData( IPDF*, PhaseSpaceBoundary*, vector< DataSetConfiguration* >, vector< IPrecalculator* > );
		~PDFWithData();

		void AddCachedData( vector<IDataSet*> );
		bool SetPhysicsParameters( vector<ParameterSet*> );
		IPDF * GetPDF();
		IDataSet * GetDataSet();
		DataSetConfiguration* GetDataSetConfig();

	private:
		//	Uncopyable!
		PDFWithData ( const PDFWithData& );
		PDFWithData& operator = ( const PDFWithData& );

		IPDF * fitPDF;
		PhaseSpaceBoundary * inputBoundary;
		bool parametersAreSet;
		vector< IPrecalculator* > dataProcessors;
		vector< DataSetConfiguration* > dataSetMakers;
		vector< IDataSet* > cached_data;
};

#endif
