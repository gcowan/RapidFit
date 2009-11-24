/**
	@class PDFWithData

	A class for creating/storing a PDF and its associated data set

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-5
*/

#ifndef PDF_WITH_DATA_H
#define PDF_WITH_DATA_H

#include "IPDF.h"
#include "IDataSet.h"
#include "PhaseSpaceBoundary.h"
#include <string>
#include <vector>

class PDFWithData
{
	public:
		PDFWithData();
		PDFWithData( IPDF*, string, long, vector<string>, PhaseSpaceBoundary* );
		PDFWithData( IPDF*, IPDF*, string, long, vector<string>, PhaseSpaceBoundary* );
		~PDFWithData();

		bool SetPhysicsParameters( ParameterSet* );
		IPDF * GetPDF();
		IDataSet * GetDataSet();

	private:
		IPDF * fitPDF;
		IPDF * generatePDF;
		PhaseSpaceBoundary * inputBoundary;
		string dataSource;
		vector<string> dataArguments;
		long dataAmount;
		bool parametersAreSet;
};

#endif
