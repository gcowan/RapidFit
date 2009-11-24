/**
        @class XMLConfigReader

        Opens an xml config file and uses it to create RapidFit data objects

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef XML_CONFIG_READER_H
#define XML_CONFIG_READER_H

#include <vector>
#include <string>
#include "XMLTag.h"
#include "ParameterSet.h"
#include "PDFWithData.h"

using namespace std;

class XMLConfigReader
{
	public:
		XMLConfigReader();
		XMLConfigReader(string);
		~XMLConfigReader();

		ParameterSet * GetFitParameters();
		string GetMinimiserName();
		string GetFitFunctionName();
		vector< PDFWithData* > GetPDFsAndData();
		int GetNumberRepeats();
		bool IsLoaded();

	private:
		ParameterSet * GetParameterSet( XMLTag* );
		PhysicsParameter * GetPhysicsParameter( XMLTag*, string& );
		PhaseSpaceBoundary * GetPhaseSpaceBoundary( XMLTag* );
		IConstraint * GetConstraint( XMLTag*, string& );
		PDFWithData * GetPDFWithData( XMLTag*, XMLTag* );
		IPDF * GetNamedPDF( XMLTag* );
		IPDF * GetSumPDF( XMLTag*, PhaseSpaceBoundary* );
		IPDF * GetNormalisedSumPDF( XMLTag*, PhaseSpaceBoundary* );
		IPDF * GetProdPDF( XMLTag*, PhaseSpaceBoundary* );
		IPDF * GetPDF( XMLTag*, PhaseSpaceBoundary* );

		vector< XMLTag* > children;
		bool isLoaded;
};

#endif
