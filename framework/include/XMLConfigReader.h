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
#include "FitFunctionConfiguration.h"
#include "MinimiserConfiguration.h"
#include "OutputConfiguration.h"
#include "DataSetConfiguration.h"
#include "ConstraintFunction.h"

using namespace std;

class XMLConfigReader
{
	public:
		XMLConfigReader();
		XMLConfigReader(string);
		~XMLConfigReader();

		ParameterSet * GetFitParameters();
		MinimiserConfiguration * GetMinimiserConfiguration();
		FitFunctionConfiguration * GetFitFunctionConfiguration();
		OutputConfiguration * GetOutputConfiguration();
		vector< PDFWithData* > GetPDFsAndData();
		vector< ConstraintFunction* > GetConstraints();
		int GetNumberRepeats();
		bool IsLoaded();

		int GetSeed();	//  Return the Random Seed
		void SetSeed( int new_seed );	//  Set Seed returned by XMLFile

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
		FitFunctionConfiguration * MakeFitFunction( XMLTag* );
		MinimiserConfiguration * MakeMinimiser( XMLTag* );
		pair< string, string > MakeContourPlot( XMLTag* );
		OutputConfiguration * MakeOutputConfiguration( XMLTag* );
		IPrecalculator * MakePrecalculator( XMLTag*, PhaseSpaceBoundary* );
		DataSetConfiguration * MakeDataSetConfiguration( XMLTag*, PhaseSpaceBoundary* );
		ConstraintFunction * GetConstraintFunction( XMLTag* );
		ExternalConstraint * GetExternalConstraint( XMLTag* );

		vector< XMLTag* > children;
		bool isLoaded;

		vector<int> seed;	//  Random Seed
};

#endif
