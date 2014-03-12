
#ifndef RAPIDFIT_XMLObjectGenerator_H
#define RAPIDFIT_XMLObjectGenerator_H

#include "XMLTag.h"
#include "PhysicsParameter.h"
#include "PhaseSpaceBoundary.h"
#include "IConstraint.h"
#include "ParameterSet.h"
#include "FitFunctionConfiguration.h"
#include "IPDF.h"
#include "DataSetConfiguration.h"
#include "MinimiserConfiguration.h"
#include "OutputConfiguration.h"
#include "ExternalConstraint.h"
#include "ExternalConstMatrix.h"
#include "ConstraintFunction.h"
#include "ScanParam.h"
#include "ComponentPlotter_config.h"
#include "PDFWithData.h"

#include <string>
#include <vector>

using namespace::std;

class XMLObjectGenerator
{
	public:

		static PhysicsParameter * GetPhysicsParameter( XMLTag*, string& );

		static PhaseSpaceBoundary* GetPhaseSpaceBoundary( XMLTag* InputTag );

		static IConstraint * GetConstraint( XMLTag*, string& Name );

		static ParameterSet * GetParameterSet( XMLTag* );

		static FitFunctionConfiguration * MakeFitFunction( XMLTag* );

		static IPDF * GetNamedPDF( XMLTag*, PhaseSpaceBoundary*, XMLTag*, XMLTag*, ParameterSet* thisParameterSet, bool print=true );

		static IPDF * GetPDF( XMLTag*, PhaseSpaceBoundary*, XMLTag* overloadConfigurator, XMLTag* common, ParameterSet* thisParameterSet, bool print=true );

		static DataSetConfiguration * MakeDataSetConfiguration( XMLTag*, PhaseSpaceBoundary*, XMLTag* common, ParameterSet* thisParameterSet, int seed );

		static MinimiserConfiguration * MakeMinimiser( XMLTag*, OutputConfiguration* );

		static ExternalConstraint * GetExternalConstraint( XMLTag* );

		static ExternalConstMatrix * GetExternalConstMatrix( XMLTag * InputTag );

		static ConstraintFunction * GetConstraintFunction( XMLTag* );

		static ScanParam * GetScanParam( XMLTag * InputTag );

		static pair<ScanParam*, ScanParam*> Get2DScanParam( XMLTag * InputTag );

		static CompPlotter_config* getCompPlotterConfigs( XMLTag* CompTag );

		static pair< string, string > MakeContourPlot( XMLTag* );

		static OutputConfiguration * MakeOutputConfiguration( XMLTag* );

		static PDFWithData * GetPDFWithData( XMLTag*, XMLTag*, int StartVal, XMLTag* overloadConfigurator, XMLTag* common, ParameterSet* thisParameterSet, int seed, PhaseSpaceBoundary* thisPhaseSpaceBoundary );

	private:
		XMLObjectGenerator();
		~XMLObjectGenerator();
};

#endif

