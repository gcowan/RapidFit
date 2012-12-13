
#include "BasePDF.h"

using namespace::std;

class PolyPDF : public BasePDF
{

	public:
		PolyPDF( PDFConfigurator* );
		~PolyPDF();

		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

		bool SetPhysicsParameters( ParameterSet* );
	private:

		//	Observable
		ObservableRef xName;

		vector<ObservableRef> parameterNames;
		vector<double> parameterValues;

		unsigned int order;
};

