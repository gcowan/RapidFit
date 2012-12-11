
#include "BasePDF.h"

using namespace::std;

class FlatPDF : public BasePDF
{

	public:
		FlatPDF( PDFConfigurator* );
		~FlatPDF();

		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

		bool SetPhysicsParameters( ParameterSet* );
	private:
		//	Observable
		ObservableRef xName;
};

