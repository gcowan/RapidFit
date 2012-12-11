
#include "BasePDF.h"

using namespace::std;

class BoxPDF : public BasePDF
{

	public:
		BoxPDF( PDFConfigurator* );
		~BoxPDF();
		
		void makePrototypes();

		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

		bool SetPhysicsParameters( ParameterSet* );
	private:
		//	Observable
		ObservableRef xName;

		ObservableRef minName, maxName;

		double minimum, maximum, range, norm;

};

