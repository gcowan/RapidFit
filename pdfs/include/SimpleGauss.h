
#include "BasePDF.h"

using namespace::std;

class SimpleGauss : public BasePDF
{

	public:
		SimpleGauss( PDFConfigurator* );
		~SimpleGauss();

		void MakePrototypes();

		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

	private:
		//	Observable
		ObservableRef xName;
		//	PhysicsParameter
		ObservableRef sigmaName;
		//	Internal object(s)
};

