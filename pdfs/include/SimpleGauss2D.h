
#include "BasePDF.h"

using namespace::std;

class SimpleGauss2D : public BasePDF
{

	public:
		SimpleGauss2D( PDFConfigurator* );
		~SimpleGauss2D();

		void MakePrototypes();

		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

	private:
		//	Observable
		ObservableRef xName;
		ObservableRef yName;
		//	PhysicsParameter
		ObservableRef sigmaName;
		ObservableRef sigma2Name;
		//	Internal object(s)
};

