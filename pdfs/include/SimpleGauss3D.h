
#include "BasePDF.h"

using namespace::std;

class SimpleGauss3D : public BasePDF
{

	public:
		SimpleGauss3D( PDFConfigurator* );
		~SimpleGauss3D();

		void MakePrototypes();

		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

	private:
		//	Observable
		ObservableRef xName;
		ObservableRef yName;
		ObservableRef zName;
		//	PhysicsParameter
		ObservableRef sigmaName;
		ObservableRef sigma2Name;
		ObservableRef sigma3Name;
		//	Internal object(s)
};

