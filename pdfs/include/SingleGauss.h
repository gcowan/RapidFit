
#include "BasePDF.h"

using namespace::std;

class SingleGauss : public BasePDF
{

	public:
		SingleGauss( PDFConfigurator* );
		~SingleGauss();

		void MakePrototypes();

		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

		bool SetPhysicsParameters( ParameterSet* );
	private:
		//	Observable
		ObservableRef xName;
		//	PhysicsParameter
		ObservableRef sigmaName;
		ObservableRef centerName;
		//	Internal object(s)

		double centerValue, sigma_den, Norm;
};

