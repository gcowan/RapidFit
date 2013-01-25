
#include "BasePDF.h"

using namespace::std;

class LandauGauss : public BasePDF
{

	public:
		LandauGauss( PDFConfigurator* );
		~LandauGauss();

		void MakePrototypes();

		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

	private:
		//	Observable
		ObservableRef xName;
		//	PhysicsParameter
		ObservableRef landauSigmaName;
		ObservableRef gaussSigmaName;
		ObservableRef mpvName;
		//	Internal object(s)
};

