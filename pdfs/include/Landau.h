
#include "BasePDF.h"

using namespace::std;

class Landau : public BasePDF
{

	public:
		Landau( PDFConfigurator* );
		~Landau();

		void MakePrototypes();

		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

	private:
		//	Observable
		ObservableRef xName;
		//	PhysicsParameter
		ObservableRef sigmaName;
		ObservableRef mpvName;
		ObservableRef CName;
		//	Internal object(s)
};

