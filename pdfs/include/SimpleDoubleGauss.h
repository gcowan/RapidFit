
#include "BasePDF.h"

using namespace::std;

class SimpleDoubleGauss : public BasePDF
{

	public:
		SimpleDoubleGauss( PDFConfigurator* );
		~SimpleDoubleGauss();

		void MakePrototypes();

		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

		vector<string> PDFComponents();
		double EvaluateComponent( DataPoint*, ComponentRef* );

	private:
		//	Observable
		ObservableRef xName;
		//	PhysicsParameter
		ObservableRef sigma1Name;
		ObservableRef sigma2Name;
		ObservableRef fracName;
		//	Internal object(s)
		bool plotComponents;
		unsigned int componentIndex;
};

