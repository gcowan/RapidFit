
#include "BasePDF.h"

using namespace::std;

class StudentT : public BasePDF
{

	public:
		StudentT( PDFConfigurator* );
		~StudentT();

		void MakePrototypes();

		double Evaluate( DataPoint* );

		bool SetPhysicsParameters( ParameterSet* );
	private:
		//	Observable
		ObservableRef massName;
		//	PhysicsParameter
		ObservableRef muName;
		ObservableRef sName;
		ObservableRef nName;
		//	Internal object(s)

		double muValue, sValue, nValue;
};

