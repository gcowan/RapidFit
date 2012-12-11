
//	BasePDF where everything is that this PDF inherits how to interface
#include "BasePDF.h"
//	The Definition of the ObservableRef object
#include "ObservableRef.h"
//	The Definition of the ComponentRef object
#include "ComponentRef.h"
//	The Definition of the DataPoint object
#include "DataPoint.h"
//	The Definition of the ParameterSet object
#include "ParameterSet.h"
//	The Definition of the PhaseSpaceBoundary object
#include "PhaseSpaceBoundary.h"

using namespace::std;

//	This PDF must inherit from the BasePDF to be a PDF in a RapidFit Context
class OptimisedDoubleGauss : public BasePDF
{

	public:
		//	Required C++ Objects
		OptimisedDoubleGauss( PDFConfigurator* );
		~OptimisedDoubleGauss();

		//	Copy Constructor.
		//	Not strictly Required in this case but here to demonstrate how one should be used
		//	Only EXPLICITLY required when you have object references within your private member variables
		OptimisedDoubleGauss( const OptimisedDoubleGauss& );

		//	Method to initialize Observables/ParameterSet
		void MakePrototypes();

		//	Methods where the main Calculations are performed
		double Evaluate( DataPoint* );
		double Normalisation( PhaseSpaceBoundary* );

		//	Method to advertise and Evaluate the Components of this PDF
		vector<string> PDFComponents();
		double EvaluateComponent( DataPoint*, ComponentRef* );

		//	Method to intercept the SetPhysicsParameters and cache some basic calculations
		bool SetPhysicsParameters( ParameterSet* );

	private:
		//	Copy by Assignment is NOT performed within RapidFit.
		//	All PDFs are expected be copied through ClassLookUp::CopyPDF and this does NOT make use of assignment operators
		//	In Principle Assignment Operators are fine, although,
		//	for derrived types they require A LOT of work/effort to write/maintain. Hence why we do NOT use them!
		OptimisedDoubleGauss& operator=( const OptimisedDoubleGauss& );


		//	Each DataPoint contains the Method GetObservable( string ) and GetObservable( ObservableRef )
		//	Using the 2nd method and constructing ObservableRef objects means that the DataPoint does NOT need to 
		//	Perform a lookup about where in the DataPoint the object is for every single event in the nTuple
		//	This replaces 1,000,000 calls to StringProcessing::VectorContains with 1,000 calls a clear improvement

		//	Observable
		ObservableRef xName;
		//	PhysicsParameter
		ObservableRef sigma1Name;
		ObservableRef sigma2Name;
		ObservableRef fracName;
		ObservableRef centerName;

		//	Internal object(s)

		//	Component Plotting
		bool plotComponents;
		unsigned int componentIndex;

		//	Storing the Result of Calculations
		double sigma1_denom, sigma2_denom, f, f2, denominator, center;
};

