/** @class Sinusoid Sinusoid.h
 *
 *  RapidFit Sinusoid 
 *
 *  @author Rob
 *  @date 20xx-yy-zz
 */

#ifndef Sinusoid_H
#define Sinusoid_H

#include "BasePDF.h"

class Sinusoid : public BasePDF
{
	public:
		Sinusoid( PDFConfigurator* );
		~Sinusoid();

		//	Calculate the PDF value
		virtual double Evaluate(DataPoint*);

		double EvaluateComponent( DataPoint*, ComponentRef* );

		vector<string> PDFComponents();

	protected:
		//	Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

		virtual bool SetPhysicsParameters( ParameterSet* );
	private:
		void MakePrototypes();
		void MakePrototypeDataPoint();
		void MakePrototypeParameterSet();

		//	Physics parameters
		ObservableRef PhaseName, FreqName, AmplitudeName;

		//	Observables
		ObservableRef TimeName;

		double A, f, phase;

		bool useCos;

		bool plotComponents;

		//	Which Component am I evalueating right now
		int componentIndex;

};

#endif

