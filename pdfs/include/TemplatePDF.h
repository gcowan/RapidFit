/** @class TemplatePDF TemplatePDF.h
 *
 *  RapidFit TemplatePDF 
 *
 *  @author Rob
 *  @date 20xx-yy-zz
 */

#ifndef TemplatePDF_H
#define TemplatePDF_H

#include "BasePDF.h"

class TemplatePDF : public BasePDF
{
	public:
		TemplatePDF( PDFConfigurator* );
		~TemplatePDF();

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
		ObservableRef ParameterName;

		//	Observables
		ObservableRef ObservableName;

		bool plotComponents;

		//	Which Component am I evalueating right now
		int componentIndex;

};

#endif

