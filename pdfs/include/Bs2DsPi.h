/** @class Bs2DsPi Bs2DsPi.h
 *
 *  RapidFit PDF for Bs2DsPi
 *
 *  @author Gemma Fardell
 */

#ifndef Bs2DsPi_H
#define Bs2DsPi_H

#include "BasePDF.h"

class Bs2DsPi : public BasePDF
{
	public:
		Bs2DsPi( PDFConfigurator* );
		~Bs2DsPi();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);
		virtual bool SetPhysicsParameters(ParameterSet*);
		virtual vector<string> GetDoNotIntegrateList();

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// These contain the strings that correspond
		// to the physics parameter names that will be
		// used in the minimiser.
		ObservableRef gammaName;		// gamma
		ObservableRef deltaGammaName;	// delta gamma
		ObservableRef deltaMName;		// delta mass

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		ObservableRef timeName;		// proper time
		ObservableRef tagName;		// B tag
		ObservableRef mistagName;		// B mistag

		ObservableRef timeconstraintName;

		void getPhysicsParameters( double&, double&, double&);
		
		void getTimeDependentFuncs( double&, double&, double&, DataPoint*);
		void getTimeDependentFuncsInt(double&, double&, double&, PhaseSpaceBoundary*);

};

#endif
