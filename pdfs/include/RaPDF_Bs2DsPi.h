/** @class RaPDF_Bs2DsPi RaPDF_Bs2DsPi.h
 *
 *  RapidFit PDF for Bs2DsPi
 *
 *  @author Gemma Fardell
 */

#ifndef RaPDF_Bs2DsPi_H
#define RaPDF_Bs2DsPi_H

#include "BasePDF.h"

class RaPDF_Bs2DsPi : public BasePDF
{
	public:
		RaPDF_Bs2DsPi();
		~RaPDF_Bs2DsPi();

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
		string gammaName;		// gamma
		string deltaGammaName;	// delta gamma
		string deltaMName;		// delta mass

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		string timeName;		// proper time
		string tagName;		// B tag
		string mistagName;		// B mistag
	
		void getPhysicsParameters( double&, double&, double&);
		
		void getTimeDependentFuncs( double&, double&, double&, DataPoint*);
		void getTimeDependentFuncsInt(double&, double&, double&, PhaseSpaceBoundary*);

};

#endif
