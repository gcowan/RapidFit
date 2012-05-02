/** @class Bs2DsPi_mistagParameter Bs2DsPi_mistagParameter.h
 *
 *  RapidFit PDF for Bs2DsPi
 *
 *  Modified by Pete to make mistag a fit parameter and add time resolution
 *
 *  @author Gemma Fardell
 */

#ifndef Bs2DsPi_mistagParameter_H
#define Bs2DsPi_mistagParameter_H

#include "BasePDF.h"

class Bs2DsPi_mistagParameter : public BasePDF
{
	public:
		Bs2DsPi_mistagParameter( PDFConfigurator* );
		~Bs2DsPi_mistagParameter();

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
		ObservableRef mistagName;		// B mistag
	       	ObservableRef timeresName;		// B mistag

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		ObservableRef timeName;		// proper time
		ObservableRef tagName;		// B tag

		ObservableRef timeconstraintName;

		// Physics parameters
	    double gamma, deltaGamma, deltaM, mistag, timeRes ;
	
	    // Observables
		double time ;
		int tag;

	    // Limits
		double tlow, thigh ;
	
		void getPhysicsParameters(  );
	    void getObservables( DataPoint* );
		
		void getTimeDependentFuncsInt(double&, double&, double&, PhaseSpaceBoundary*);
	
		// Widths
		double gamma_l() const ;
		double gamma_h() const ;
		double gambar() const ;
	
		// Time primitives
		double expL() const ;
		double expH() const ;
		double expCos() const ;
		double expLint() const ;
		double expHint() const ;
		double expCosInt() const ;
		
		// Functions to help convolve a single gaussian into time primitives
		// DIDNT APPEAR TO BE USED RooComplex evalCerf( double, double, double ) const ;
		//RooComplex evalCerfApprox( double, double, double ) const ;
		//double evalCerfRe( double, double, double ) const ;
		//double evalCerfIm( double, double, double ) const ;
	

};

#endif
