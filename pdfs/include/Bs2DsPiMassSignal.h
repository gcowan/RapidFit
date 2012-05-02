/** @class Bs2JpsiPhiMassSignal Bs2DsPiMassSignal.h
 *
 *  RapidFit PDF for Bs2DsPi mass signal
 *
 *  @author Gemma Fardell gfardell@cern.ch
 *  @date 2010-01-29
 */

#ifndef Bs2DsPiMassSignal_H
#define Bs2DsPiMassSignal_H

#include "BasePDF.h"

class Bs2DsPiMassSignal : public BasePDF
{
	public:
		Bs2DsPiMassSignal( PDFConfigurator* );
		~Bs2DsPiMassSignal();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		ObservableRef f_sig_m1Name;	// fraction
		ObservableRef sigma_m1Name;	// width 1
		ObservableRef alpha_m1Name;	// width 1
		ObservableRef eta_m1Name;	// width 1
		ObservableRef sigma_m2Name;	// width 2 
		ObservableRef alpha_m2Name;	// width 1
		ObservableRef eta_m2Name;	// width 1	
		ObservableRef m_BsName;	// Bs mass

		// Observables
		ObservableRef recoMassName;	// reconstructed Bs mass
		ObservableRef constraint_recoMassName;
	
};

#endif

