/** @class Bs2JpsiPhiMassSignal Bs2DsPiMassSignal.h
 *
 *  RapidFit PDF for Bs2DsPi mass signal
 *
 *  @author Gemma Fardell gfardell@cern.ch
 *  @date 2010-01-29
 */

#ifndef Bs2DsPiMassSignal_H
#define Bs2DsPiMassSignal_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif

class Bs2DsPiMassSignal : public BasePDF
{
	public:
		Bs2DsPiMassSignal();
		~Bs2DsPiMassSignal();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		ObservableRef f_sig_m1Name;	// fraction
		ObservableRef sigma_m1Name;	// width 1
		ObservableRef sigma_m2Name;	// width 2 
		ObservableRef m_BsName;	// Bs mass

		// Observables
		ObservableRef recoMassName;	// reconstructed Bs mass
};

#endif
