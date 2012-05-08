// $Id: Bs2JpsiPhiMassSignal.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhiMassSignal Bs2JpsiPhiMassSignal.h
 *
 *  RapidFit PDF for Bs2JpsiPhi long lived background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef Bs2JpsiPhiMassSignal_H
#define Bs2JpsiPhiMassSignal_H

#include "BasePDF.h"

class Bs2JpsiPhiMassSignal : public BasePDF
{
	public:
		Bs2JpsiPhiMassSignal( PDFConfigurator* );
		~Bs2JpsiPhiMassSignal();

		//Calculate the PDF value
		double Evaluate(DataPoint*);

		vector<string> PDFComponents();
		double EvaluateComponent( DataPoint*, ComponentRef* );

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		ObservableRef f_sig_m1Name;	// fraction
		ObservableRef sigma_m1Name;	// width 1
		ObservableRef sigma_m2Name;	// width 2 
		ObservableRef m_BsName;	// Bs mass

		// Observables
		ObservableRef recoMassName;	// reconstructed Bs mass

		int componentIndex;
		bool plotComponents;
};

#endif

