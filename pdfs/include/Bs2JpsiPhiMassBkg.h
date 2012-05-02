// $Id: Bs2JpsiPhiMassBkg.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhiMassBkg Bs2JpsiPhiMassBkg.h
 *
 *  RapidFit PDF for Bs2JpsiPhi long lived background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef Bs2JpsiPhiMassBkg_H
#define Bs2JpsiPhiMassBkg_H

#include "BasePDF.h"

class Bs2JpsiPhiMassBkg : public BasePDF
{
	public:
		Bs2JpsiPhiMassBkg( PDFConfigurator* );
		~Bs2JpsiPhiMassBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		ObservableRef alphaM_prName;	// mass decay constant
		// Observables
		ObservableRef recoMassName;	// reconstructed Bs mass
		ObservableRef constraint_recoMassName;
};

#endif

