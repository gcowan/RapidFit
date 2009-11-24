// $Id: RaPDF_Bs2JpsiPhiMassBkg.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class RaPDF_Bs2JpsiPhiMassBkg RaPDF_Bs2JpsiPhiMassBkg.h
 *
 *  RapidFit PDF for Bs2JpsiPhi long lived background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef RaPDF_Bs2JpsiPhiMassBkg_H
#define RaPDF_Bs2JpsiPhiMassBkg_H

#include "BasePDF.h"

class RaPDF_Bs2JpsiPhiMassBkg : public BasePDF
{
	public:
		RaPDF_Bs2JpsiPhiMassBkg();
		~RaPDF_Bs2JpsiPhiMassBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		string alphaM_prName;	// mass decay constant
		// Observables
		string recoMassName;	// reconstructed Bs mass
};

#endif
