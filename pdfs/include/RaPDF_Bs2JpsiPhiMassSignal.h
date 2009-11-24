// $Id: RaPDF_Bs2JpsiPhiMassSignal.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class RaPDF_Bs2JpsiPhiMassSignal RaPDF_Bs2JpsiPhiMassSignal.h
 *
 *  RapidFit PDF for Bs2JpsiPhi long lived background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef RaPDF_Bs2JpsiPhiMassSignal_H
#define RaPDF_Bs2JpsiPhiMassSignal_H

#include "BasePDF.h"

class RaPDF_Bs2JpsiPhiMassSignal : public BasePDF
{
	public:
		RaPDF_Bs2JpsiPhiMassSignal();
		~RaPDF_Bs2JpsiPhiMassSignal();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		string f_sig_m1Name;	// fraction
		string sigma_m1Name;	// width 1
		string sigma_m2Name;	// width 2 
		string m_BsName;	// Bs mass

		// Observables
		string recoMassName;	// reconstructed Bs mass
};

#endif
