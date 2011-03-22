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

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif

class Bs2JpsiPhiMassSignal : public BasePDF
{
	public:
		Bs2JpsiPhiMassSignal();
		~Bs2JpsiPhiMassSignal();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		pair<string,int> f_sig_m1Name;	// fraction
		pair<string,int> sigma_m1Name;	// width 1
		pair<string,int> sigma_m2Name;	// width 2 
		pair<string,int> m_BsName;	// Bs mass

		// Observables
		pair<string,int> recoMassName;	// reconstructed Bs mass
};

#endif
