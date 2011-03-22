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

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif

class Bs2JpsiPhiMassBkg : public BasePDF
{
	public:
		Bs2JpsiPhiMassBkg();
		~Bs2JpsiPhiMassBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		pair<string,int> alphaM_prName;	// mass decay constant
		// Observables
		pair<string,int> recoMassName;	// reconstructed Bs mass
		pair<string,int> constraint_recoMassName;
};

#endif
