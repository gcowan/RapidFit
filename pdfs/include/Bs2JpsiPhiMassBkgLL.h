// $Id: Bs2JpsiPhiMassBkgLL.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhiMassBkgLL Bs2JpsiPhiMassBkgLL.h
 *
 *  RapidFit PDF for Bs2JpsiPhi long lived background
 *
 *  @author Pete Clarke
 *  @date 2011 10 23
 */

#ifndef Bs2JpsiPhiMassBkgLL_H
#define Bs2JpsiPhiMassBkgLL_H

#include "BasePDF.h"

class Bs2JpsiPhiMassBkgLL : public BasePDF
{
	public:
		Bs2JpsiPhiMassBkgLL( PDFConfigurator* );
		~Bs2JpsiPhiMassBkgLL();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		ObservableRef alphaM_llName;	// mass decay constant
		// Observables
		ObservableRef recoMassName;	// reconstructed Bs mass
		ObservableRef constraint_recoMassName;
};

#endif

