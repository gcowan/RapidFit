// $Id: Bs2JpsiPhiLongLivedBkg.h,v 1.3 2009/11/11 18:30:09 bwynne Exp $
/** @class Bs2JpsiPhiLongLivedBkg Bs2JpsiPhiLongLivedBkg.h
 *
 *  RapidFit PDF for Bs2JpsiPhi long lived background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef Bs2JpsiPhiLongLivedBkg_H
#define Bs2JpsiPhiLongLivedBkg_H

#include "BasePDF.h"

class Bs2JpsiPhiLongLivedBkg : public BasePDF
{
	public:
		Bs2JpsiPhiLongLivedBkg( PDFConfigurator* );
		~Bs2JpsiPhiLongLivedBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		ObservableRef tau1Name;		// decay constant 1
		ObservableRef tau2Name;		// decay constant 2
		ObservableRef f_LL1Name;		// fraction

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		ObservableRef timeName;	// proper time
		ObservableRef constraint_timeName;
};

#endif
