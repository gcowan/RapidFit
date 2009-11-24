// $Id: RaPDF_Bs2JpsiPhiLongLivedBkg.h,v 1.3 2009/11/11 18:30:09 bwynne Exp $
/** @class RaPDF_Bs2JpsiPhiLongLivedBkg RaPDF_Bs2JpsiPhiLongLivedBkg.h
 *
 *  RapidFit PDF for Bs2JpsiPhi long lived background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef RaPDF_Bs2JpsiPhiLongLivedBkg_H
#define RaPDF_Bs2JpsiPhiLongLivedBkg_H

#include "BasePDF.h"

class RaPDF_Bs2JpsiPhiLongLivedBkg : public BasePDF
{
	public:
		RaPDF_Bs2JpsiPhiLongLivedBkg();
		~RaPDF_Bs2JpsiPhiLongLivedBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		string tau1Name;		// decay constant 1
		string tau2Name;		// decay constant 2
		string f_LL1Name;		// fraction

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		string timeName;	// proper time
};

#endif
