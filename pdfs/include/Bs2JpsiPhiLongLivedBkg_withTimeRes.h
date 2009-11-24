// $Id: Bs2JpsiPhiLongLivedBkg_withTimeRes.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Bs2JpsiPhiLongLivedBkg_withTimeRes Bs2JpsiPhiLongLivedBkg_withTimeRes.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution.
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef Bs2JpsiPhiLongLivedBkg_withTimeRes_H
#define Bs2JpsiPhiLongLivedBkg_withTimeRes_H

#include "BasePDF.h"

class Bs2JpsiPhiLongLivedBkg_withTimeRes : public BasePDF
{
	public:
		Bs2JpsiPhiLongLivedBkg_withTimeRes();
		~Bs2JpsiPhiLongLivedBkg_withTimeRes();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		double erfc( double, double, double);
		double erfcInt( double, double, double);
		
		// Physics parameters
		string tau1Name;		// decay constant 1
		string tau2Name;		// decay constant 2
		string f_LL1Name;		// fraction
		string sigmaLL1Name;		// time res sigma 1
		string sigmaLL2Name;		// time res sigma 2

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		string timeName;	// proper time
};

#endif
