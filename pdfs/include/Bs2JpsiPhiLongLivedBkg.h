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
		Bs2JpsiPhiLongLivedBkg();
		~Bs2JpsiPhiLongLivedBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		pair<string,int> tau1Name;		// decay constant 1
		pair<string,int> tau2Name;		// decay constant 2
		pair<string,int> f_LL1Name;		// fraction

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		pair<string,int> timeName;	// proper time
		pair<string,int> constraint_timeName;
};

#endif
