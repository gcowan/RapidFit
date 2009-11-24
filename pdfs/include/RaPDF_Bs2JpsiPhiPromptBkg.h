// $Id: RaPDF_Bs2JpsiPhiPromptBkg.h,v 1.2 2009/11/13 09:57:06 gcowan Exp $
/** @class RaPDF_Bs2JpsiPhiPromptBkg RaPDF_Bs2JpsiPhiPromptBkg.h
 *
 *  RapidFit PDF for Bs2JpsiPhi prompt Jpsi background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef RaPDF_Bs2JpsiPhiPromptBkg_H
#define RaPDF_Bs2JpsiPhiPromptBkg_H

#include "BasePDF.h"

class RaPDF_Bs2JpsiPhiPromptBkg : public BasePDF
{
	public:
		RaPDF_Bs2JpsiPhiPromptBkg();
		~RaPDF_Bs2JpsiPhiPromptBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

                // Physics parameters
                string meanName;       	// mean of the prompt dist
                string sigmaName;    	// sigma of the prompt dist

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		string timeName;	// proper time
		string tagName;		// B tag
		string mistagName;	// B mistag
};

#endif
