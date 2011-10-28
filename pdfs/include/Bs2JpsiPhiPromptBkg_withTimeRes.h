// $Id: Bs2JpsiPhiPromptBkg_withTimeRes.h,v 1.1 2009/11/13 09:57:06 gcowan Exp $
/** @class Bs2JpsiPhiPromptBkg_withTimeRes Bs2JpsiPhiPromptBkg_withTimeRes.h
 *
 *  PDF for Bs2JpsiPhi prompt Jpsi background with time resolution
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#ifndef Bs2JpsiPhiPromptBkg_withTimeRes_H
#define Bs2JpsiPhiPromptBkg_withTimeRes_H

#include "BasePDF.h"

class Bs2JpsiPhiPromptBkg_withTimeRes : public BasePDF
{
	public:
		Bs2JpsiPhiPromptBkg_withTimeRes( PDFConfigurator* );
		~Bs2JpsiPhiPromptBkg_withTimeRes();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

                // Physics parameters
                ObservableRef sigmaPrName; // sigma of the prompt dist

		// These contain the strings that correspond
		// to the observable names that are used in the PDF
		ObservableRef timeName;    // proper time

		ObservableRef timeconstraintName;
};

#endif
