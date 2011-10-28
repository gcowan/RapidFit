/** @class Bs2JpsiPhiPromptBkg_withTimeResDouble Bs2JpsiPhiPromptBkg_withTimeResDouble.h
 *
 *  PDF for Bs2JpsiPhi prompt Jpsi background with time resolution - double gaussian
 *
 *  @author Pete Clarke
 *  @date 2019-10-28
 */

#ifndef Bs2JpsiPhiPromptBkg_withTimeResDouble_H
#define Bs2JpsiPhiPromptBkg_withTimeResDouble_H

#include "BasePDF.h"

class Bs2JpsiPhiPromptBkg_withTimeResDouble : public BasePDF
{
	public:
		Bs2JpsiPhiPromptBkg_withTimeResDouble( PDFConfigurator* );
		~Bs2JpsiPhiPromptBkg_withTimeResDouble();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		ObservableRef frac_sigmaPrName; // sigma of the prompt dist
		ObservableRef sigmaPrName; // sigma of the prompt dist
		ObservableRef sigmaPr2Name; // sigma of the prompt dist
		// These contain the ObservableRefs that correspond
		// to the observable names that are used in the PDF
		ObservableRef timeName;    // proper time
		ObservableRef timeconstraintName;
};

#endif
