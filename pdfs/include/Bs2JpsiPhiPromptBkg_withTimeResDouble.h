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
		Bs2JpsiPhiPromptBkg_withTimeResDouble();
		~Bs2JpsiPhiPromptBkg_withTimeResDouble();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		string frac_sigmaPrName; // sigma of the prompt dist
		string sigmaPrName; // sigma of the prompt dist
		string sigmaPr2Name; // sigma of the prompt dist
		// These contain the strings that correspond
		// to the observable names that are used in the PDF
		string timeName;    // proper time
};

#endif
