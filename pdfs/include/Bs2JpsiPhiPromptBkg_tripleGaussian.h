/** @class Bs2JpsiPhiPromptBkg_tripleGaussian Bs2JpsiPhiPromptBkg_tripleGaussian.h
 *
 *  PDF for Bs2JpsiPhi prompt Jpsi background with time resolution - double gaussian
 *
 *  @author Pete Clarke
 *  @date 2011-01-22
 */

#ifndef Bs2JpsiPhiPromptBkg_tripleGaussian_H
#define Bs2JpsiPhiPromptBkg_tripleGaussian_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif

class Bs2JpsiPhiPromptBkg_tripleGaussian : public BasePDF
{
	public:
		Bs2JpsiPhiPromptBkg_tripleGaussian();
		~Bs2JpsiPhiPromptBkg_tripleGaussian();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		string frac_sigmaPr1Name;    // Fraction of first Gaussian : second two Gaussians
		string frac_sigmaPr23Name;   // Fraction of second : third Gaussian
		string sigmaPr1Name;			 // Width of Gaussian 1
		string sigmaPr2Name;		 // Width of Gaussian 2
		string sigmaPr3Name;		 // Width of Gaussian 3

		// These contain the strings that correspond to the observable names that are used in the PDF
		string timeName;    // proper time
};

#endif
