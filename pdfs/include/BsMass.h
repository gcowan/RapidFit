// $Id: BsMass.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class BsMass BsMass.h
 *
 *  RapidFit PDF 
 *
 *  @author Pete Clarke
 *  @date 2011-07-30
 */

#ifndef BsMass_H
#define BsMass_H

#include "BasePDF.h"

class BsMass : public BasePDF
{
	public:
		BsMass( PDFConfigurator* );
		~BsMass();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

		double EvaluateComponent( DataPoint*, ComponentRef* );

		vector<string> PDFComponents();
	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		ObservableRef f_sig_m1Name;	// fraction
		ObservableRef sigma_m1Name;	// width 1
		ObservableRef ratio_21Name;	// width 2 
		ObservableRef m_BsName;	// Bs mass

		// Observables
		ObservableRef recoMassName;	// reconstructed Bs mass
	
		//Limits
		double mlow, mhigh;

		int componentIndex;

		bool plotComponents;
};

#endif

