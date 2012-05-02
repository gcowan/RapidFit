// $Id: CrystalBall.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class CrystalBall CrystalBall.h
 *
 *  RapidFit PDF 
 *
 *  @author Pete Clarke
 *  @date 2011-07-30
 */

#ifndef CrystalBall_H
#define CrystalBall_H

#include "BasePDF.h"

class CrystalBall : public BasePDF
{
	public:
		CrystalBall( PDFConfigurator* );
		~CrystalBall();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		double ApproxErf( double ) const;

		// Physics parameters
		ObservableRef m0Name;	// fraction
		ObservableRef sigmaName;	// width 1
		ObservableRef alphaName;	// width 2 
		ObservableRef nName;	// Bs mass

		// Observables
		ObservableRef recoMassName;	// reconstructed Bs mass
};

#endif

