// $Id: Bs2DsPiBkg_withTimeRes.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Bs2DsPiBkg_withTimeRes Bs2DsPiBkg_withTimeRes.h
 *
 *  PDF for Bs2DsPi  background from Bd->Dpi
 *
 *  @author Greig A Cowan, Pete Clarke
 *  @date 2009-10-04
 */

#ifndef Bs2DsPiBkg_withTimeRes_H
#define Bs2DsPiBkg_withTimeRes_H

#include "BasePDF.h"

class Bs2DsPiBkg_withTimeRes : public BasePDF
{
	public:
		Bs2DsPiBkg_withTimeRes( PDFConfigurator* );
		~Bs2DsPiBkg_withTimeRes();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:

		void MakePrototypes();
		
		// Physics parameters
		ObservableRef lifetimeBdName;		// decay constant 1
		ObservableRef timeResName ;		// decay constant 2

		// Observables
		ObservableRef timeName;	// proper time
		ObservableRef timeconstraintName;
};

#endif
