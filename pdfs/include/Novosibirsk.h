// $Id: Novosibirsk.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Novosibirsk Novosibirsk.h
 *
 *  RapidFit PDF 
 *
 *  @author Pete Clarke
 *  @date 2011-07-30
 */

#ifndef Novosibirsk_H
#define Novosibirsk_H

#include "BasePDF.h"

class Novosibirsk : public BasePDF
{
	public:
		Novosibirsk( PDFConfigurator* );
		Novosibirsk( const Novosibirsk & );
		~Novosibirsk();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		ObservableRef widthName;
		ObservableRef peakName;
		ObservableRef tailName;

		double width;
		double peak;
		double tail;

		// Observables
		ObservableRef xName;
};

#endif

