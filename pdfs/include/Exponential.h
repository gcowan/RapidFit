// $Id: Exponential.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Exponential Exponential.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution.
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef Exponential_H
#define Exponential_H

#include "BasePDF.h"
#include "IResolutionModel.h"

class Exponential : public BasePDF
{
	public:
		Exponential( PDFConfigurator* );
		Exponential( const Exponential & );
		~Exponential();

		//Calculate the PDF value
		double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		Exponential operator = ( const Exponential& );

		bool SetPhysicsParameters(ParameterSet*);

		// Physics parameters
		ObservableRef tauName;		// decay constant 1
		ObservableRef timeName;		// proper time
		ObservableRef timeConst;

		double tau;
		double gamma;

		// These contain the ObservableRefs that correspond
		// to the observable names that are used in the
		// PDF.
		double time;

		IResolutionModel* resolutionModel;

		void MakePrototypes();
};

#endif
