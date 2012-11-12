// $Id: MistagDistribution.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class MistagDistribution MistagDistribution.h
 *
 *  RapidFit PDF 
 *
 *  @author Pete Clarke
 *  @date 2011-07-30
 */

#ifndef MistagDistribution_H
#define MistagDistribution_H

#include "BasePDF.h"
//#include "TF1.h"

class MistagDistribution : public BasePDF
{
	public:
		MistagDistribution( PDFConfigurator* );
		~MistagDistribution();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*,PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);

		// Physics parameters
		double gamma ;
		double mu ;
		double beta ;
		double shoulder ;
		ObservableRef GFgammaName;
		ObservableRef GFmuName;
		ObservableRef GFbetaName;
		ObservableRef GFshoulderName;

		// Observables
		double x ;
		ObservableRef GFxName;	
		ObservableRef tagName;

		bool NumericallyIntegrate;
};

#endif

