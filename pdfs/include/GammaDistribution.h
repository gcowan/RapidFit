// $Id: GammaDistribution.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class GammaDistribution GammaDistribution.h
 *
 *  RapidFit PDF 
 *
 *  @author Pete Clarke
 *  @date 2011-07-30
 */

#ifndef GammaDistribution_H
#define GammaDistribution_H

#include "BasePDF.h"
#include "ObservableRef.h"

using namespace::std;

class GammaDistribution : public BasePDF
{
	public:
		GammaDistribution( PDFConfigurator* );
		~GammaDistribution();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
	
		//TF1 * PDF ;
	
		// Physics parameters
		double gamma ;
		double mu ;
		double beta ;
		ObservableRef GFgammaName;
		ObservableRef GFmuName;
		ObservableRef GFbetaName;

		// Observables
		double x ;
		ObservableRef GFxName;	
};

#endif

