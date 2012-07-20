// $Id: LogNormalDistribution.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class LogNormalDistribution LogNormalDistribution.h
 *
 *  RapidFit PDF 
 *
 *  @author Pete Clarke
 *  @date 2011-07-30
 */

#ifndef LogNormalDistribution_H
#define LogNormalDistribution_H

#include "BasePDF.h"
#include "ObservableRef.h"

using namespace::std;

class LogNormalDistribution : public BasePDF
{
	public:
		LogNormalDistribution( PDFConfigurator* );
		~LogNormalDistribution();

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
		double sigma ;
		double theta ;
		double m ;
		ObservableRef LNsigmaName;
		ObservableRef LNthetaName;
		ObservableRef LNmName;

		// Observables
		double x ;
		ObservableRef LNxName;	
};

#endif

