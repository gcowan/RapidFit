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

                double EvaluateComponent( DataPoint*, ComponentRef* );

                vector<string> PDFComponents();

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
	
		//TF1 * PDF ;
	
		// Physics parameters
		double sigma1 ;
		double theta1 ;
		double m1 ;
		double sigma2 ;
		double theta2 ;
		double m2 ;
		double f ;
		ObservableRef LNsigma1Name;
		ObservableRef LNtheta1Name;
		ObservableRef LNm1Name;
		ObservableRef LNsigma2Name;
		ObservableRef LNtheta2Name;
		ObservableRef LNm2Name;
		ObservableRef LNfName;

		// Observables
		double x ;
		ObservableRef LNxName;	

		bool plotComponents;
};

#endif

