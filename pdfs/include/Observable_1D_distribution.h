// $Id: Observable_1D_distribution.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Observable_1D_distribution Observable_1D_distribution.h
 *
 *  RapidFit PDF 
 *
 *  @author Pete Clarke
 *  @date 2011-07-30
 */

#ifndef Observable_1D_distribution_H
#define Observable_1D_distribution_H

#include "TH1.h"

#include "ObservableRef.h"
#include "IDataSet.h"
#include "BasePDF.h"

class Observable_1D_distribution : public BasePDF
{
	public:
		Observable_1D_distribution( PDFConfigurator* );
		~Observable_1D_distribution();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		// Physics parameters
		ObservableRef wantedObservable;

		IDataSet* givenDataSet;

		string wanted_poly;

		bool use_function;

		TH1* normalised_histogram;

		TF1* functional_form;
};

#endif

