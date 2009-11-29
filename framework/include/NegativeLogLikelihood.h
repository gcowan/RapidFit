/**
        @class NegativeLogLikelihood

        A fit function with evaulate methods for an NLL calculation

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef NEGATIVE_LOG_LIKELIHOOD_H
#define NEGATIVE_LOG_LIKELIHOOD_H

#include "FitFunction.h"

class NegativeLogLikelihood : public FitFunction
{
	public:
		NegativeLogLikelihood();
		NegativeLogLikelihood(string);
		~NegativeLogLikelihood();

		virtual double UpErrorValue(int);

	protected:
		virtual double EvaluateDataSet( IPDF*, IDataSet*, RapidFitIntegrator*, int );
		virtual double EvaluateParameterSet( ParameterSet*, vector<string> );
};
#endif
