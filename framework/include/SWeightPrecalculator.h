/**
  @class SWeightPrecalculator

  A Precalculator for producing the sWeights for a data set for backgorund subtraction

  @author Benjamin Wynne bwynne@cern.ch
  @date 2009-12-14
  */

#ifndef S_WEIGHT_PRECALCULATOR_H
#define S_WEIGHT_PRECALCULATOR_H

//	RapidFit Headers
#include "IPrecalculator.h"
#include "IPDF.h"
//	System Headers
#include <vector>
#include <string>

class SWeightPrecalculator : public IPrecalculator
{
	public:
		SWeightPrecalculator();
		SWeightPrecalculator( IPDF*, IPDF*, ParameterSet*, string WeightName = "sWeight" );
		~SWeightPrecalculator();

		virtual IDataSet * ProcessDataSet( IDataSet* );
		pair< double, double > CalculateMatrixElements( long, long, IDataSet*, vector<double>&, vector<double>& );

	private:
		//	Uncopyable!
		SWeightPrecalculator ( const SWeightPrecalculator& );
		SWeightPrecalculator& operator = ( const SWeightPrecalculator& );
		IPDF * signalPDF;
		IPDF * backgroundPDF;
		ParameterSet * fitParameters;
		string weightName;
};

#endif
