/**
  @class SWeightPrecalculator

  A Precalculator for producing the sWeights for a data set for backgorund subtraction

  @author Benjamin Wynne bwynne@cern.ch
  @date 2009-12-14
  */

#pragma once
#ifndef S_WEIGHT_PRECALCULATOR_H
#define S_WEIGHT_PRECALCULATOR_H

//	RapidFit Headers
#include "IPrecalculator.h"
#include "IPDF.h"
#include "FitResult.h"
//	System Headers
#include <vector>
#include <string>

class SWeightPrecalculator : public IPrecalculator
{
	public:
		SWeightPrecalculator( FitResult* inputResult, string WeightName = "sWeight", unsigned int config=1 );
		~SWeightPrecalculator();

		virtual IDataSet * ProcessDataSet( IDataSet*, IPDF* );
		pair< double, double > CalculateMatrixElements( long, long, IDataSet*, vector<double>&, vector<double>&, double&, vector<double>& );

		virtual void SetApplyAlphaCorrection( bool useAlpha );

	private:
		void ConfigurePDFs( IPDF* InputPDF );

		//	Uncopyable!
		SWeightPrecalculator ( const SWeightPrecalculator& );
		SWeightPrecalculator& operator = ( const SWeightPrecalculator& );

		void ApplyAlphaCorrection( IDataSet* );

		FitResult* inputResult;
		IPDF * signalPDF;
		IPDF * backgroundPDF;
		string weightName;
		string fractionName;
		unsigned int config;
		bool useAlpha;
};

#endif

