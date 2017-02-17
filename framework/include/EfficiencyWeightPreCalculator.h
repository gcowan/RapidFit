/**
  @class EfficiencyWeightPreCalculator

  A Precalculator for producing the sWeights for a data set for backgorund subtraction

  @author Benjamin Wynne bwynne@cern.ch
  @date 2009-12-14
  */

#pragma once
#ifndef EFFICIENCY_WEIGHT_PRECALCULATOR_H
#define EFFICIENCY_WEIGHT_PRECALCULATOR_H

//	RapidFit Headers
#include "IPrecalculator.h"
#include "IPDF.h"
#include "FitResult.h"
//	System Headers
#include <vector>
#include <string>

class EfficiencyWeightPreCalculator : public IPrecalculator
{
	public:
		/*!
		 * @brief Constructor
		 *
		 * @Param inputResult
		 *
		 * @Param WeightName
		 *
		 * @Param config
		 */
		EfficiencyWeightPreCalculator( FitResult* inputResult, string WeightName = "sWeight", unsigned int config=1 );

		/*!
		 * @brief Destructor Function
		 */
		~EfficiencyWeightPreCalculator();

		/*!
		 * @brief
		 */
		virtual IDataSet * ProcessDataSet( IDataSet*, IPDF* );
        virtual void SetApplyAlphaCorrection( bool );

	private:
		/*!
		 * @brief
		 */
		void ConfigurePDFs( IPDF* InputPDF );

		//	Uncopyable!
		EfficiencyWeightPreCalculator ( const EfficiencyWeightPreCalculator& );
		EfficiencyWeightPreCalculator& operator = ( const EfficiencyWeightPreCalculator& );

		/*!
		 * Internal Objects used in the calculation of the sWeights
		 */

		FitResult* inputResult;
		IPDF * signalPDF;
		string weightName;
		unsigned int config;
};

#endif

