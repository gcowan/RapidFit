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
		/*!
		 * @brief Constructor
		 *
		 * @Param inputResult
		 *
		 * @Param WeightName
		 *
		 * @Param config
		 */
		SWeightPrecalculator( FitResult* inputResult, string WeightName = "sWeight", unsigned int config=1 );

		/*!
		 * @brief Destructor Function
		 */
		~SWeightPrecalculator();

		/*!
		 * @brief 
		 */
		virtual IDataSet * ProcessDataSet( IDataSet*, IPDF* );

		/*!
		 * @brief This calculates useful values in order to calculate the sWeights
		 *
		 * @param numSignal  This is the number of signal events in the dataset
		 *
		 * @param numBack  This is the number of background events in the dataset
		 *
		 * @param DataSet  This is the Dataset to calclate the Matrix Elements for
		 *
		 * @param signalValues
		 *
		 * @param backValues
		 *
		 * @param denom2
		 *
		 * @param numer2
		 *
		 * @return pair of SignalSignal and BackgroundBackground values
		 */
		pair< double, double > CalculateMatrixElements( double numSignal, double numBack, IDataSet* DataSet,
							vector<double>& signalValues, vector<double>& backValues, double& denom2, vector<double>& numer2 );

		/*!
		 * @brief Function to decide whether the alpha correction should be applied to these weights before they're stored to disk
		 * (This modifies the weights with a linear factor!)
		 *
		 * @param useAlpha This is the alpha correction to be applied to the calculated sWeights
		 *
		 * @return void
		 */
		virtual void SetApplyAlphaCorrection( bool useAlpha );

	private:
		/*!
		 * @brief 
		 */
		void ConfigurePDFs( IPDF* InputPDF );

		//	Uncopyable!
		SWeightPrecalculator ( const SWeightPrecalculator& );
		SWeightPrecalculator& operator = ( const SWeightPrecalculator& );

		/*!
		 * @brief
		 */
		void ApplyAlphaCorrection( IDataSet* );

		/*!
		 * Internal Objects used in the calculation of the sWeights
		 */

		FitResult* inputResult;
		IPDF * signalPDF;
		IPDF * backgroundPDF;
		string weightName;
		string fractionName;
		unsigned int config;
		bool useAlpha;
};

#endif

