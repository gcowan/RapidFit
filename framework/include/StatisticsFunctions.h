/**
        @class StatisticsFunctions

        A collection of static methods for statistics calculations

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef STATISTICS_FUNCTIONS_H
#define STATISTICS_FUNCTIONS_H

//	RapidFit Headers
#include "PhaseSpaceBoundary.h"
#include "IPDF.h"
#include "IDataSet.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class StatisticsFunctions
{
	public:
		/*!
		 * @brief Returns the mean of the vector that has been passed
		 *
		 * @param inputList	The vector of elements to calculate the mean for
		 *
		 * @return Returns the mean of the values passed in inputList
		 */
		static double Mean( const vector<double> inputList );

		/*!
		 * @brief Returns the standard variance of the vector of elements that has been passed
		 *
		 * @param inputList	The vector of elements to calculate the variance for
		 *
		 * @return Returns the variance of the values passed in inputList
		 */
		static double Variance( const vector<double> inputList );

		/*!
		 * @brief This calculates the 'Optimal' bin number according to D. Scott's method, published 1979
		 *
		 * @param inputList	This is the vector of elements to calculate the 
		 *
		 * @return Returns the optimal bin number assuming the input values are doubles
		 */
		static int OptimumBinNumber( const vector<double> inputList );

		/*!
		 * @brief This returns the maximum element from the passed vector
		 *
		 * @param inputList	This is the vector to find the maximum element from
		 *
		 * @return Returns the maximum element from inputList
		 */
		static double Maximum( const vector<double> inputList );

		/*!
		 * @brief This returns the minimum element from the passed vector
		 *
		 * @param inputList	This is the vector to find the minimum element from
		 *
		 * @return Returns the minimum element from inputList
		 */
		static double Minimum( const vector<double> inputList );

		/*!
		 * @brief Returns all possible combinations of discrete observable values from a PhaseSpaceBoundary
		 *
		 * @param allNames
		 *
		 * @param inputPhaseSpace
		 *
		 * @param discreteNames
		 *
		 * @param continuousNames
		 *
		 * @param discreteValues
		 *
		 * @return
		 */
		static vector< vector<double> > DiscreteCombinations( const vector<string>* allNames, const PhaseSpaceBoundary* inputPhaseSpace,
				vector<string>& discreteNames, vector<string>& continuousNames, vector< vector<double> >& discreteValues );

		/*!
		 * @brief
		 *
		 * @param inputPDF
		 *
		 * @param inputPhaseSpace
		 *
		 * @param DontIntegrateThese
		 *
		 * @param DoIntegrateList
		 *
		 * @param DontIntegrateList
		 *
		 * @return
		 */
		static void DoDontIntegrateLists( IPDF* inputPDF, const PhaseSpaceBoundary* inputPhaseSpace, const vector<string>* DontIntegrateThese,
				vector<string>& DoIntegrateList, vector<string>& DontIntegrateList );

		/*!
		 * @brief
		 *
		 * @param inputDataSet
		 *
		 * @param DiscreteCombinations
		 *
		 * @param DiscreteValues
		 *
		 * @param DiscreteNames
		 *
		 * @param ContinuousNames
		 *
		 * @param DataPointDescriptions
		 *
		 * @param DataPointWeights
		 *
		 * @return
		 */
		static vector<DataPoint*> DataAverage( const IDataSet* inputDataSet, const vector< vector<double> > DiscreteCombinations, const vector< vector<double> > DiscreteValues,
				const vector<string> DiscreteNames, const vector<string> ContinuousNames, vector<string>& DataPointDescriptions, vector<double>& DataPointWeights );

		/*!
		 * @brief
		 *
		 * @param input_values
		 *
		 * @return
		 */
		static vector<vector<double> > Combinatorix( const vector<vector<double> > input_values );

		/*!
		 * @brief Returns the number of possible combinations of discrete constraints, ignoring a pre-defined constraint
		 *
		 * @param input		These are the possible values associated with different discrete constraints
		 *
		 * @param index		This is the optional discrete constraint to ignore for this calculation
		 *
		 * @return Returns the number of different constraints due to the different combinations
		 */
		static unsigned int NumberOfCombinations( const vector<vector<double> > input, const int index=-1 );

	private:

		StatisticsFunctions();
};

#endif

