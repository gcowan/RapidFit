/**
  @class StatisticsFunctions

  A collection of static methods for statistics calculations

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

#include "StatisticsFunctions.h"
#include "math.h"
#include "StringProcessing.h"
#include <iostream>

//Return the mean of a vector of doubles
double StatisticsFunctions::Mean( vector<double> Numbers )
{
	double sum = 0.0;
	for ( int index = 0; index < Numbers.size(); index++ )
	{
		sum += Numbers[index];
	}

	return sum / Numbers.size();
}

//Return the variance of a vector of doubles
double StatisticsFunctions::Variance( vector<double> Numbers )
{
	double sum = 0.0;
	double mean = Mean(Numbers);
	for ( int index = 0; index < Numbers.size(); index++ )
	{
		sum += pow( ( Numbers[index] - mean ), 2 );
	}

	return sum / ( Numbers.size() - 1 );
}

//Return the ideal number of bins for a histogram of a vector of doubles
//Uses D. Scott's method, published 1979
int StatisticsFunctions::OptimumBinNumber( vector<double> Numbers )
{
	double width = 3.49 * sqrt( Variance(Numbers) ) * pow( Numbers.size(), -(1.0/3.0) );
	double range = Maximum(Numbers) - Minimum(Numbers);
	return (int)ceil( range / width );
}

//Returns the maximum and minimum of a vector of doubles 
double StatisticsFunctions::Maximum( vector<double> Numbers )
{
	if ( Numbers.size() > 0 )
	{
		double maximum = Numbers[0];
		for ( int index = 0; index < Numbers.size(); index++ )
		{
			if ( Numbers[index] > maximum )
			{
				maximum = Numbers[index];
			}
		}

		return maximum;
	}
	else
	{
		return 0.0;
	}
}
double StatisticsFunctions::Minimum( vector<double> Numbers )
{       
	if ( Numbers.size() > 0 )
	{
		double minimum = Numbers[0];
		for ( int index = 0; index < Numbers.size(); index++ )
		{
			if ( Numbers[index] < minimum )
			{
				minimum = Numbers[index];
			}
		}

		return minimum;
	}
	else
	{
		return 0.0;
	}
}

//Returns all possible combinations of discrete observable values from a PhaseSpaceBoundary
vector< vector<double> > StatisticsFunctions::DiscreteCombinations( vector<string> * AllNames, PhaseSpaceBoundary * InputBoundary, vector<string> & DiscreteNames, vector<string> & ContinuousNames, vector< vector<double> > & discreteValues )
{
	//Construct a vector<vector> containing all discrete values. List the names of discrete and continuous observables.
	vector<string>::iterator nameIterator;
	for ( nameIterator = AllNames->begin(); nameIterator != AllNames->end(); nameIterator++ )
	{
		IConstraint * observableConstraint = InputBoundary->GetConstraint( *nameIterator );
		if ( observableConstraint->IsDiscrete() )
		{
			discreteValues.push_back( observableConstraint->GetValues() );
			DiscreteNames.push_back( *nameIterator );
		}
		else
		{
			ContinuousNames.push_back( *nameIterator );
		}
	}

	//Make a vector of combinations of discrete values
	vector< vector<double> > discreteCombinations;

	//Check if there actually are any discrete observables
	if ( DiscreteNames.size() > 0 )
	{
		bool finished = false;
		vector<int> indices;
		for ( int observableIndex = 0; observableIndex < discreteValues.size(); observableIndex++ )
		{
			//Just initialise the counter
			indices.push_back(0);
		}
		while( !finished )
		{
			//Store one combination, based on the indices
			vector<double> oneCombination;
			for ( int observableIndex = 0; observableIndex < discreteValues.size(); observableIndex++ )
			{
				oneCombination.push_back( discreteValues[observableIndex][ indices[observableIndex] ] );
			}
			discreteCombinations.push_back(oneCombination);

			//Increment the indices
			for ( int observableIndex = discreteValues.size() - 1; observableIndex >= 0; observableIndex-- )
			{
				indices[observableIndex]++;

				//Check if the index has reached its maximum
				if ( indices[observableIndex] == discreteValues[observableIndex].size() )
				{
					if ( observableIndex == 0 )
					{
						//If the most significant index has reached its maximum, all combinations have been found
						finished = true;
					}
					else
					{
						//Zero this index, and examine the next most significant
						indices[observableIndex] = 0;
					}
				}
				else
				{
					break;
				}
			}
		}
	}
	else
	{
		//Just make an empty entry
		vector<double> noDiscreteVariables;
		discreteCombinations.push_back(noDiscreteVariables);
	}

	return discreteCombinations;
}

//Return names of observables to integrate or not
void StatisticsFunctions::DoDontIntegrateLists( IPDF * InputPDF, PhaseSpaceBoundary * InputBoundary, vector<string> * DontIntegrateThese, vector<string> & DoIntegrateList, vector<string> & DontIntegrateList )
{
	//Make lists of observables to integrate and not to integrate
	vector<string> observableNames = InputPDF->GetPrototypeDataPoint();
	for ( int observableIndex = 0; observableIndex < observableNames.size(); observableIndex++ )
	{
		bool continuous = !( InputBoundary->GetConstraint( observableNames[observableIndex] )->IsDiscrete() );
		bool integrate = ( StringProcessing::VectorContains( DontIntegrateThese, &( observableNames[observableIndex] ) ) == -1 );

		if ( continuous && integrate )
		{
			DoIntegrateList.push_back( observableNames[observableIndex] );
		}
		else
		{
			DontIntegrateList.push_back( observableNames[observableIndex] );
		}
	}
}
