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
vector< vector<double> > StatisticsFunctions::DiscreteCombinations( vector<string> * AllNames, PhaseSpaceBoundary * InputBoundary, vector<string> & DiscreteNames,
		vector<string> & ContinuousNames, vector< vector<double> > & discreteValues )
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

//Perform data averaging
vector<DataPoint> StatisticsFunctions::DataAverage( IDataSet * InputData, vector< vector<double> > DiscreteCombinations, vector< vector<double> > DiscreteValues, vector<string> DiscreteNames, vector<string> ContinuousNames,
	       vector<string> & DataPointDescriptions, vector<double> & DataPointWeights )
{
	//Initialise the data averaging
	vector<double> continuousSums;
	vector<long> combinationCounts;
	for ( int continuousIndex = 0; continuousIndex < ContinuousNames.size(); continuousIndex++ )
	{
		continuousSums.push_back(0.0);
	}
	for ( int combinationIndex = 0; combinationIndex < DiscreteCombinations.size(); combinationIndex++ )
	{
		combinationCounts.push_back(0);
	}

	//Examine the data set. Find the average value for each continuous observable, and the weight for each discrete combination
	for ( int dataIndex = 0; dataIndex < InputData->GetDataNumber(); dataIndex++ )
	{
		DataPoint * readDataPoint = InputData->GetDataPoint(dataIndex);

		//Sum the continuous values, in preparation for taking the average
		for ( int continuousIndex = 0; continuousIndex < ContinuousNames.size(); continuousIndex++ )
		{
			continuousSums[continuousIndex] += readDataPoint->GetObservable( ContinuousNames[continuousIndex] )->GetValue();
		}

		//Calculate the index for the discrete combination, and increment the corresponding count
		int combinationIndex = 0;
		int incrementValue = 1;
		for ( int discreteIndex = DiscreteNames.size() - 1; discreteIndex >= 0; discreteIndex-- )
		{
			double currentValue = readDataPoint->GetObservable( DiscreteNames[discreteIndex] )->GetValue();

			for ( int valueIndex = 0; valueIndex < DiscreteValues[discreteIndex].size(); valueIndex++ )
			{
				if ( DiscreteValues[discreteIndex][valueIndex] == currentValue )
				{
					combinationIndex += ( incrementValue * valueIndex );
					incrementValue *= DiscreteValues[discreteIndex].size();
					break;
				}
			}
		}
		combinationCounts[combinationIndex]++;
	}

	//Calculate averages and weights
	vector<double> combinationWeights;
	double dataNumber = (double)InputData->GetDataNumber();
	for ( int continuousIndex = 0; continuousIndex < ContinuousNames.size(); continuousIndex++ )
	{
		continuousSums[continuousIndex] /= dataNumber;
	}

	for ( int combinationIndex = 0; combinationIndex < DiscreteCombinations.size(); combinationIndex++ )
	{
		combinationWeights.push_back( (double)combinationCounts[combinationIndex] / dataNumber );
	}

	//Create the data points to return
	vector<DataPoint> newDataPoints;
	vector<string> allDescriptions;
	DataPoint templateDataPoint = *( InputData->GetDataPoint(0) );
	for ( int continuousIndex = 0; continuousIndex < ContinuousNames.size(); continuousIndex++ )
	{
		Observable * newValue = templateDataPoint.GetObservable( ContinuousNames[continuousIndex] );
		newValue->SetValue( continuousSums[continuousIndex] );
		templateDataPoint.SetObservable( ContinuousNames[continuousIndex], newValue );
	}
	for ( int combinationIndex = 0; combinationIndex < DiscreteCombinations.size(); combinationIndex++ )
	{
		string description = "(";

		//Output the discrete values for this combination
		for ( int discreteIndex = 0; discreteIndex < DiscreteNames.size(); discreteIndex++ )
		{
			//Set the data point
			Observable * newValue = templateDataPoint.GetObservable( DiscreteNames[discreteIndex] );
			newValue->SetValue( DiscreteCombinations[combinationIndex][discreteIndex] );
			templateDataPoint.SetObservable( DiscreteNames[discreteIndex], newValue );

			//Make the description
			char value[100];
			sprintf( value, "%f", DiscreteCombinations[combinationIndex][discreteIndex] );
			string addToDescription;
			if ( discreteIndex == DiscreteNames.size() - 1 )
			{
				addToDescription = DiscreteNames[discreteIndex] + "=" + value;
			}
			else
			{
				addToDescription = DiscreteNames[discreteIndex] + "=" + value + "_";
			}
			description.append(addToDescription);
		}

		description.append(")");
		allDescriptions.push_back(description);
		newDataPoints.push_back(templateDataPoint);
	}

	//Output the results
	DataPointDescriptions = allDescriptions;
	DataPointWeights = combinationWeights;
	return newDataPoints;
}
