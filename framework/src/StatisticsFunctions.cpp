/**
  @class StatisticsFunctions

  A collection of static methods for statistics calculations

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

//	RapidFit Headers
#include "StatisticsFunctions.h"
#include "StringProcessing.h"
//	System Headers
#include <math.h>
#include <iostream>
#include <float.h>

//#define DOUBLE_TOLERANCE DBL_MIN
#define DOUBLE_TOLERANCE 1E-6

//Return the mean of a vector of doubles
double StatisticsFunctions::Mean( const vector<double> Numbers )
{
	double sum = 0.0;
	for (unsigned int index = 0; index < Numbers.size(); ++index )
	{
		sum += Numbers[index];
	}

	return sum / double(Numbers.size());
}

//Return the variance of a vector of doubles
double StatisticsFunctions::Variance( const vector<double> Numbers )
{
	double sum = 0.0;
	double mean = Mean(Numbers);
	for (unsigned int index = 0; index < Numbers.size(); ++index )
	{
		sum += pow( ( Numbers[index] - mean ), 2 );
	}

	return sum / ( double(Numbers.size()) - 1 );
}

//Return the ideal number of bins for a histogram of a vector of doubles
//Uses D. Scott's method, published 1979
int StatisticsFunctions::OptimumBinNumber( const vector<double> Numbers )
{
	double width = 3.49 * sqrt( Variance(Numbers) ) * pow( double(Numbers.size()), -(1.0/3.0) );
	double range = Maximum(Numbers) - Minimum(Numbers);
	return (int)ceil( range / width );
}

//Returns the maximum and minimum of a vector of doubles 
double StatisticsFunctions::Maximum( const vector<double> Numbers )
{
	if ( Numbers.size() > 0 )
	{
		double maximum = Numbers[0];
		for (unsigned int index = 0; index < Numbers.size(); ++index )
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

double StatisticsFunctions::Minimum( const vector<double> Numbers )
{       
	if ( Numbers.size() > 0 )
	{
		double minimum = Numbers[0];
		for (unsigned int index = 0; index < Numbers.size(); ++index )
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
vector< vector<double> > StatisticsFunctions::DiscreteCombinations( const vector<string> * AllNames, const PhaseSpaceBoundary * InputBoundary,
	vector<string> & DiscreteNames, vector<string> & ContinuousNames, vector< vector<double> > & discreteValues )
{
	//Construct a vector<vector> containing all discrete values. List the names of discrete and continuous observables.
	vector<string>::const_iterator nameIterator;
	for( nameIterator = AllNames->begin(); nameIterator != AllNames->end(); ++nameIterator )
	{
		IConstraint * observableConstraint = InputBoundary->GetConstraint( *nameIterator );
		if( observableConstraint != NULL )
		{
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
		else return vector<vector<double> >();
	}

	//Make a vector of combinations of discrete values
	vector< vector<double> > discreteCombinations;

	//Check if there actually are any discrete observables
	if ( DiscreteNames.size() > 0 )
	{
		bool finished = false;
		vector<int> indices;
		for (unsigned int observableIndex = 0; observableIndex < discreteValues.size(); ++observableIndex )
		{
			//Just initialise the counter
			indices.push_back(0);
		}
		while( !finished )
		{
			//Store one combination, based on the indices
			vector<double> oneCombination;
			for (unsigned int observableIndex = 0; observableIndex < discreteValues.size(); ++observableIndex )
			{
				oneCombination.push_back( discreteValues[unsigned(observableIndex)][ unsigned(indices[unsigned(observableIndex)]) ] );
			}
			discreteCombinations.push_back(oneCombination);

			//Increment the indices
			for ( int observableIndex = int(discreteValues.size()) - 1; observableIndex >= 0; --observableIndex )
			{
				++indices[unsigned(observableIndex)];

				//Check if the index has reached its maximum
				if ( indices[unsigned(observableIndex)] == int(discreteValues[unsigned(observableIndex)].size()) )
				{
					if ( observableIndex == 0 )
					{
						//If the most significant index has reached its maximum, all combinations have been found
						finished = true;
					}
					else
					{
						//Zero this index, and examine the next most significant
						indices[unsigned(observableIndex)] = 0;
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
void StatisticsFunctions::DoDontIntegrateLists( IPDF * InputPDF, const PhaseSpaceBoundary * InputBoundary, const vector<string> * DontIntegrateThese,
		vector<string> & DoIntegrateList, vector<string> & DontIntegrateList )
{
	//Make lists of observables to integrate and not to integrate
	vector<string> observableNames = InputPDF->GetPrototypeDataPoint();
	for (unsigned int observableIndex = 0; observableIndex < observableNames.size(); ++observableIndex )
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
vector<DataPoint*> StatisticsFunctions::DataAverage( const IDataSet * InputData, const vector< vector<double> > DiscreteCombinations,
	const vector< vector<double> > DiscreteValues, const vector<string> DiscreteNames, const vector<string> ContinuousNames,
	vector<string> & DataPointDescriptions, vector<double> & DataPointWeights )
{
	//Initialise the data averaging
	vector<double> continuousSums;
	vector<long> combinationCounts;
	for (unsigned int continuousIndex = 0; continuousIndex < ContinuousNames.size(); ++continuousIndex )
	{
		continuousSums.push_back(0.0);
	}
	for (unsigned int combinationIndex = 0; combinationIndex < DiscreteCombinations.size(); ++combinationIndex )
	{
		combinationCounts.push_back(0);
	}

	//Examine the data set. Find the average value for each continuous observable, and the weight for each discrete combination
	for ( int dataIndex = 0; dataIndex < InputData->GetDataNumber(); ++dataIndex )
	{
		DataPoint * readDataPoint = InputData->GetDataPoint(dataIndex);

		//Sum the continuous values, in preparation for taking the average
		for (unsigned int continuousIndex = 0; continuousIndex < ContinuousNames.size(); ++continuousIndex )
		{
			continuousSums[continuousIndex] += readDataPoint->GetObservable( ContinuousNames[continuousIndex] )->GetValue();
		}

		//Calculate the index for the discrete combination, and increment the corresponding count
		int combinationIndex = 0;
		int incrementValue = 1;
		for ( int discreteIndex = int(DiscreteNames.size()) - 1; discreteIndex >= 0; --discreteIndex )
		{
			double currentValue = readDataPoint->GetObservable( DiscreteNames[unsigned(discreteIndex)] )->GetValue();

			for (unsigned int valueIndex = 0; valueIndex < DiscreteValues[unsigned(discreteIndex)].size(); ++valueIndex )
			{
				if ( fabs( DiscreteValues[unsigned(discreteIndex)][unsigned(valueIndex)] - currentValue ) < DOUBLE_TOLERANCE )
				{
					combinationIndex += ( incrementValue * int(valueIndex) );
					incrementValue *= int(DiscreteValues[unsigned(discreteIndex)].size());
					break;
				}
			}
		}
		++combinationCounts[unsigned(combinationIndex)];
	}

	//Calculate averages and weights
	vector<double> combinationWeights;
	double dataNumber = (double)InputData->GetDataNumber();
	for (unsigned int continuousIndex = 0; continuousIndex < ContinuousNames.size(); ++continuousIndex )
	{
		continuousSums[continuousIndex] /= dataNumber;
	}

	for (unsigned int combinationIndex = 0; combinationIndex < DiscreteCombinations.size(); ++combinationIndex )
	{
		combinationWeights.push_back( (double)combinationCounts[combinationIndex] / dataNumber );
	}

	//Create the data points to return
	vector<DataPoint*> newDataPoints;
	vector<string> allDescriptions;
	DataPoint* templateDataPoint =  new DataPoint( *InputData->GetDataPoint(0) );
	for (unsigned int continuousIndex = 0; continuousIndex < ContinuousNames.size(); ++continuousIndex )
	{
		Observable * newValue = templateDataPoint->GetObservable( ContinuousNames[continuousIndex] );
		Observable* newValue2 = new Observable( newValue->GetName(), continuousSums[continuousIndex], newValue->GetUnit() );
		templateDataPoint->SetObservable( ContinuousNames[continuousIndex], newValue2 );
		delete newValue2;
	}
	for (unsigned int combinationIndex = 0; combinationIndex < DiscreteCombinations.size(); ++combinationIndex )
	{
		string description = "(";

		//Output the discrete values for this combination
		for (unsigned int discreteIndex = 0; discreteIndex < DiscreteNames.size(); ++discreteIndex )
		{
			//Set the data point
			Observable * newValue = templateDataPoint->GetObservable( DiscreteNames[discreteIndex] );
			Observable* newValue2 = new Observable( newValue->GetName(), DiscreteCombinations[combinationIndex][discreteIndex], newValue->GetUnit() );
			templateDataPoint->SetObservable( DiscreteNames[discreteIndex], newValue2 );
			delete newValue2;

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

vector<vector<double> > StatisticsFunctions::Combinatorix( const vector<vector<double> > input_values )
{
	unsigned int output_size=1;
	for( unsigned int i=0; i< input_values.size(); ++i )
	{
		output_size *= (unsigned)input_values[i].size();
	}

	//cout << output_size << "\t" << input_values.size() << endl;
	vector<vector<double> > output( output_size, vector<double>( input_values.size(), 0. ) );

	/*cout << "input:" << endl;
	for( unsigned int i=0; i<input_values.size(); ++i )
	{
		for( unsigned int j=0; j<input_values[i].size(); ++j )
		{
			cout << "\t" << input_values[i][j];// << endl;
		}
		cout << endl;
	}*/

	for( unsigned int i=0; i< input_values.size(); ++i )
	{
		//unsigned int modulo = NumberOfCombinations( input_values, i );
		unsigned int repeats = (unsigned)output_size/(unsigned)input_values[i].size();
		//cout << "i: " << i << "\tmodulo: " << modulo << "\trepeats: " << repeats << endl << endl;

		for( unsigned int repeat_num=0; repeat_num < repeats; ++repeat_num )
		{
			//cout << input_values[i].size() << endl;
			for( unsigned int k=0; k< input_values[i].size(); ++k )
			{
				//cout << repeat_num*(input_values[i].size()-1)+k << "\t" << i << "\t\t" << repeat_num*(input_values[i].size()-1) << "\t" << k << endl;
				//output[repeat_num*(input_values[i].size()-1)+k][i] = 0.;
				output[repeat_num*input_values[i].size()+k][i] = input_values[i][k];
				//cout << (repeat_num*input_values[i].size()+k) << "\t" << i << "\t" << output[repeat_num*input_values[i].size()+k][i] << endl;
			}
		}
	}

	return output;
}

unsigned int StatisticsFunctions::NumberOfCombinations( const vector<vector<double> > input, const int index )
{
	unsigned int output=1;
	for( unsigned int i=0; i< input.size(); ++i )
	{
		//cout << index << "\t" << i << "\tsize: " << input[i].size() << endl;
		if( (int)i != index ) output *= (unsigned)input[i].size();
	}
	return output;
}

