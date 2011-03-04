#include "FoamIntegrator.h"
#include "StatisticsFunctions.h"
#include <math.h>

#define DOUBLE_TOLERANCE 1E-6

//Default constructor
FoamIntegrator::FoamIntegrator()
{
}

//Constructor with correct arguments
FoamIntegrator::FoamIntegrator( IPDF * InputPDF, IDataSet * InputData )
{
	//Calculate all possible combinations of discrete observables
	vector<string> continuousNames;
	vector<string> allNames = InputPDF->GetPrototypeDataPoint();
	vector< vector<double> > discreteCombinations = StatisticsFunctions::DiscreteCombinations( &allNames, InputData->GetBoundary(), discreteNames, continuousNames, discreteValues );

	//Create data points for each combination (using data averaged values for continuous, non-integrable functions)
	vector<string> dataPointDescriptions;
	vector<double> dataPointWeights;
	vector<DataPoint> combinationPoints =  StatisticsFunctions::DataAverage( InputData, discreteCombinations, discreteValues, discreteNames, continuousNames, dataPointDescriptions, dataPointWeights );

	//Make a foam for each discrete combination
	for (unsigned int combinationIndex = 0; combinationIndex < combinationPoints.size(); combinationIndex++ )
	{
		MakeFoam combinationFoam( InputPDF, InputData->GetBoundary(), &( combinationPoints[combinationIndex] ) );
		allIntegrators.push_back(combinationFoam);
	}	
}

//Destructor
FoamIntegrator::~FoamIntegrator()
{
}

//Select and run the correct integrator
double FoamIntegrator::Integral( DataPoint * InputPoint, PhaseSpaceBoundary * InputBoundary )
{
	PhaseSpaceBoundary* null_p = InputBoundary;
	null_p = NULL;
	//The integral won't work if the boundary has changed, but you might want a check that it's the same

	//Use the data point to find the index of the correct foam
	int combinationIndex = 0;
	int incrementValue = 1;
	for ( int discreteIndex = discreteNames.size() - 1; discreteIndex >= 0; discreteIndex-- )
	{
		//Retrieve the observable value
		Observable * temporaryObservable = InputPoint->GetObservable( discreteNames[discreteIndex] );
		double currentValue = temporaryObservable->GetValue();

		//Calculate the index
		for (unsigned int valueIndex = 0; valueIndex < discreteValues[discreteIndex].size(); valueIndex++ )
		{
			if ( fabs(discreteValues[discreteIndex][valueIndex] - currentValue ) < DOUBLE_TOLERANCE )
			{
				combinationIndex += ( incrementValue * valueIndex );
				incrementValue *= discreteValues[discreteIndex].size();
				break;
			}
		}
	}

	//Use the foam to integrate
	return allIntegrators[combinationIndex].Integral();
}
