/**
  @class Foam

  Class for generating toy data from a PDF.
  Just a wrapper for the Root TFoam class.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-11-10
 */

#include "Foam.h"
#include "InvalidObject.h"
#include "StatisticsFunctions.h"
#include <iostream>

//Default constructor
Foam::Foam()
{
}

//Constructor with correct argument
Foam::Foam( PhaseSpaceBoundary * NewBoundary, IPDF * NewPDF ) : generationBoundary(NewBoundary), dataNumber(0)
{
	newDataSet = new MemoryDataSet(generationBoundary);
	rootRandom = new TRandom3(0);

	//Retrieve all combinations of discrete variables
	allNames = generationBoundary->GetAllNames();
	vector< vector<double> > discreteCombinations = StatisticsFunctions::DiscreteCombinations( &allNames, generationBoundary, discreteNames, continuousNames, discreteValues );
	
	//Retrieve the maxima and minima, to use in the coordinate transforms
	for ( int continuousIndex = 0; continuousIndex < continuousNames.size(); continuousIndex++ )
	{
		double maximum = generationBoundary->GetConstraint( continuousNames[continuousIndex] )->GetMaximum();
		double minimum = generationBoundary->GetConstraint( continuousNames[continuousIndex] )->GetMinimum();

		minima.push_back(minimum);
		ranges.push_back( maximum - minimum );
	}

	//Make a Foam generator for each combination
	for ( int combinationIndex = 0; combinationIndex < discreteCombinations.size(); combinationIndex++ )
	{
		//Make the data point for the combination
		DataPoint * temporaryDataPoint = new DataPoint(allNames);
		for ( int discreteIndex = 0; discreteIndex < discreteNames.size(); discreteIndex++ )
		{
			string unit = generationBoundary->GetConstraint( discreteNames[discreteIndex] )->GetUnit();
			temporaryDataPoint->SetObservable( discreteNames[discreteIndex], discreteCombinations[combinationIndex][discreteIndex], 0.0, unit );
			cout << discreteNames[discreteIndex] << " value: " << discreteCombinations[combinationIndex][discreteIndex] << " unit: " << unit << endl;
		}
		for ( int continuousIndex = 0; continuousIndex < continuousNames.size(); continuousIndex++ )
		{
			string unit = generationBoundary->GetConstraint( continuousNames[continuousIndex] )->GetUnit();
			temporaryDataPoint->SetObservable( continuousNames[continuousIndex], 0.0, 0.0, unit );
			cout << continuousNames[continuousIndex] << " value: 0.0 unit: " << unit << endl;
		}

		//Make the function wrapper
		IntegratorFunction * combinationFunction = new IntegratorFunction( NewPDF, temporaryDataPoint, continuousNames, discreteNames, minima, ranges );

		//Initialise Foam
		TFoam * foamGenerator = new TFoam();
		foamGenerator->SetkDim( continuousNames.size() );
		foamGenerator->SetPseRan(rootRandom);
		foamGenerator->SetRho(combinationFunction);
		foamGenerator->SetnCells(1000);
		foamGenerator->SetnSampl(200);
		foamGenerator->SetnBin(8);
		foamGenerator->SetOptRej(1);
		foamGenerator->SetOptDrive(2);
		foamGenerator->SetEvPerBin(25);
		foamGenerator->SetChat(1);
		foamGenerator->SetMaxWtRej(1.1);
		foamGenerator->Initialize();

		//Store the foam generator
		foamGenerators.push_back(foamGenerator);
	}
}

//Destructor
Foam::~Foam()
{
	for ( int generatorIndex = 0; generatorIndex < foamGenerators.size(); generatorIndex++ )
	{
		delete foamGenerators[generatorIndex];
	}
	delete rootRandom;
}

//Use accept/reject method to create data
int Foam::GenerateData( int DataAmount )
{
	for ( int dataIndex = 0; dataIndex < DataAmount; dataIndex++ )
	{
		//Generate the discrete observables, and select the correct Foam generator to use
		int combinationIndex = 0;
		int incrementValue = 1;
		DataPoint * temporaryDataPoint = new DataPoint(allNames);
		for ( int discreteIndex = discreteNames.size() - 1; discreteIndex >= 0; discreteIndex-- )
		{
			//Create the discrete observable
			Observable * temporaryObservable = generationBoundary->GetConstraint( discreteNames[discreteIndex] )->CreateObservable(rootRandom);
			double currentValue = temporaryObservable->GetValue();
			temporaryDataPoint->SetObservable( discreteNames[discreteIndex], temporaryObservable );

			//Calculate the index
			for ( int valueIndex = 0; valueIndex < discreteValues[discreteIndex].size(); valueIndex++ )
			{
				if ( discreteValues[discreteIndex][valueIndex] == currentValue )
				{
					combinationIndex += ( incrementValue * valueIndex );
					incrementValue *= discreteValues[discreteIndex].size();
					break;
				}
			}
		}

		//Use the index calculated to select a Foam generator and generate an event with it
		Double_t generatedEvent[ continuousNames.size() ];
		foamGenerators[combinationIndex]->MakeEvent();
		foamGenerators[combinationIndex]->GetMCvect(generatedEvent);

		//Store the continuous observables
		for ( int continuousIndex = 0; continuousIndex < continuousNames.size(); continuousIndex++ )
		{
			string unit = generationBoundary->GetConstraint( continuousNames[continuousIndex] )->GetUnit();
			double newValue = minima[continuousIndex] + ( ranges[continuousIndex] * generatedEvent[continuousIndex] );
			temporaryDataPoint->SetObservable( continuousNames[continuousIndex], newValue, 0.0, unit );
		}

		//Store the event
		newDataSet->AddDataPoint(temporaryDataPoint);
		dataNumber++;
	}

	return DataAmount;
}

//Return data set
IDataSet * Foam::GetDataSet()
{
	return newDataSet;
}
