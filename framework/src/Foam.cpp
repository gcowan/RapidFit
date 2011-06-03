/**
  @class Foam

  Class for generating toy data from a PDF.
  Just a wrapper for the Root TFoam class.

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-11-10
 */

//	ROOT Headers
#include "TFile.h"
#include "TFoam.h"
//	RapidFit Headers
#include "Foam.h"
#include "StatisticsFunctions.h"
//	System Headers
#include <iostream>
#include <math.h>
#include <limits.h>

#define DOUBLE_TOLERANCE 1E-6

//Default constructor
Foam::Foam() : Open_Files(), InputPDF(NULL), generationFunction(), generationBoundary(), newDataSet(), rootRandom(), foamGenerators(), storedIntegrator(), storedDatapoint(), dataNumber(0), discreteCombinations(), allNames(), discreteNames(), continuousNames(), discreteValues(), minima(), ranges()
{
}

//Constructor with correct argument
Foam::Foam( PhaseSpaceBoundary * NewBoundary, IPDF * NewPDF ) :  Open_Files(), InputPDF(NewPDF), generationFunction(), generationBoundary(NewBoundary), newDataSet(), rootRandom(), foamGenerators(), storedIntegrator(), storedDatapoint(), dataNumber(0), discreteCombinations(), allNames(), discreteNames(), continuousNames(), discreteValues(), minima(), ranges()
{
	rootRandom = NewPDF->GetRandomFunction();

	//Retrieve all combinations of discrete variables
	allNames = generationBoundary->GetAllNames();
	discreteCombinations = StatisticsFunctions::DiscreteCombinations( &allNames, generationBoundary, discreteNames, continuousNames, discreteValues );

        newDataSet = new MemoryDataSet(generationBoundary);
        cout << "Initializing Generator(s)" << endl;

        Foam::Init();
}

void Foam::Init()
{
//Retrieve the maxima and minima, to use in the coordinate transforms
	for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
	{
		double maximum = generationBoundary->GetConstraint( continuousNames[continuousIndex] )->GetMaximum();
		double minimum = generationBoundary->GetConstraint( continuousNames[continuousIndex] )->GetMinimum();
		minima.push_back(minimum);
		ranges.push_back( maximum - minimum );
	}

	//Make a Foam generator for each combination
	for (unsigned int combinationIndex = 0; combinationIndex < discreteCombinations.size(); ++combinationIndex )
	{
		//Make the data point for the combination
		DataPoint * init_temporaryDataPoint = new DataPoint(allNames);
		for (unsigned int discreteIndex = 0; discreteIndex < discreteNames.size(); ++discreteIndex )
		{
			string unit = generationBoundary->GetConstraint( discreteNames[discreteIndex] )->GetUnit();
			init_temporaryDataPoint->SetObservable( discreteNames[discreteIndex], discreteCombinations[combinationIndex][discreteIndex], 0.0, unit );
		}
		for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
		{
			string unit = generationBoundary->GetConstraint( continuousNames[continuousIndex] )->GetUnit();
			init_temporaryDataPoint->SetObservable( continuousNames[continuousIndex], 0.0, 0.0, unit );
		}

		//Make the function wrapper
		IntegratorFunction * combinationFunction = new IntegratorFunction( InputPDF, init_temporaryDataPoint, continuousNames, discreteNames, minima, ranges );

		//cout << InputPDF->GET_ID() << endl;
		//if( InputPDF->GetCacheStatus() )	cout << "HAS CACHED FOAM" << endl;
		//else	cout << "HAS NOT CACHED FOAM" << endl;

		//	Generate unique name for this foam instance
		TString Name("Foam-");
		Name.Append(InputPDF->GET_ID());
		Name.Append("-");
		Name+=combinationIndex;
		//	Create a root file to store the TFoam universe
		TString RootName(Name);
		RootName.Append(".root");

		TFile* MC_Cache = NULL;
		TFoam * foamGenerator = NULL;

		if( !InputPDF->GetMCCacheStatus() )
		{
			cout << "FOAM NOT Cached, Generating Foam:\t" << Name << endl;
			MC_Cache = new TFile( RootName, "RECREATE" );
			//Initialise Foam
			foamGenerator = new TFoam( Name );
			foamGenerator->SetkDim( Int_t(continuousNames.size()) );
			foamGenerator->SetPseRan( rootRandom );
			foamGenerator->SetRho( combinationFunction );	//	Can afford to Boot Foam's ability if we're using just one cached instance :D
			foamGenerator->SetnCells( 1000 );	//	1000	Total number of bins to construct
			foamGenerator->SetnSampl( 200 );	//	200	Samples to take when constructing bins
			foamGenerator->SetnBin( 8 );		//	8	Bins along each axis
			foamGenerator->SetOptRej( 1 );		//	1/0	Don't/Use Weighted Distribution
			foamGenerator->SetOptDrive( 2 );	//	1/2	Best Varience/Weights
			foamGenerator->SetEvPerBin( 25 );	//	25	Weights per bin... This doesn't Saturate as object is written before generating events
			foamGenerator->SetChat( 0 );		//	0	verbosity
			foamGenerator->SetMaxWtRej( 1.1 );	//	1.1	Unknown what effect this has, something to do with weights
			foamGenerator->Initialize();
//			foamGenerator->ResetPseRan( rootRandom );
//			foamGenerator->ResetRho( combinationFunction );	//	Can afford to Boot Foam's ability if we're using just one cached instance :D
//			foamGenerator->Initialize();
			//	As we haven't cached yet, write to file
			foamGenerator->Write(Name);
			cout << "Storing TFOAM TObject in:\t\t" << RootName << endl;
			InputPDF->AddCacheObject( RootName );
			MC_Cache->Write();
		} else {
			cout << "FOAM Cached, Re-Aquiring TFOAM TObject:\t" << Name << endl;
			MC_Cache = new TFile( RootName, "UPDATE" );
			//gDirectory->ls();
			//Cached_Files.back()->Map();
			//Cached_Files.back()->ShowStreamerInfo();
			//Cached_Files.back()->GetListOfKeys()->Print();
			foamGenerator = (TFoam*) MC_Cache->Get( Name );
			//	This code has been production testeted enough such that I think this is pointlessly scaring users
			//cout << "Checking Consistancy, will crash if this fails... this is unavoidable" << endl;
			foamGenerator->SetPseRan( rootRandom );
			foamGenerator->SetRho( combinationFunction );
			foamGenerator->SetChat( 0 );
			//	Make one event here to check everything was processed correctly
			foamGenerator->MakeEvent();
			cout << "Check OK" << endl;
		}

		Open_Files.push_back( MC_Cache );
		//Debug
		//char filename[100];
		//sprintf( filename, "DebugFoamPlot.%d.c", combinationIndex );
		//foamGenerator->RootPlot2dim(filename);
		//foamGenerator->PrintCells();

		//Store the foam generator
		foamGenerators.push_back( foamGenerator );
		storedIntegrator.push_back( combinationFunction );
		storedDatapoint.push_back( init_temporaryDataPoint );
	}
	InputPDF->SetMCCacheStatus( true );

}

//Destructor
Foam::~Foam()
{
	//Foam::RemoveGenerator();
	//delete rootRandom;
}

void Foam::RemoveGenerator()
{
	while( !foamGenerators.empty() )
	{
		Double_t one, two;
		foamGenerators.back()->Finalize( one, two );
		delete foamGenerators.back();
		foamGenerators.pop_back();
	}
	while( !storedIntegrator.empty() )
	{
		delete storedIntegrator.back();
		storedIntegrator.pop_back();
	}
	while( !storedDatapoint.empty() )
	{
		delete storedDatapoint.back();
		storedDatapoint.pop_back();
	}

	while( !Open_Files.empty() )
	{
		Open_Files.back()->Close();
		delete Open_Files.back();
		Open_Files.pop_back();
	}
}

//Use accept/reject method to create data
int Foam::GenerateData( int DataAmount )
{
	newDataSet->ReserveDataSpace( DataAmount );
	for (int dataIndex = 0; dataIndex < DataAmount; ++dataIndex )
	{
		//Generate the discrete observables, and select the correct Foam generator to use
		int combinationIndex = 0;
		int incrementValue = 1;
		DataPoint * temporaryDataPoint = new DataPoint(allNames);
		for ( int discreteIndex = int(discreteNames.size() - 1); discreteIndex >= 0; --discreteIndex )
		{
			//Create the discrete observable
			Observable * temporaryObservable = generationBoundary->GetConstraint( discreteNames[unsigned(discreteIndex)] )->CreateObservable(rootRandom);
			double currentValue = temporaryObservable->GetValue();
			temporaryDataPoint->SetObservable( discreteNames[unsigned(discreteIndex)], temporaryObservable );

			//Calculate the index
			for (unsigned int valueIndex = 0; valueIndex < discreteValues[unsigned(discreteIndex)].size(); ++valueIndex )
			{
				if ( fabs( discreteValues[unsigned(discreteIndex)][valueIndex] - currentValue ) < DOUBLE_TOLERANCE )
				{
					combinationIndex += ( incrementValue * int(valueIndex) );
					incrementValue *= int(discreteValues[unsigned(discreteIndex)].size());
					break;
				}
			}
		}

		//Use the index calculated to select a Foam generator and generate an event with it
		Double_t* generatedEvent = new Double_t[ continuousNames.size() ];
		foamGenerators[unsigned(combinationIndex)]->MakeEvent();
		foamGenerators[unsigned(combinationIndex)]->GetMCvect(generatedEvent);

		//Store the continuous observables
		for (unsigned int continuousIndex = 0; continuousIndex < continuousNames.size(); ++continuousIndex )
		{
			string unit = generationBoundary->GetConstraint( continuousNames[continuousIndex] )->GetUnit();
			double newValue = minima[continuousIndex] + ( ranges[continuousIndex] * generatedEvent[continuousIndex] );
			temporaryDataPoint->SetObservable( continuousNames[continuousIndex], newValue, 0.0, unit );
			//cout << continuousNames[continuousIndex] << "\t" << newValue << "\t" << 0.0 << "\t" << unit << endl;
		}

		delete generatedEvent;
		//	Store the event
		newDataSet->AddDataPoint(temporaryDataPoint);
		//	Passed by reference but copied into memory
		delete temporaryDataPoint;
		++dataNumber;
	}

	//cout << "Destroying Generator(s)" << endl;
	Foam::RemoveGenerator();
	return DataAmount;
}

//Return data set
IDataSet * Foam::GetDataSet()
{
	return newDataSet;
}
