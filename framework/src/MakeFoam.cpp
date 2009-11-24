#include "MakeFoam.h"
#include <stack>
#include "StatisticsFunctions.h"
#include <iostream>
#include <cmath>

using namespace std;

const int MAXIMUM_SAMPLES = 100;
const double MAXIMUM_GRADIENT_TOLERANCE = 0.01;
const int MAXIMUM_CELLS = 100;

//Default constructor
MakeFoam::MakeFoam()
{
}

//Constructor with correct argument
MakeFoam::MakeFoam( IPDF * InputPDF, PhaseSpaceBoundary * InputBoundary, DataPoint * InputPoint )
{
	//Make the container to hold the possible cells
	stack<PhaseSpaceBoundary> possibleCells;
	PhaseSpaceBoundary firstCell = *InputBoundary;
	possibleCells.push(firstCell);

	//Make a list of observables to integrate over
	vector<string> doIntegrate, dontIntegrate;
	vector<string> pdfDontIntegrate = InputPDF->GetDoNotIntegrateList();
	StatisticsFunctions::DoDontIntegrateLists( InputPDF, InputBoundary, &(pdfDontIntegrate), doIntegrate, dontIntegrate );

	//Continue until all possible cells have been examined
	while ( !possibleCells.empty() && finishedCells.size() < MAXIMUM_CELLS )
	{
		//cout << endl << endl << "Start loop" << endl;

		//Remove the next possible cell
		PhaseSpaceBoundary currentCell = possibleCells.top();
		possibleCells.pop();

		//MC sample the cell, store the DataPoints producing the highest and lowest values
		double maximumValue, minimumValue;
		DataPoint maximumPoint, minimumPoint;
		for ( int sampleIndex = 0; sampleIndex < MAXIMUM_SAMPLES; sampleIndex++ )
		{
			//Create a data point within the current cell
			DataPoint samplePoint( InputPoint->GetAllNames() );
			for ( int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
			{
				//Generate random values to explore integrable observables
				IConstraint * temporaryConstraint = currentCell.GetConstraint( doIntegrate[observableIndex] );
				samplePoint.SetObservable( doIntegrate[observableIndex], temporaryConstraint->CreateObservable() );
			}
			for ( int observableIndex = 0; observableIndex < dontIntegrate.size(); observableIndex++ )
			{
				//Use given values for unintegrable observables
				Observable * temporaryObservable = InputPoint->GetObservable( dontIntegrate[observableIndex] );
				samplePoint.SetObservable( dontIntegrate[observableIndex], temporaryObservable );
			}

			//Evaluate the function at this point
			double sampleValue = InputPDF->Evaluate( &samplePoint );

			//Update maximum and minimum
			if ( sampleIndex == 0 )
			{
				//Populate maximum and minimum
				maximumValue = sampleValue;
				minimumValue = sampleValue;
				maximumPoint = samplePoint;
				minimumPoint = samplePoint;
			}
			else
			{
				//Find maximum
				if ( sampleValue > maximumValue )
				{
					maximumValue = sampleValue;
					maximumPoint = samplePoint;
					//cout << "Maximum value: " << maximumValue << endl;
				}

				//Find minimum
				if ( sampleValue < minimumValue )
				{
					minimumValue = sampleValue;
					minimumPoint = samplePoint;
					//cout << "Minimum value: " << minimumValue << endl;
				}
			}
		}

		//Find the maximum gradient
		vector<double> midPoints;
		string maximumGradientObservable, unit;
		double maximumGradient, lowPoint, midPoint, highPoint;
		for ( int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
		{
			//cout << "Examining: " << doIntegrate[observableIndex] << endl;

			//Find the size of the cell in this observable
			IConstraint * temporaryConstraint = currentCell.GetConstraint( doIntegrate[observableIndex] );
			double cellMaximum = temporaryConstraint->GetMaximum();
			double cellMinimum = temporaryConstraint->GetMinimum();

			//cout << "Cell maximum: " << cellMaximum << endl;
			//cout << "Cell minimum: " << cellMinimum << endl;

			//Find the distance between the maximum and minimum points in this observable
			double pointMaximum = maximumPoint.GetObservable( doIntegrate[observableIndex] )->GetValue();
			double pointMinimum = minimumPoint.GetObservable( doIntegrate[observableIndex] )->GetValue();
			double pointRange = pointMaximum - pointMinimum;

			//cout << "Point maximum: " << pointMaximum << endl;
			//cout << "Point minimum: " << pointMinimum << endl;

			//Find the gradient
			double gradient = abs( ( maximumValue - minimumValue ) * pointRange );//( cellMaximum - cellMinimum ) / pointRange );

			//cout << "Gradient: " << gradient << endl;

			//Store the mid point
			double observableMiddle = cellMinimum + ( ( cellMaximum - cellMinimum ) / 2.0 );
			midPoints.push_back(observableMiddle);

			//Update maximum
			if ( observableIndex == 0 || gradient > maximumGradient )
			{
				//cout << "Is latest maximum" << endl;
				maximumGradient = gradient;
				maximumGradientObservable = doIntegrate[observableIndex];
				unit = temporaryConstraint->GetUnit();
				highPoint = cellMaximum;
				midPoint = observableMiddle;
				lowPoint = cellMinimum;
			}
		}

		//cout << "Maximum gradient found: " << maximumGradient << endl;

		//If the maximum gradient is within tolerance, the cell is finished
		if ( maximumGradient < MAXIMUM_GRADIENT_TOLERANCE )
		{
			//Store the finished cell
			finishedCells.push_back(currentCell);

			//Create a data point at the center of the current cell
			DataPoint cellCenter( InputPoint->GetAllNames() );
			for ( int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
			{
				//Use the mid points for the integrable values
				Observable * temporaryObservable = cellCenter.GetObservable( doIntegrate[observableIndex] );
				temporaryObservable->SetValue( midPoints[observableIndex] );
				cellCenter.SetObservable( doIntegrate[observableIndex], temporaryObservable );
			}
			for ( int observableIndex = 0; observableIndex < dontIntegrate.size(); observableIndex++ )
			{
				//Use given values for unintegrable observables
				Observable * temporaryObservable = InputPoint->GetObservable( dontIntegrate[observableIndex] );
				cellCenter.SetObservable( dontIntegrate[observableIndex], temporaryObservable );
			}

			//Store the center point
			centerPoints.push_back(cellCenter);

			//cout << "Stored finished cell" << endl;
		}
		else
		{
			//Create two cells to replace the current cell
			PhaseSpaceBoundary daughterCell1( currentCell.GetAllNames() );
			PhaseSpaceBoundary daughterCell2( currentCell.GetAllNames() );

			for ( int observableIndex = 0; observableIndex < doIntegrate.size(); observableIndex++ )
			{
				IConstraint * temporaryConstraint = currentCell.GetConstraint( doIntegrate[observableIndex] );
				if ( doIntegrate[observableIndex] == maximumGradientObservable )
				{
					//Split the cells on the observable with the greatest gradient
					if ( lowPoint == midPoint )
					{
						cout << "lowPoint == midPoint" << endl;
						return;
					}
					if ( midPoint == highPoint )
					{
						cout << "midPoint == highPoint" << endl;
						return;
					}
					daughterCell1.SetConstraint( doIntegrate[observableIndex], lowPoint, midPoint, unit );
					daughterCell2.SetConstraint( doIntegrate[observableIndex], midPoint, highPoint, unit );
				}
				else
				{
					//Copy the continuous constraint (if it can be integrated, it must be continuous)
					daughterCell1.SetConstraint( doIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
					daughterCell2.SetConstraint( doIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
				}
			}
			for ( int observableIndex = 0; observableIndex < dontIntegrate.size(); observableIndex++ )
			{
				IConstraint * temporaryConstraint = currentCell.GetConstraint( dontIntegrate[observableIndex] );
				if ( temporaryConstraint->IsDiscrete() )
				{
					//Copy the discrete constraint
					daughterCell1.SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetValues(), temporaryConstraint->GetUnit() );
					daughterCell2.SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetValues(), temporaryConstraint->GetUnit() );
				}
				else
				{
					//Copy the continuous constraint
					daughterCell1.SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
					daughterCell2.SetConstraint( dontIntegrate[observableIndex], temporaryConstraint->GetMinimum(), temporaryConstraint->GetMaximum(), temporaryConstraint->GetUnit() );
				}
			}

			//Add the two new cells to the possibles
			possibleCells.push(daughterCell1);
			possibleCells.push(daughterCell2);

			//cout << "Made new possible cells" << endl;
		}

		//cout << "Number finished: " << finishedCells.size() << endl;
		//cout << "Number possible: " << possibleCells.size() << endl;
		//return;
	}
}

//Destructor
MakeFoam::~MakeFoam()
{
}

//Output all cells to console
void MakeFoam::Debug()
{
	cout << centerPoints.size() << " cells generated" << endl << endl;
	vector<string> allNames = centerPoints[0].GetAllNames();
	vector<string>::iterator observableIterator;
	for ( observableIterator = allNames.begin(); observableIterator != allNames.end(); observableIterator++ )
	{
		cout << *observableIterator << ", ";
	}
	cout << endl;
	for ( int pointIndex = 0; pointIndex < centerPoints.size(); pointIndex++ )
	{
		for ( observableIterator = allNames.begin(); observableIterator != allNames.end(); observableIterator++ )
		{
			cout << centerPoints[pointIndex].GetObservable( *observableIterator )->GetValue() << ", ";
		}
		cout << endl;
	}
}
