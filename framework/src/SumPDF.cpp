/**
  @class SumPDF

  An implementation of IPDF for adding the values of two other IPDFs

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

#include "SumPDF.h"
#include "StringProcessing.h"
#include <iostream>

using namespace std;

//Default constructor
SumPDF::SumPDF()
{
}

//Constructor not specifying fraction parameter name
SumPDF::SumPDF( IPDF * FirstPDF, IPDF * SecondPDF, PhaseSpaceBoundary * InputBoundary ) : firstPDF(FirstPDF), secondPDF(SecondPDF), fractionName("FirstPDFFraction"), firstFraction(0.5)
{
	MakePrototypes(InputBoundary);
}

//Constructor specifying fraction parameter name
SumPDF::SumPDF( IPDF * FirstPDF, IPDF * SecondPDF, PhaseSpaceBoundary * InputBoundary, string FractionName ) : firstPDF(FirstPDF), secondPDF(SecondPDF), fractionName(FractionName), firstFraction(0.5)
{
	MakePrototypes(InputBoundary);
}

//Assemble the vectors of parameter/observable names needed
void SumPDF::MakePrototypes( PhaseSpaceBoundary * InputBoundary )
{
	//Make sure the ratio of the two PDFs is included
	vector<string> secondParameterSet = secondPDF->GetPrototypeParameterSet();
	secondParameterSet.push_back(fractionName);

	//Make the prototype parameter set
	prototypeParameterSet = StringProcessing::CombineUniques( firstPDF->GetPrototypeParameterSet(), secondParameterSet );

	//Make the prototype data point
	vector<string> firstObservables = firstPDF->GetPrototypeDataPoint();
	vector<string> secondObservables = secondPDF->GetPrototypeDataPoint();
	prototypeDataPoint = StringProcessing::CombineUniques( firstObservables, secondObservables );

	//Make the do not integrate list
	doNotIntegrateList = StringProcessing::CombineUniques( firstPDF->GetDoNotIntegrateList(), secondPDF->GetDoNotIntegrateList() );

	//Make the correctionss to the integrals for observables unused by only one PDF
	vector<string>::iterator observableIterator;
	IConstraint * inputConstraint;
	firstIntegralCorrection = 1.0;
	secondIntegralCorrection = 1.0;
	for ( observableIterator = firstObservables.begin(); observableIterator != firstObservables.end(); observableIterator++ )
	{
		if ( StringProcessing::VectorContains( &secondObservables, &(*observableIterator) ) == -1 )
		{
			//The first PDF uses this observable, the second doesn't
			inputConstraint = InputBoundary->GetConstraint( *observableIterator );
			bool doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrateList, &(*observableIterator) ) == -1 );
			
			//Update this integral correction
			if ( !inputConstraint->IsDiscrete() && doIntegrate )
			{
				secondIntegralCorrection *= ( inputConstraint->GetMaximum() - inputConstraint->GetMinimum() );
			}
		}
	}
	for ( observableIterator = secondObservables.begin(); observableIterator != secondObservables.end(); observableIterator++ )
	{
		if ( StringProcessing::VectorContains( &firstObservables, &(*observableIterator) ) == -1 )
		{
			//The second PDF uses this observable, the first doesn't
			inputConstraint = InputBoundary->GetConstraint( *observableIterator );
                        bool doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrateList, &(*observableIterator) ) == -1 );

			//Update this integral correction
			if ( !inputConstraint->IsDiscrete() && doIntegrate )
			{
				firstIntegralCorrection *= ( inputConstraint->GetMaximum() - inputConstraint->GetMinimum() );
			}
		}
	}
}

//Destructor
SumPDF::~SumPDF()
{
}

//Indicate whether the function has been set up correctly
bool SumPDF::IsValid()
{
	return firstPDF->IsValid() && secondPDF->IsValid();
}

//Set the function parameters
bool SumPDF::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	PhysicsParameter * newFraction = NewParameterSet->GetPhysicsParameter(fractionName);
	if ( newFraction->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Parameter \"" << fractionName << "\" expected but not found" << endl;
		return false;
	}
	else
	{
		double newFractionValue = newFraction->GetValue();

		//Stupidity check
		if ( newFractionValue > 1.0 || newFractionValue < 0.0 )
		{
			cerr << "Requested impossible fraction: " << newFractionValue << endl;
			return false;
		}
		else
		{
			firstFraction = newFractionValue;
			return firstPDF->SetPhysicsParameters( NewParameterSet ) && secondPDF->SetPhysicsParameters( NewParameterSet );
		}
	}
}

//Return the integral of the function over the given boundary
double SumPDF::Integral( DataPoint* NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	//Get the PDFs' integrals, weighted by firstFraction and corrected for unused observables
	double termOne = firstPDF->Integral( NewDataPoint, NewBoundary ) * firstFraction * firstIntegralCorrection;
	double termTwo = secondPDF->Integral( NewDataPoint, NewBoundary ) * ( 1 - firstFraction ) * secondIntegralCorrection;
	return termOne + termTwo;
}

//Return the function value at the given point
double SumPDF::Evaluate( DataPoint * NewDataPoint )
{
	//Get the PDFs' values, weighted by firstFraction
	double termOne = firstPDF->Evaluate( NewDataPoint ) * firstFraction;
	double termTwo = secondPDF->Evaluate( NewDataPoint ) * ( 1 - firstFraction );
	return termOne + termTwo;
}

//Return a prototype data point
vector<string> SumPDF::GetPrototypeDataPoint()
{
	return prototypeDataPoint;
}

//Return a prototype set of physics parameters
vector<string> SumPDF::GetPrototypeParameterSet()
{
	return prototypeParameterSet;
}

//Return a list of parameters not to be integrated
vector<string> SumPDF::GetDoNotIntegrateList()
{
	return doNotIntegrateList;
}
