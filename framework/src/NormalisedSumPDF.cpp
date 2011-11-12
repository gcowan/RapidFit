/**
  @class NormalisedSumPDF

  An implementation of IPDF for adding the values of two other IPDFs, normalised relative to each other

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-11-12
  */

//	RapidFit Headers
#include "ClassLookUp.h"
#include "NormalisedSumPDF.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include <math.h>
#include <cstdlib>

pthread_mutex_t pdf_lock_1, pdf_lock_2;

using namespace std;

//Default constructor
NormalisedSumPDF::NormalisedSumPDF() : prototypeDataPoint(), prototypeParameterSet(), doNotIntegrateList(), firstPDF(NULL), secondPDF(NULL), firstIntegrator(NULL), secondIntegrator(NULL), firstFraction(), firstIntegralCorrection(), secondIntegralCorrection(), fractionName(), integrationBoundary(NULL)
{
}

//Constructor not specifying fraction parameter name
NormalisedSumPDF::NormalisedSumPDF( IPDF * FirstPDF, IPDF * SecondPDF, PhaseSpaceBoundary * InputBoundary ) : prototypeDataPoint(), prototypeParameterSet(), doNotIntegrateList(), firstPDF( ClassLookUp::CopyPDF(FirstPDF) ), secondPDF( ClassLookUp::CopyPDF(SecondPDF) ), firstIntegrator( new RapidFitIntegrator(FirstPDF) ), secondIntegrator( new RapidFitIntegrator(SecondPDF) ), firstFraction(0.5), firstIntegralCorrection(), secondIntegralCorrection(), fractionName("FirstPDFFraction"), integrationBoundary(InputBoundary)
{
	firstIntegrator->SetPDF( firstPDF );
	secondIntegrator->SetPDF( secondPDF );
	MakePrototypes(InputBoundary);
}

NormalisedSumPDF::NormalisedSumPDF( const NormalisedSumPDF& input ) : BasePDF( (BasePDF) input ),
	prototypeDataPoint( input.prototypeDataPoint ), prototypeParameterSet( input.prototypeParameterSet ), firstPDF( ClassLookUp::CopyPDF( input.firstPDF ) ),
	secondPDF( ClassLookUp::CopyPDF( input.secondPDF ) ), doNotIntegrateList( StringProcessing::CombineUniques( input.firstPDF->GetDoNotIntegrateList(), input.secondPDF->GetDoNotIntegrateList() ) ),
	firstIntegrator( new RapidFitIntegrator( *(input.firstIntegrator) ) ), secondIntegrator( new RapidFitIntegrator( *(input.secondIntegrator) ) ),
	firstFraction( input.firstFraction ), firstIntegralCorrection( input.firstIntegralCorrection ), secondIntegralCorrection( input.secondIntegralCorrection ),
	fractionName( input.fractionName ), integrationBoundary( input.integrationBoundary )
{
	firstIntegrator->SetPDF( firstPDF );
	secondIntegrator->SetPDF( secondPDF );

	ParameterSet tempSet = input.allParameters;
	if( ! firstPDF->GetActualParameterSet()->GetAllNames().empty() ) firstPDF->SetPhysicsParameters( &tempSet );
	if( ! secondPDF->GetActualParameterSet()->GetAllNames().empty() ) secondPDF->SetPhysicsParameters( &tempSet );
}

//Constructor specifying fraction parameter name
NormalisedSumPDF::NormalisedSumPDF( IPDF * FirstPDF, IPDF * SecondPDF, PhaseSpaceBoundary * InputBoundary, string FractionName ) : prototypeDataPoint(), prototypeParameterSet(), doNotIntegrateList(), firstPDF(FirstPDF), secondPDF(SecondPDF), firstIntegrator( new RapidFitIntegrator(FirstPDF) ), secondIntegrator( new RapidFitIntegrator(SecondPDF) ), firstFraction(0.5), firstIntegralCorrection(), secondIntegralCorrection(), fractionName(FractionName), integrationBoundary(InputBoundary)
{
	MakePrototypes(InputBoundary);
}

//Assemble the vectors of parameter/observable names needed
void NormalisedSumPDF::MakePrototypes( PhaseSpaceBoundary * InputBoundary )
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
	IConstraint * inputConstraint=NULL;
	firstIntegralCorrection = 1.0;
	secondIntegralCorrection = 1.0;
	bool doIntegrate=false;
	for ( observableIterator = firstObservables.begin(); observableIterator != firstObservables.end(); ++observableIterator )
	{
		if ( StringProcessing::VectorContains( &secondObservables, &(*observableIterator) ) == -1 )
		{
			//The first PDF uses this observable, the second doesn't
			inputConstraint = InputBoundary->GetConstraint( *observableIterator );
			doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrateList, &(*observableIterator) ) == -1 );

			//Update this integral correction
			if ( !inputConstraint->IsDiscrete() && doIntegrate )
			{
				secondIntegralCorrection *= ( inputConstraint->GetMaximum() - inputConstraint->GetMinimum() );
			}
		}
	}
	for ( observableIterator = secondObservables.begin(); observableIterator != secondObservables.end(); ++observableIterator )
	{
		if ( StringProcessing::VectorContains( &firstObservables, &(*observableIterator) ) == -1 )
		{
			//The second PDF uses this observable, the first doesn't
			inputConstraint = InputBoundary->GetConstraint( *observableIterator );
			doIntegrate = ( StringProcessing::VectorContains( &doNotIntegrateList, &(*observableIterator) ) == -1 );

			//Update this integral correction
			if ( !inputConstraint->IsDiscrete() && doIntegrate )
			{
				firstIntegralCorrection *= ( inputConstraint->GetMaximum() - inputConstraint->GetMinimum() );
			}
		}
	}
}

//Destructor
NormalisedSumPDF::~NormalisedSumPDF()
{
	//cout << "Hello from Normalised destructor" << endl;
	if( firstPDF != NULL ) delete firstPDF;
	if( secondPDF != NULL ) delete secondPDF;
	if( firstIntegrator != NULL ) delete firstIntegrator;
	if( secondIntegrator != NULL ) delete secondIntegrator;
}

//Indicate whether the function has been set up correctly
bool NormalisedSumPDF::IsValid()
{
	return firstPDF->IsValid() && secondPDF->IsValid();
}

//Set the function parameters
bool NormalisedSumPDF::SetPhysicsParameters( ParameterSet * NewParameterSet )
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
			bool output = ( firstPDF->SetPhysicsParameters( NewParameterSet ) && secondPDF->SetPhysicsParameters( NewParameterSet ) );
			output = output && allParameters.AddPhysicsParameters( NewParameterSet );
			return output;
		}
	}
}

//Return the integral of the function over the given boundary
double NormalisedSumPDF::Integral( DataPoint* NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	//	Stupid gcc
	(void)NewDataPoint;
	(void)NewBoundary;
	//The evaluate method already returns a normalised value
	return 1.0;
}

//Return the function value at the given point
double NormalisedSumPDF::Evaluate( DataPoint * NewDataPoint )
{
	//Calculate the integrals of the PDFs
	double firstIntegral = firstIntegrator->Integral( NewDataPoint, integrationBoundary, true ) * firstIntegralCorrection;
	double secondIntegral = secondIntegrator->Integral( NewDataPoint, integrationBoundary, true ) * secondIntegralCorrection;

	//Get the PDFs' values, normalised and weighted by firstFrsction
	double termOne = ( firstPDF->Evaluate( NewDataPoint ) * firstFraction ) / firstIntegral;
	double termTwo = ( secondPDF->Evaluate( NewDataPoint ) * ( 1 - firstFraction ) ) / secondIntegral;

	//Return the sum
	return termOne + termTwo;
}

//Return the function value at the given point
double NormalisedSumPDF::EvaluateForNumericIntegral( DataPoint * NewDataPoint )
{
	//Calculate the integrals of the PDFs
	double firstIntegral = firstIntegrator->Integral( NewDataPoint, integrationBoundary, true ) * firstIntegralCorrection;
	double secondIntegral = secondIntegrator->Integral( NewDataPoint, integrationBoundary, true ) * secondIntegralCorrection;
	
	//Get the PDFs' values, normalised and weighted by firstFrsction
	double termOne = ( firstPDF->EvaluateForNumericIntegral( NewDataPoint ) * firstFraction ) / firstIntegral;
	double termTwo = ( secondPDF->EvaluateForNumericIntegral( NewDataPoint ) * ( 1 - firstFraction ) ) / secondIntegral;
	
	//Return the sum
	return termOne + termTwo;
}

//Return the function value at the given point
vector<double> NormalisedSumPDF::EvaluateComponents( DataPoint * NewDataPoint )
{
	//Calculate the integrals of the PDFs
	double firstIntegral = firstIntegrator->Integral( NewDataPoint, integrationBoundary, true ) * firstIntegralCorrection;
	double secondIntegral = secondIntegrator->Integral( NewDataPoint, integrationBoundary, true ) * secondIntegralCorrection;

	//Get the components of each term
	vector<double> termOneComponents = firstPDF->EvaluateComponents( NewDataPoint ) ;
	vector<double> termTwoComponents = secondPDF->EvaluateComponents( NewDataPoint );

	//Insert components in output vector with correct weights.	
	vector<double> components ;
	for(unsigned int ii=0; ii<termOneComponents.size(); ++ii ) components.push_back( termOneComponents[ii]*firstFraction/firstIntegral ) ;
	for(unsigned int ii=0; ii<termTwoComponents.size(); ++ii ) components.push_back( termTwoComponents[ii]*(1.-firstFraction)/secondIntegral ) ;

	// Return the complete set of components
	return components;
}


//Return a prototype data point
vector<string> NormalisedSumPDF::GetPrototypeDataPoint()
{
	return prototypeDataPoint;
}

//Return a prototype set of physics parameters
vector<string> NormalisedSumPDF::GetPrototypeParameterSet()
{
	return prototypeParameterSet;
}

//Return a list of parameters not to be integrated
vector<string> NormalisedSumPDF::GetDoNotIntegrateList()
{
	return doNotIntegrateList;
}

// Update the integral cache for the two RapidFitIntegrators
void NormalisedSumPDF::UpdateIntegralCache()
{
	if( integrationBoundary !=NULL )
	{
		firstIntegrator->UpdateIntegralCache(integrationBoundary);
		secondIntegrator->UpdateIntegralCache(integrationBoundary);
	}
}
