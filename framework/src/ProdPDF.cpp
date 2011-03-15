/**
        @class ProdPDF

        An implementation of IPDF for multiplying the values of two other IPDFs

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "ProdPDF.h"
#include <iostream>
#include "StringProcessing.h"

using namespace std;

//Default constructor
ProdPDF::ProdPDF()
{
}

//Constructor not specifying fraction parameter name
ProdPDF::ProdPDF( IPDF * FirstPDF, IPDF * SecondPDF ) : firstPDF(FirstPDF), secondPDF(SecondPDF)
{
	MakePrototypes();
}

//Assemble the vectors of parameter/observable names needed
void ProdPDF::MakePrototypes()
{
	prototypeParameterSet = StringProcessing::CombineUniques( firstPDF->GetPrototypeParameterSet(), secondPDF->GetPrototypeParameterSet() );
	prototypeDataPoint = StringProcessing::CombineUniques( firstPDF->GetPrototypeDataPoint(), secondPDF->GetPrototypeDataPoint() );
	doNotIntegrateList = StringProcessing::CombineUniques( firstPDF->GetDoNotIntegrateList(), secondPDF->GetDoNotIntegrateList() );
}

//Destructor
ProdPDF::~ProdPDF()
{
}

//Indicate whether the function has been set up correctly
bool ProdPDF::IsValid()
{
	return firstPDF->IsValid() && secondPDF->IsValid();
}

//Set the function parameters
bool ProdPDF::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	return firstPDF->SetPhysicsParameters( NewParameterSet ) && secondPDF->SetPhysicsParameters( NewParameterSet );
}

//Return the integral of the function over the given boundary
double ProdPDF::Integral( DataPoint* NewDataPoint, PhaseSpaceBoundary * NewBoundary )
{
	//Note that this is almost certainly wrong. However, I don't know a good analytical solution.
	//In cases that the formula is incorrect, it will be caught by the numerical integration check.
	double termOne = firstPDF->Integral( NewDataPoint, NewBoundary );
	double termTwo = secondPDF->Integral( NewDataPoint, NewBoundary );
	return termOne * termTwo;
}

//Return the function value at the given point
double ProdPDF::Evaluate( DataPoint * NewDataPoint )
{
	double termOne = firstPDF->Evaluate( NewDataPoint );
	double termTwo = secondPDF->Evaluate( NewDataPoint );
	return termOne * termTwo;
}

//Return the function value at the given point
vector<double> ProdPDF::EvaluateComponents( DataPoint * NewDataPoint )
{

	//Get the components of each term
	vector<double> termOneComponents = firstPDF->EvaluateComponents( NewDataPoint ) ;
	vector<double> termTwoComponents = secondPDF->EvaluateComponents( NewDataPoint );
	
	//Insert components in output vector with correct weights.	
	vector<double> components ;
	for(unsigned int ii=0; ii<termOneComponents.size(); ++ii ) {
		for(unsigned int jj=0; jj<termTwoComponents.size(); ++jj ) {
			components.push_back( termOneComponents[ii]*termTwoComponents[jj] ) ;
		}
	}

	// Return the complete set of components
	return components;
}


//Return a prototype data point
vector<string> ProdPDF::GetPrototypeDataPoint()
{
	return prototypeDataPoint;
}

//Return a prototype set of physics parameters
vector<string> ProdPDF::GetPrototypeParameterSet()
{
	return prototypeParameterSet;
}

//Return a list of parameters not to be integrated
vector<string> ProdPDF::GetDoNotIntegrateList()
{
        return doNotIntegrateList;
}

// Update integral cache
void ProdPDF::UpdateIntegralCache()
{
}
