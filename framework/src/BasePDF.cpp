/**
        @class BasePDF

        Class that provides a general implementation of IPDF.
        Can inherit from this to make a PDF without worrying about the details.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "BasePDF.h"
#include "ObservableRef.h"
//	System Headers
#include <iostream>

//Constructor
BasePDF::BasePDF() : cachedIntegral(-1.0), cacheValid(false), allParameters(), allObservables(), valid(false), observables()
{
}

//Destructor
BasePDF::~BasePDF()
{
}

//Indicate whether the function has been set up correctly
bool BasePDF::IsValid()
{
	return valid;
}

//Set the function parameters
bool BasePDF::SetPhysicsParameters(ParameterSet * NewParameterSet)
{
	bool success = allParameters.SetPhysicsParameters(NewParameterSet);

	//Invalidate the cache
	cacheValid = false;
	cachedIntegral = -1.0;

	return success;
}

//Return the integral of the function over the given boundary
double BasePDF::Integral(DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary)
{
	if( observables.size() != this->GetPrototypeDataPoint().size() )
	{
		observables=ObservableRef( this->GetPrototypeDataPoint() );
	}
	//Check the boundary is within the correct phase space
	for( unsigned int nameIndex=0; nameIndex < allObservables.size(); ++nameIndex )
	{
		IConstraint * testConstraint = NewBoundary->GetConstraint( observables[nameIndex] );
		if (testConstraint->GetUnit() == "NameNotFoundError")
		{
			cerr << "PDF cannot integrate over phase space: observable \"" << observables[nameIndex].Name().c_str() << "\" not found" << endl;
			return -1.0;
		}
	}

	//Use the cached integral, so long as it makes sense
	if (cacheValid)
	{
		return cachedIntegral;
	}
	else
	{
		//If the PDF normalisation does not depend on the DataPoint, cache the value to save time
		cachedIntegral = Normalisation(NewBoundary);

		//Check the cache is valid
		if ( cachedIntegral > 0.0 )
		{
			cacheValid = true;
			return cachedIntegral;
		}
		else
		{
			//Integral must depend on the DataPoint
			return Normalisation(NewDataPoint, NewBoundary);
		}
	}
	return -1.;
}

//Do the integration
double BasePDF::Normalisation(PhaseSpaceBoundary * NewBoundary)
{
	//	Stupid gcc
	(void)NewBoundary;
	//Just a default value
	return -1.0;
}

//Do the integration
double BasePDF::Normalisation(DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary)
{
	//	Stupid gcc
	(void)NewDataPoint;
	(void)NewBoundary;
	//Just a default value
	return -1.0;
}

/*//Return the function value at the given point
double BasePDF::Evaluate(DataPoint * NewDataPoint)
{
	//Check the datapoint is within the PDF boundary
	if ( xRange.Evaluate(NewDataPoint) > 0.0 )
	{
		return  Value(NewDataPoint);
	}
	else
	{
		cerr << "Datapoint is outside PDF boundary" << endl;
		return 0.0;
	}
}*/

//Calculate the function value
double BasePDF::Evaluate(DataPoint * NewDataPoint)
{
	//	Stupid gcc
	(void)NewDataPoint;
	//Just a default value
	return  -1.0;
}

//Calculate the function value for numerical integration
double BasePDF::EvaluateForNumericIntegral(DataPoint * NewDataPoint)
{
	return this->Evaluate( NewDataPoint ) ;
}

//Calculate the function value
vector<double> BasePDF::EvaluateComponents(DataPoint * NewDataPoint)
{
	//Just assume a single component equal to the ordinary evaluate method
	vector<double> components ;
	components.push_back( this->Evaluate( NewDataPoint ) ) ;
	return  components ;
}

//Return a prototype data point
vector<string> BasePDF::GetPrototypeDataPoint()
{
	return allObservables;
}

//Return a prototype set of physics parameters
vector<string> BasePDF::GetPrototypeParameterSet()
{
	return allParameters.GetAllNames();
}

//Return a list of parameters not to be integrated
vector<string> BasePDF::GetDoNotIntegrateList()
{
	vector<string> emptyVector;
	return emptyVector;
}

// Update integral cache
void BasePDF::UpdateIntegralCache()
{
}

