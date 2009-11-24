/**
        @class BasePDF

        Class that provides a general implementation of IPDF.
        Can inherit from this to make a PDF without worrying about the details.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#include "BasePDF.h"
#include <iostream>

//Constructor
BasePDF::BasePDF() : valid(false), cacheValid(false), cachedIntegral(-1.0)
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
	//Check the boundary is within the correct phase space
	vector<string>::iterator nameIterator;
	for(nameIterator = allObservables.begin(); nameIterator != allObservables.end(); nameIterator++)
	{
		IConstraint * testConstraint = NewBoundary->GetConstraint( *nameIterator );
		if (testConstraint->GetUnit() == "NameNotFoundError")
		{
			cerr << "PDF cannot integrate over phase space: observable \"" << *nameIterator << "\" not found" << endl;
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
}

//Do the integration
double BasePDF::Normalisation(PhaseSpaceBoundary * NewBoundary)
{
	//Just a default value
	return -1.0;
}

//Do the integration
double BasePDF::Normalisation(DataPoint * NewDataPoint, PhaseSpaceBoundary * NewBoundary)
{
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
	//Just a default value
	return  -1.0;
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
