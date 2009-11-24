/**
        @class PhysicsBottle

        A collection of PDF-DataSet pairs to be fitted simultaneously with a given parameter set.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "PhysicsBottle.h"
#include <iostream>
#include "InvalidObject.h"
#include "StringProcessing.h"

//Default constructor
PhysicsBottle::PhysicsBottle() : finalised(false)
{
}

//Constructor with correct argument
PhysicsBottle::PhysicsBottle(ParameterSet * NewParameters) : bottleParameters(NewParameters), finalised(false)
{
}

//Destructor
PhysicsBottle::~PhysicsBottle()
{
}

//Store a PDF/dataset pair
void PhysicsBottle::AddResult( IPDF * NewPDF, IDataSet * NewDataSet )
{
	if ( finalised )
	{
		cerr << "Bottle finalised - cannot add result" << endl;
	}
	else
	{
		allPDFs.push_back( NewPDF );
		allDataSets.push_back( NewDataSet );
	}
}

//Retrieve the number of results stored
int PhysicsBottle::NumberResults()
{
	return allPDFs.size();
}

//Retrieve the PDF corresponding to a particular result number
IPDF * PhysicsBottle::GetResultPDF(int Index)
{
	if ( Index < allPDFs.size() )
	{
		return allPDFs[Index];
	}
	else
	{
		cerr << "Index (" << Index << ") out of range in PhysicsBottle" << endl;
		return new InvalidObject( "Index out of range in PhysicsBottle" );
	}
}

//Retrieve the data set corresponding to a particular result number
IDataSet * PhysicsBottle::GetResultDataSet(int Index)
{
	if ( Index < allDataSets.size() )
	{
		return allDataSets[ Index ];
	}
	else
	{
		cerr << "Index (" << Index << ") out of range in PhysicsBottle" << endl;
		return new InvalidObject( "Index out of range in PhysicsBottle" );
	}

}

//Retrieve the parameter set
ParameterSet * PhysicsBottle::GetParameterSet()
{
	return bottleParameters;
}

//Change the parameter values
bool PhysicsBottle::SetParameterSet(ParameterSet * NewParameters)
{
	if ( bottleParameters->SetPhysicsParameters( NewParameters ) )
	{
		//Propagate the change to all stored PDFs
		for (int pdfIndex = 0; pdfIndex < allPDFs.size(); pdfIndex++)
		{
			allPDFs[pdfIndex]->SetPhysicsParameters( bottleParameters );
		}

		return true;
	}
	else
	{
		cerr << "Bottle parameters not successfully updated" << endl;
		return false;
	}
}

//Make the bottle read only, and cull unused parameters: "Seal the bottle"
void PhysicsBottle::Finalise()
{
	if (finalised)
	{
		cerr << "Bottle already finalised" << endl;
	}
	else
	{
		//Create a unique vector of parameters that are used
		vector<string> usedParameters;
		for ( int pdfIndex = 0; pdfIndex < allPDFs.size(); pdfIndex++ )
		{
			usedParameters = StringProcessing::CombineUniques( usedParameters, allPDFs[pdfIndex]->GetPrototypeParameterSet() );
		}

		//Create a new parameter set, containing only those that are used
		ParameterSet * culledParameters = new ParameterSet(usedParameters);
		for ( int usedIndex = 0; usedIndex < usedParameters.size(); usedIndex++ )
		{
			PhysicsParameter * usedParameter = bottleParameters->GetPhysicsParameter( usedParameters[usedIndex] );
			culledParameters->SetPhysicsParameter( usedParameters[usedIndex], usedParameter );
		}

		//Replace the bottle's parameter set pointer
		bottleParameters = culledParameters;

		finalised = true;
	}
}
