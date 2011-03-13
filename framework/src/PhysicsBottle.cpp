/**
        @class PhysicsBottle

        A collection of PDF-DataSet pairs to be fitted simultaneously with a given parameter set.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#include "PhysicsBottle.h"
#include <iostream>
#include "StringProcessing.h"
#include <stdlib.h>

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
		exit(1);
	}
	else
	{
		allPDFs.push_back( NewPDF );
		allDataSets.push_back( NewDataSet );
	}
}

//Store a ConstraintFunction
void PhysicsBottle::AddConstraint( ConstraintFunction * NewConstraint )
{
	allConstraints.push_back(NewConstraint);
}
vector< ConstraintFunction* > PhysicsBottle::GetConstraints()
{
	return allConstraints;
}

//Retrieve the number of results stored
int PhysicsBottle::NumberResults()
{
	return int(allPDFs.size());
}

//Retrieve the PDF corresponding to a particular result number
IPDF * PhysicsBottle::GetResultPDF(int Index)
{
	if ( Index < int(allPDFs.size()) )
	{
		return allPDFs[Index];
	}
	else
	{
		cerr << "PDF index (" << Index << ") out of range in PhysicsBottle" << endl;
		exit(1);
	}
}

//Retrieve the data set corresponding to a particular result number
IDataSet * PhysicsBottle::GetResultDataSet(int Index)
{
	if ( Index < int(allDataSets.size()) )
	{
		return allDataSets[ Index ];
	}
	else
	{
		cerr << "DataSet index (" << Index << ") out of range in PhysicsBottle" << endl;
		exit(1);
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
		for (unsigned int pdfIndex = 0; pdfIndex < allPDFs.size(); ++pdfIndex)
		{
			allPDFs[pdfIndex]->SetPhysicsParameters( bottleParameters );
		}

		return true;
	}
	else
	{
		cerr << "Bottle parameters not successfully updated" << endl;
		exit(1);
	}
}

//Make the bottle read only, and cull unused parameters: "Seal the bottle"
void PhysicsBottle::Finalise()
{
	if (finalised)
	{
		cout << "Bottle already finalised" << endl;
	}
	else
	{
		//Create a unique vector of parameters that are used
		vector<string> usedParameters;
		for (unsigned int pdfIndex = 0; pdfIndex < allPDFs.size(); ++pdfIndex )
		{
			usedParameters = StringProcessing::CombineUniques( usedParameters, allPDFs[pdfIndex]->GetPrototypeParameterSet() );
		}

		//Create a new parameter set, containing only those that are used
		ParameterSet * culledParameters = new ParameterSet(usedParameters);
		for (unsigned int usedIndex = 0; usedIndex < usedParameters.size(); ++usedIndex )
		{
			PhysicsParameter * usedParameter = bottleParameters->GetPhysicsParameter( usedParameters[usedIndex] );
			culledParameters->SetPhysicsParameter( usedParameters[usedIndex], usedParameter );
		}

		//Replace the bottle's parameter set pointer
		bottleParameters = culledParameters;

		finalised = true;
	}
}
