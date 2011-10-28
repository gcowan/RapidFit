/**
        @class PhysicsBottle

        A collection of PDF-DataSet pairs to be fitted simultaneously with a given parameter set.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

//	RapidFit Headers
#include "ClassLookUp.h"
#include "PhysicsBottle.h"
#include "StringProcessing.h"
//	System Headers
#include <iostream>
#include 	<stdlib.h>

//Default constructor
PhysicsBottle::PhysicsBottle() : allPDFs(), allDataSets(), allConstraints(), bottleParameters(), finalised(false)
{
}

//Constructor with correct argument
PhysicsBottle::PhysicsBottle(ParameterSet * NewParameters) : allPDFs(), allDataSets(), allConstraints(), bottleParameters(NewParameters), finalised(false)
{
}

PhysicsBottle::PhysicsBottle(const PhysicsBottle& newParameters ) : allPDFs(), allDataSets(newParameters.allDataSets), allConstraints(newParameters.allConstraints), bottleParameters(newParameters.bottleParameters), finalised(newParameters.finalised)
{
	for( unsigned int i=0; i< newParameters.allPDFs.size(); ++i )
	{
		allPDFs.push_back( ClassLookUp::CopyPDF(newParameters.allPDFs[i]) );
	}
}

//Destructor
PhysicsBottle::~PhysicsBottle()
{
	while( !allPDFs.empty() )
	{
		if( allPDFs.back() != NULL ) delete allPDFs.back();
		allPDFs.pop_back();
	}
	//cout << "Hello from PhysicsBottle destructor" << endl;
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
		allPDFs.push_back( ClassLookUp::CopyPDF(NewPDF) );
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
		return allPDFs[unsigned(Index)];
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
		return allDataSets[unsigned(Index)];
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
