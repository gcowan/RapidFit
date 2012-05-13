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
#include <stdlib.h>

//Constructor with correct argument
PhysicsBottle::PhysicsBottle( ParameterSet * NewParameters ) : allPDFs(), allDataSets(), allConstraints(), bottleParameters(new ParameterSet( *NewParameters) )
{
}

PhysicsBottle::PhysicsBottle(const PhysicsBottle& newParameters ) :
	allPDFs(), allDataSets(newParameters.allDataSets), allConstraints(newParameters.allConstraints), bottleParameters()
{
	for( unsigned int i=0; i< newParameters.allPDFs.size(); ++i )
	{
		allPDFs.push_back( ClassLookUp::CopyPDF( newParameters.allPDFs[i] ) );
	}
	if( newParameters.bottleParameters == NULL ) bottleParameters = NULL;
	else bottleParameters = new ParameterSet( *(newParameters.bottleParameters) );
}

//Destructor
PhysicsBottle::~PhysicsBottle()
{
	while( !allPDFs.empty() )
	{
		if( allPDFs.back() != NULL ) delete allPDFs.back();
		allPDFs.pop_back();
	}
	if( bottleParameters != NULL ) delete bottleParameters;
}

//Store a PDF/dataset pair
void PhysicsBottle::AddResult( IPDF * NewPDF, IDataSet * NewDataSet )
{
	allPDFs.push_back( ClassLookUp::CopyPDF(NewPDF) );
	allDataSets.push_back( NewDataSet );
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
	bottleParameters->SetPhysicsParameters( NewParameters );

	//Propagate the change to all stored PDFs
	for (unsigned int pdfIndex = 0; pdfIndex < allPDFs.size(); ++pdfIndex)
	{
		allPDFs[pdfIndex]->UpdatePhysicsParameters( bottleParameters );
		allPDFs[pdfIndex]->UnsetCache();
	}

	return true;
}


