/**
  @class FitAssembler

  The intention is for this class to formalise the process of assembling the components of a fit
  Ideally it will be a set of nested static methods, starting from more and more rudimentary components

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

#include "FitAssembler.h"
#include "ClassLookUp.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

//The final stage - do the minimisation
FitResult * FitAssembler::DoFit( IMinimiser * Minimiser, FitFunction * TheFunction )
{
	Minimiser->Minimise(TheFunction);
	return Minimiser->GetFitResult();
}

//Create the minimiser and fit function
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, PhysicsBottle * Bottle )
{
	IMinimiser * minimiser = MinimiserConfig->GetMinimiser( Bottle->GetParameterSet()->GetAllNames().size() );
	FitFunction * theFunction = FunctionConfig->GetFitFunction();
	theFunction->SetPhysicsBottle(Bottle);

	FitResult * result = DoFit( minimiser, theFunction );

	delete theFunction;
	delete minimiser;
	return result;
}

//Create the physics bottle
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,
		vector< PDFWithData* > BottleData, vector< ConstraintFunction* > BottleConstraints )
{
	PhysicsBottle * bottle = new PhysicsBottle(BottleParameters);

	//Fill the bottle - data generation occurs in this step
	for ( int resultIndex = 0; resultIndex < BottleData.size(); resultIndex++ )
	{
		BottleData[resultIndex]->SetPhysicsParameters(BottleParameters);
		bottle->AddResult( BottleData[resultIndex]->GetPDF(), BottleData[resultIndex]->GetDataSet() );
	}

	//Add the constraints
	for ( int constraintIndex = 0; constraintIndex < BottleConstraints.size(); constraintIndex++ )
	{
		bottle->AddConstraint( BottleConstraints[constraintIndex] );
	}

	bottle->Finalise();
	FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );

	delete bottle;
	return result;
}

//Create the physics bottle with pre-made data
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters,
		vector< IPDF* > AllPDFs, vector< IDataSet* > AllData, vector< ConstraintFunction* > BottleConstraints )
{
	if ( AllPDFs.size() == AllData.size() )
	{
		PhysicsBottle * bottle = new PhysicsBottle(BottleParameters);

		//Fill the bottle - data already generated
		for ( int resultIndex = 0; resultIndex < AllData.size(); resultIndex++ )
		{
			AllPDFs[resultIndex]->SetPhysicsParameters(BottleParameters);
			bottle->AddResult( AllPDFs[resultIndex], AllData[resultIndex] );
		}

		//Add the constraints
		for ( int constraintIndex = 0; constraintIndex < BottleConstraints.size(); constraintIndex++ )
		{
			bottle->AddConstraint( BottleConstraints[constraintIndex] );
		}  

		bottle->Finalise();
		FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );

		delete bottle;
		return result;
	}
	else
	{
		cerr << "Mismatched number of PDFs and DataSets" << endl;
		exit(1);
	}
}
