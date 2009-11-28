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
FitResult * FitAssembler::DoFit( MinimiserConfiguration * MinimiserConfig, FitFunctionConfiguration * FunctionConfig, ParameterSet * BottleParameters, vector< PDFWithData* > BottleData )
{
	PhysicsBottle * bottle = new PhysicsBottle(BottleParameters);

	//Fill the bottle - data generation occurs in this step
	for ( int resultIndex = 0; resultIndex < BottleData.size(); resultIndex++ )
	{
		BottleData[resultIndex]->SetPhysicsParameters(BottleParameters);
		bottle->AddResult( BottleData[resultIndex]->GetPDF(), BottleData[resultIndex]->GetDataSet() );
	}

	bottle->Finalise();
	FitResult * result = DoFit( MinimiserConfig, FunctionConfig, bottle );

	delete bottle;
	return result;
}
