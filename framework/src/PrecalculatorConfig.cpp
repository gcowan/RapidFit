

#include "PrecalculatorConfig.h"
#include "IPrecalculator.h"
#include "ClassLookUp.h"
#include "FitResult.h"

#include <string>

using namespace::std;

PrecalculatorConfig::PrecalculatorConfig() : name("undefined"), weightName("undefined"), filename("WeightedFile"), config(0), useAlpha(false)
{
}

void PrecalculatorConfig::SetCalculatorName( string input )
{
	name = input;
}

void PrecalculatorConfig::SetWeightName( string input )
{
	weightName = input;
}

void PrecalculatorConfig::SetConfig( unsigned int input )
{
	config = input;
}

void PrecalculatorConfig::SetAlpha( bool input )
{
	useAlpha = input;
}

void PrecalculatorConfig::SetFileName( string input )
{
	size_t found = input.find(".root");
	if( found == string::npos )
	{
		filename = input;
	}
	else
	{
		filename = input.substr(0,(input.length()-5));
	}
}

string PrecalculatorConfig::GetFileName()
{
	return filename;
}

IPrecalculator* PrecalculatorConfig::GetPreCalculator( FitResult* inputResult )
{
	IPrecalculator* thisCalculator = ClassLookUp::LookUpPrecalculator( name, weightName, inputResult, config );
	thisCalculator->SetApplyAlphaCorrection( useAlpha );
	return thisCalculator;
}

