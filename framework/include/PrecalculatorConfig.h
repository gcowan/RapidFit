
#pragma once
#ifndef RAPIDFIT_PRECALC_H
#define RAPIDFIT_PRECALC_H

#include "IPrecalculator.h"
#include "FitResult.h"
#include <string>
#include <cmath>

using namespace::std;

class PrecalculatorConfig
{

	public:
		PrecalculatorConfig();

		void SetCalculatorName( string );
		void SetWeightName( string );
		void SetConfig( unsigned int );
		void SetFileName( string );
		void SetAlpha( bool );
		string GetFileName();

		IPrecalculator* GetPreCalculator( FitResult* input );

	private:
		string name;
		string weightName;
		string filename;
		unsigned int config;
		bool useAlpha;
};

#endif

