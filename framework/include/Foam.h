/**
        @class Foam

        Class for generating toy data from a PDF.
        Just a wrapper for the Root TFoam generator

 	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-11-10
*/


#ifndef FOAM_H
#define FOAM_H

#include "IDataGenerator.h"
#include "MemoryDataSet.h"
#include "TRandom3.h"
#include "IntegratorFunction.h"
#include "TFoam.h"

class Foam : public IDataGenerator
{
	public:
		Foam();
		Foam( PhaseSpaceBoundary*, IPDF* );
		~Foam();

		//Interface functions
		virtual int GenerateData(int);
		virtual IDataSet * GetDataSet();

	protected:
		IntegratorFunction * generationFunction;
		PhaseSpaceBoundary * generationBoundary;
		MemoryDataSet * newDataSet;
		TRandom3 * rootRandom;
		vector< TFoam* > foamGenerators;
		int dataNumber;
		vector<string> allNames, discreteNames, continuousNames;
		vector< vector<double> > discreteValues;
		vector<double> minima, ranges;
};

#endif
