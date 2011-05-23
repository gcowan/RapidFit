/**
        @class Foam

        Class for generating toy data from a PDF.
        Just a wrapper for the Root TFoam generator

 	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-11-10
*/


#ifndef FOAM_H
#define FOAM_H

//	ROOT Headers
#include "TRandom3.h"
#include "TFoam.h"
#include "TFile.h"
//	RapidFit Headers
#include "IDataGenerator.h"
#include "MemoryDataSet.h"
#include "IntegratorFunction.h"
//	System Headers
#include <vector>

class Foam : public IDataGenerator
{
	public:
		Foam();
		Foam( PhaseSpaceBoundary*, IPDF* );
		virtual ~Foam();

		//Interface functions
		virtual int GenerateData(int);
		virtual IDataSet * GetDataSet();

	protected:
		//	Uncopyable!
		Foam ( const Foam& );
		Foam& operator = ( const Foam& );

		void Init();
		void RemoveGenerator();

		vector<TFile*> Open_Files;
		IPDF * InputPDF;
		IntegratorFunction * generationFunction;
		PhaseSpaceBoundary * generationBoundary;
		MemoryDataSet * newDataSet;
		TRandom3 * rootRandom;
		vector< TFoam* > foamGenerators;
		vector< IntegratorFunction* > storedIntegrator;
		vector< DataPoint* > storedDatapoint;
		int dataNumber;
		vector< vector<double> > discreteCombinations;
		vector<string> allNames, discreteNames, continuousNames;
		vector< vector<double> > discreteValues;
		vector<double> minima, ranges;
};

#endif
