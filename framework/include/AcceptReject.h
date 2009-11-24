/**
        @class AcceptReject

        Class for generating toy data from a PDF.
        Can inherit from this to implement preselection for a particular PDF.

 	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef ACCEPT_REJECT_H
#define ACCEPT_REJECT_H

#include "IDataGenerator.h"
#include "MemoryDataSet.h"
#include "TRandom3.h"
#include "IPDF.h"

class AcceptReject : public IDataGenerator
{
	public:
		AcceptReject();
		AcceptReject( PhaseSpaceBoundary*, IPDF* );
		~AcceptReject();

		//Interface functions
		virtual int GenerateData(int);
		virtual IDataSet * GetDataSet();

	protected:
		virtual bool Preselection( DataPoint*, double );

		IPDF * generationFunction;
		PhaseSpaceBoundary * generationBoundary;
		int dataNumber;
		MemoryDataSet * newDataSet;
		TRandom3 * rootRandom;
		double moreThanMaximum;
		int numberAttempts;
};

#endif
