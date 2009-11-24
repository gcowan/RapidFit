/**
        @class PhysicsBottle

        A collection of PDF-DataSet pairs to be fitted simultaneously with a given parameter set.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef PHYSICS_BOTTLE_H
#define PHYSICS_BOTTLE_H

#include "IPDF.h"
#include "IDataSet.h"
#include "ParameterSet.h"
#include <vector>

class PhysicsBottle
{
	public:
		PhysicsBottle();
		PhysicsBottle( ParameterSet* );
		~PhysicsBottle();

		void AddResult( IPDF*, IDataSet* );
		int NumberResults();
		IPDF * GetResultPDF(int);
		IDataSet * GetResultDataSet(int);
		ParameterSet * GetParameterSet();
		bool SetParameterSet( ParameterSet* );
		void Finalise();

	private:
		vector< IPDF* > allPDFs;
		vector< IDataSet* > allDataSets;
		ParameterSet * bottleParameters;
		bool finalised;
};

#endif
