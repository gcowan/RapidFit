/**
        @class PhysicsBottle

        A collection of PDF-DataSet pairs to be fitted simultaneously with a given parameter set.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef PHYSICS_BOTTLE_H
#define PHYSICS_BOTTLE_H

//	RapidFit Headers
#include "IPDF.h"
#include "IDataSet.h"
#include "ParameterSet.h"
#include "ConstraintFunction.h"
//	System Headers
#include <vector>

class PhysicsBottle
{
	public:
		PhysicsBottle();
		PhysicsBottle( ParameterSet* );
		~PhysicsBottle();
		PhysicsBottle(const PhysicsBottle&);

		void AddResult( IPDF*, IDataSet* );
		void AddConstraint( ConstraintFunction* );
		int NumberResults();
		IPDF* GetResultPDF(int);
		IDataSet* GetResultDataSet(int);
		vector< ConstraintFunction* > GetConstraints();
		ParameterSet * GetParameterSet();
		bool SetParameterSet( ParameterSet* );
		void Finalise();

	private:
		//	Uncopyable!
		PhysicsBottle& operator = ( const PhysicsBottle& );

		vector< IPDF* > allPDFs;
		vector< IDataSet* > allDataSets;
		vector< ConstraintFunction* > allConstraints;
		ParameterSet * bottleParameters;
		bool finalised;
};

#endif
