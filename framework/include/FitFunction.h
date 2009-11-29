/**
        @class FitFunction

        Parent class for the function to minimise
        Overload the evaluate methods and UP value for Chi2, NLL, etc.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef FIT_FUNCTION_H
#define FIT_FUNCTION_H

#include "PhysicsBottle.h"
#include "RapidFitIntegrator.h"

class FitFunction
{
	public:
		FitFunction();
		~FitFunction();

		void SetPhysicsBottle( PhysicsBottle* );
		PhysicsBottle * GetPhysicsBottle();
		bool SetParameterSet( ParameterSet* );
		ParameterSet * GetParameterSet();
		double Evaluate();
		void Finalise();
		void UseEventWeights(string);

		//Overload this function in child classes
		virtual double UpErrorValue(int);

	protected:
		//Overload these functions in child classes
		virtual double EvaluateDataSet( IPDF*, IDataSet*, RapidFitIntegrator*, int );
		virtual double EvaluateParameterSet( ParameterSet*, vector<string> );
		virtual void Precalculation();

		PhysicsBottle * allData;
		vector<string> interestingParameters;
		vector<string>::iterator nameIterator;
		vector< RapidFitIntegrator* > allIntegrators;
		double testDouble;
		bool useWeights;
		string weightObservableName;
};

#endif
