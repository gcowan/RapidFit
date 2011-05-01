/**
        @class FitFunction

        Parent class for the function to minimise
        Overload the evaluate methods and UP value for Chi2, NLL, etc.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef FIT_FUNCTION_H
#define FIT_FUNCTION_H

//	RapidFit Headers
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
		//	Uncopyable!
		FitFunction ( const FitFunction& );
		FitFunction& operator = ( const FitFunction& );
		//Overload these functions in child classes
		virtual double EvaluateDataSet( IPDF*, IDataSet*, RapidFitIntegrator* );

		PhysicsBottle * allData;
		vector< RapidFitIntegrator* > allIntegrators;
		double testDouble;
		bool useWeights;
		string weightObservableName;
};

#endif
