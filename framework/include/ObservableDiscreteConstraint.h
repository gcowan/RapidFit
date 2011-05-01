/**
        @class ObservableDiscreteConstraint

        A constraint defining an observable that can take any of a list of discrete values.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef OBSERVABLE_DISCRETE_CONSTRAINT_H
#define OBSERVABLE_DISCRETE_CONSTRAINT_H

//	RapidFit Headers
#include "IConstraint.h"
//	System Headers
#include <string>
#include <vector>

using namespace std;

class ObservableDiscreteConstraint : public IConstraint
{
	public:
		ObservableDiscreteConstraint();
		ObservableDiscreteConstraint( string, vector<double>, string );
		~ObservableDiscreteConstraint();

		void AddValue(double);
		void SetValues( vector<double> );

		//Interface functions
		virtual bool CheckObservable( Observable* );
		virtual Observable* CreateObservable();
		virtual Observable* CreateObservable( TRandom3* );
		virtual string GetUnit();
		virtual double GetMaximum();
		virtual double GetMinimum();
		virtual vector<double> GetValues();
		virtual bool IsDiscrete();

	private:
		vector<double> allValues;
		string unit;
};

#endif
