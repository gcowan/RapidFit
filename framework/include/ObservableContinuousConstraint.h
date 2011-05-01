/**
        @class ObservableContinuousConstraint

        A constraint defining an observable that can take any value between a given maximum and minimum.

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef OBSERVABLE_CONTINUOUS_CONSTRAINT_H
#define OBSERVABLE_CONTINUOUS_CONSTRAINT_H

//	RapidFit Headers
#include "IConstraint.h"
//	System Headers
#include <string>

using namespace std;

class ObservableContinuousConstraint : public IConstraint
{
	public:
		ObservableContinuousConstraint();
		ObservableContinuousConstraint( string, double, double, string );
		~ObservableContinuousConstraint();

		void SetMinimum(double);
		void SetMaximum(double);
		void SetLimits(double, double);

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
		double minimum;
		double maximum;
		string unit;
};

#endif
