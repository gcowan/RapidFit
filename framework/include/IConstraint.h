/**
        @interface IConstraint

        Interface for all constraints defining a phase space boundary

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef I_CONSTRAINT_H
#define I_CONSTRAINT_H

#include "Observable.h"
#include "TRandom3.h"
#include <string>
#include <vector>

using namespace std;

class IConstraint
{
	public:
		//Return true if an observable value is compatible with the constraint
		virtual bool CheckObservable( Observable* ) = 0;

		//Return a pointer to an observable randomly generated from the constraint
		virtual Observable * CreateObservable() = 0;
		virtual Observable * CreateObservable( TRandom3* ) = 0;

		//Return the unit
		virtual string GetUnit() = 0;

		//If the observable is continuous, return the maximum and minimum
		virtual double GetMaximum() = 0;
		virtual double GetMinimum() = 0;
		//If the observable is discrete, return the possible values
		virtual vector<double> GetValues() = 0;

		//Return whether the observable is discrete
		virtual bool IsDiscrete() = 0;
};

#endif
