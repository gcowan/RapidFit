/*!
 * @class ObservableContinuousConstraint
 *
 * A constraint defining an observable that can take any value between a given maximum and minimum.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef OBSERVABLE_CONTINUOUS_CONSTRAINT_H
#define OBSERVABLE_CONTINUOUS_CONSTRAINT_H

//	RapidFit Headers
#include "IConstraint.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class ObservableContinuousConstraint : public IConstraint
{
	public:
		ObservableContinuousConstraint( string, double, double, string, string="" );
		ObservableContinuousConstraint( const IConstraint* input );
		~ObservableContinuousConstraint();

		void SetMinimum(double);
		void SetMaximum(double);
		void SetLimits(double, double);

		//Interface functions
		virtual bool CheckObservable( Observable* ) const;
		virtual Observable* CreateObservable() const;
		virtual Observable* CreateObservable( TRandom3* ) const;
		virtual string GetName() const;
		virtual string GetUnit() const;
		virtual double GetMaximum() const;
		virtual double GetMinimum() const;
		virtual vector<double> GetValues() const;
		virtual bool IsDiscrete() const;

		virtual void Print() const;

		virtual Observable* GetMidRangeValue() const;

		virtual string GetTF1() const;
		virtual void SetTF1( const string );

		string XML() const;
	private:
		string name;
		double minimum;
		double maximum;
		string unit;
		string tf1;
};

#endif

