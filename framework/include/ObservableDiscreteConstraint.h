/**
 * @class ObservableDiscreteConstraint
 *
 * A constraint defining an observable that can take any of a list of discrete values.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef OBSERVABLE_DISCRETE_CONSTRAINT_H
#define OBSERVABLE_DISCRETE_CONSTRAINT_H

//	RapidFit Headers
#include "IConstraint.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class ObservableDiscreteConstraint : public IConstraint
{
	public:
		ObservableDiscreteConstraint( string, vector<double>, string, string="" );
		ObservableDiscreteConstraint( const IConstraint* input );
		~ObservableDiscreteConstraint();

		void AddValue(double);
		void SetValues( vector<double> );

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
		vector<double> allValues;
		string unit;
		string tf1;
};

#endif

