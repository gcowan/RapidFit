
#pragma once
#ifndef _PseudoObservable_H
#define _PseudoObservable_H

///	RapidFit Headers
#include "ObservableRef.h"
#include "Observable.h"
///	System Headers
#include <vector>
#include <string>

using namespace::std;

class PseudoObservable
{
	public:
		PseudoObservable();

		PseudoObservable( string Name );

		PseudoObservable( const PseudoObservable& );

		PseudoObservable& operator= ( const PseudoObservable& );

		~PseudoObservable();

		void AddFunction( double (*pseudoRelation)(vector<double>) );

		void AddDependency( ObservableRef Input );

		void AddDependencies( vector<ObservableRef> Input );

		/*!
		 * @brief This controls whether the PseudoObservable is to be Recalulated
		 * 
		 * @param Input      true: yes PseudoObservable will be recalculated
		 *                   false: no PseudoObservable will NOT be recalculated
		 */
		void SetValid( bool Input ) const;

		bool GetValid() const;

		bool GetValid( const vector<double> input ) const;

		vector<ObservableRef>* GetDependencies() const;

		void SetInput( const vector<double> );

		double GetValue() const;

		string GetName() const;

		double GetPseudoObservable() const;

		int GetIndex() const;

		void SetIndex( const int ) const;

		void Print() const;

		vector<double> GetCacheInput() const;

		double Eval( const vector<double> ) const;

	private:
		double (*internal_function)(vector<double>);

		mutable vector<ObservableRef> Dependencies;

		mutable double Value;
		mutable bool valid;

		string Name;
		mutable int Index;

		vector<double> internal_Input;
};

#endif

