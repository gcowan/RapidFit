
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

		//PseudoObservable& operator= ( const PseudoObservable& );

		~PseudoObservable();

		void AddFunction( double (*pseudoRelation)(vector<double>) );
		void AddFunction( double (*pseudoRelation)(vector<double>,vector<double>) );

		void AddErrorFunction( double (*pseudoRelation)(vector<double>) );
		void AddErrorFunction( double (*pseudoRelation)(vector<double>,vector<double>) );

		void AddDependency( ObservableRef Input );

		void AddDependencies( vector<ObservableRef> Input );

		/*!
		 * @brief This controls whether the PseudoObservable is to be Recalulated
		 * 
		 * @param Input      true: yes PseudoObservable will be recalculated
		 *                   false: no PseudoObservable will NOT be recalculated
		 */
		void SetValid( bool Input );

		bool GetValid();

		vector<ObservableRef>* GetDependencies();

		void SetInput( vector<double> );

		double GetValue() const;

		double GetError() const;

		string GetName() const;

		void SetUnit( string );

		string GetUnit() const;

		Observable* GetPseudoObservable();

		int GetIndex() const;

		void SetIndex( int );

		void Print() const;

		void ExtraInput() const;

	private:
		double (*internal_function)(vector<double>);
		double (*internal_function_extra)(vector<double>,vector<double>);
		double (*internal_error_function)(vector<double>);
		double (*internal_error_function_extra)(vector<double>,vector<double>);

		vector<ObservableRef> Dependencies;

		mutable double Value;
		mutable double Error;

		mutable bool valid;

		string Name;
		string Unit;

		int Index;

		vector<double> internal_Input;

		Observable* internalObservable;
};

#endif

