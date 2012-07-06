/**
        @class Observable

        The measured value of a particular variable for an event

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef OBSERVABLE_H
#define OBSERVABLE_H

//	System Headers
#include <string>

using namespace::std;

class Observable
{
	public:

		Observable( string );

		Observable( string, double, double, string );

		/*!
		 * @brief Destruction
		 */
		~Observable();

		double GetValue() const;

		double GetError() const;

		string GetUnit() const;

		/*!
		 * @brief An ugly lie, but it keeps the compiler happy for more abstract functions which should be treated as const!
		 */
		void SetBinNumber( int ) const;
		int GetBinNumber() const;

		double GetAcceptance() const;
		void SetAcceptance( double ) const;

		string GetName() const;

		void Print() const;

		void SetObservable( Observable* );

	private:
		string name;
		double value;
		double error;
		string unit;
		mutable int bin_num;
		mutable double acceptance;
};

#endif

