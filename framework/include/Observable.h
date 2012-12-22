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

		Observable( const string );

		Observable( const string, const double, const string );

		Observable( const Observable& input );

		/*!
		 * @brief Destruction
		 */
		~Observable();

		double GetValue() const;

		string GetUnit() const;

		/*!
		 * @brief An ugly lie, but it keeps the compiler happy for more abstract functions which should be treated as const!
		 */
		void SetBinNumber( const int ) const;
		void SetBkgBinNumber( const int ) const;
		int GetBinNumber() const;
		int GetBkgBinNumber() const;

		double GetAcceptance() const;
		double GetBkgAcceptance() const;
		void SetAcceptance( const double ) const;
		void SetBkgAcceptance( const double ) const;

		void SetOffSet( const double ) const;
		double GetOffSet() const;

		string GetName() const;

		void Print() const;

		void SetObservable( const Observable* );
		void ExternallySetValue( const double );

	private:
		string name;
		double value;
		string unit;
		mutable int bin_num;
		mutable int bkg_bin_num;
		mutable double acceptance;
		mutable double bkg_acceptance;
		mutable double offset;
};

#endif

