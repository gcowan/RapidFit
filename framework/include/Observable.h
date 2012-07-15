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
		void SetBkgBinNumber( int ) const;
		int GetBinNumber() const;
		int GetBkgBinNumber() const;

		double GetAcceptance() const;
		double GetBkgAcceptance() const;
		void SetAcceptance( double ) const;
		void SetBkgAcceptance( double ) const;

		void SetOffSet( double ) const;
		double GetOffSet() const;

		string GetName() const;

		void Print() const;

		void SetObservable( Observable* );
		void ExternallySetValue( double );
		void ExternallySetError( double );

	private:
		string name;
		double value;
		double error;
		string unit;
		mutable int bin_num;
		mutable int bkg_bin_num;
		mutable double acceptance;
		mutable double bkg_acceptance;
		mutable double offset;
};

#endif

