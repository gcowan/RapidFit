/**
  @class SlicedAcceptance

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */

#pragma once
#ifndef SLICED_ACCEPTANCE_H
#define SLICED_ACCEPTANCE_H

//	RapidFir Headers
#include "Observable.h"
//	System Headers
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

using namespace::std;

//=======================================

class AcceptanceSlice
{
	public:
		/*!
		 * @brief Constructor
		 *
		 * @param tl
		 *
		 * @param th
		 *
		 * @param h
		 *
		 */
		AcceptanceSlice( double tl, double th, double h ) : _tlow(tl), _thigh(th), _height(h) {}

		/*!
		 * @brief Copy Constructor
		 */
		AcceptanceSlice( const AcceptanceSlice& input );

		/*!
		 * @brief
		 */
		double tlow()   const { return _tlow ; }

		/*!
		 * @brief
		 */
		double thigh()  const { return _thigh ; }

		/*!
		 * @brief
		 */
		 double height() const { return _height ; }

		void Print() const { cout << "Min: " << _tlow << endl << "Max: " << _thigh << endl << "Height: " << _height << endl; }

	private:
		/*!
		 * @brief Don't Copy This Way
		 */
		AcceptanceSlice& operator = ( const AcceptanceSlice& );

		/*!
		 * Internal Objects to an Acceptance Slice
		 */
		double _tlow;
		double _thigh;
		double _height;
};


//=======================================
class SlicedAcceptance
{
	public:

		/*!
		 * @brief Destructor
		 */
		~SlicedAcceptance();

		/*!
		 * @brief
		 *
		 * @param tlow
		 *
		 * @param thigh
		 *
		 */
		SlicedAcceptance( double tlow, double thigh, bool quiet=false );

		/*!
		 * @brief 
		 *
		 * @param tlow
		 *
		 * @param thigh
		 *
		 */
		SlicedAcceptance( double tlow, double thigh, double beta, bool quiet=false  );

		/*!
		 * @brief
		 *
		 * @param s
		 *
		 */
		SlicedAcceptance( string s, bool quiet=false );

		/*!
		 * @brief
		 *
		 * @param s1
		 *
		 * @param s2
		 *
		 */
		SlicedAcceptance( string s1, string s2, bool quiet=false  );

		/*!
		 * @brief Copy Constructor
		 */
		SlicedAcceptance( const SlicedAcceptance& input );

		/*!
		 * @brief Method for numerator of PDF to return acceptance for event
		 *
		 * @param time
		 *
		 * @return
		 */
		double getValue( const double time ) const;

		/*!
		 * @brief Method for numerator of PDF to return acceptance for event
		 *
		 * @param time
		 *
		 * @param timeOffset
		 *
		 * @return
		 */
		double getValue( const Observable* time, const double timeOffset=0. ) const;

		/*!
		 * @brief Method for the normalisation integral in slices
		 *
		 * @return
		 */
		unsigned int numberOfSlices() const;

		/*!
		 * @brief Method for the normalisation integral in slices
		 *
		 * @param slice
		 *
		 * @return
		 */
		AcceptanceSlice* getSlice( const unsigned int slice ) const;

		unsigned int findSliceNum( const Observable* time, const double offSet=0. ) const;

		bool GetIsSorted() const;

		double GetMax() const;
		double GetMin() const;

		void Print() const;
	private:
		//      Uncopyable!
		SlicedAcceptance& operator = ( const SlicedAcceptance& );

		double stream(ifstream& stream);

		vector <AcceptanceSlice*> slices;
		AcceptanceSlice* nullSlice;
		double tlow;
		double thigh;
		double beta;
		bool _sortedSlices;

		bool isSorted() const;

		void FindMaxMin() const;

		mutable double t_max;
		mutable double t_min;
		mutable bool maxminset;

		mutable bool _hasChecked;
		mutable bool _storedDecision;
};

#endif

