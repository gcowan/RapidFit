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
		AcceptanceSlice( double tl, double th, double h ) : _tlow(tl), _thigh(th), _height(h) {}
		AcceptanceSlice( const AcceptanceSlice& input );
		double tlow()   const { return _tlow ; }
		double thigh()  const { return _thigh ; }
		double height() const { return _height ; }
	private:
		AcceptanceSlice& operator = ( const AcceptanceSlice& );
		double _tlow;
		double _thigh;
		double _height;
};


//=======================================
class SlicedAcceptance
{
	public:

		//Constructors
		~SlicedAcceptance();
		SlicedAcceptance( double tlow, double thigh );
		SlicedAcceptance( double tlow, double thigh, double beta  );
		SlicedAcceptance( string s  );
		SlicedAcceptance( string s1, string s2  );
		//	Copy Constructor
		SlicedAcceptance( const SlicedAcceptance& input );

		// Methods for numerator of PDF to return acceptance for event
		double getValue( const double time ) const;

		double getValue( const Observable* time, const double timeOffset=0. ) const;

		// Methods for the normalisation integral in slices
		unsigned int numberOfSlices() const;
		AcceptanceSlice * getSlice( const unsigned int slice ) const;

	private:	
		//      Uncopyable!
		SlicedAcceptance& operator = ( const SlicedAcceptance& );

		double stream(ifstream& stream);

		vector <AcceptanceSlice*> slices;
		AcceptanceSlice* nullSlice;
		double tlow;
		double thigh;
		double beta;

		mutable size_t uniqueID;
};

#endif

