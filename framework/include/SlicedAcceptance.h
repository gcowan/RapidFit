/**
  @class SlicedAcceptance

  A class for holding a sliced propertime acceptance

  @author Pete Clarke
  @data 2011-06-07
  */

#ifndef SLICED_ACCEPTANCE_H
#define SLICED_ACCEPTANCE_H

//	System Headers
#include <iostream>
#include<fstream>
#include <cstdlib>
#include <string>
#include <vector>
using namespace std;

//=======================================

class AcceptanceSlice{

	public:
		AcceptanceSlice( double tl, double th, double h ) : _tlow(tl), _thigh(th), _height(h) {}
		double tlow()   const { return _tlow ; }
		double thigh()  const { return _thigh ; }
		double height() const { return _height ; }
	private:
		double _tlow;
		double _thigh;
		double _height;
};


//=======================================
class SlicedAcceptance
{
	public:

		//Constructors
		SlicedAcceptance( double tlow, double thigh ) ;
		SlicedAcceptance( double tlow, double thigh, double beta  ) ;
		SlicedAcceptance( string s  ) ;
		SlicedAcceptance( string s, string s  ) ;

		// Methods for numerator of PDF to return acceptance for event
		double getValue( double time ) const ;

		// Methods for the normalisation integral in slices
		int numberOfSlices() const ;
		AcceptanceSlice * getSlice( int slice ) const ;

	private:	
		//      Uncopyable!
		SlicedAcceptance ( const SlicedAcceptance& );
		SlicedAcceptance& operator = ( const SlicedAcceptance& );

		double stream(ifstream& stream) ;

		vector <AcceptanceSlice*> slices ;
		AcceptanceSlice* nullSlice ;	
		double tlow;
		double thigh;
		double beta;	
};

#endif

