#ifndef RAPIDFIT_UTILS_TEMPLATE_H
#define RAPIDFIT_UTILS_TEMPLATE_H

#include "TList.h"
#include "TClass.h"
#include "TObject.h"
#include "TString.h"

#include <vector>
#include <iostream>
#include <cmath>

using namespace::std;

class Template_Functions : public TObject
{
	public:

		//      Print a vector of variables which are compatible with cout with a single print command
		template<class T> static void print( const vector<T>& output, bool new_line=true );

		//	Print a vector of vectors to screen in some almost human reable format
		template<class T> static void print( const vector<vector<T> >& output );

		//      Print a vector of pairs of variables which are compatible with cout with a single command
		template<class T, class U> static void print( const vector<pair<T,U> >& output, bool new_line=true );

		//      Return a vector of the first objects in a vector of pairs
		template<class T, class U > static vector<T> return_first( const vector<pair<T,U> >& input );

		//      Return a vector of the first objects in a vector of pairs
		template<class T, class U > static vector<vector<T> > return_first( const vector<pair<vector<T>,U> >& input );

		//      return a vector of the second objects in a vector of pairs
		template<class T, class U > static vector<U> return_second( const vector<pair<T,U> >& input );

		//      return a vector of the second objects in a vector of pairs
		template<class T, class U > static vector<vector<U> > return_second( const vector<pair<T,vector<U> > >& input );

		//      empty a vector
		template<class T> static void clear( vector<T>* input );

		//      empty a vector of pairs
		template<class T, class U> static void clear( vector<pair<T,U> >* input );

		//	Convert a vector of pairs to a pair of vectors of independent type
		template<class T, class U > static pair<vector<T>,vector<U> > reparam ( const vector<pair<T,U> > input );

		//	Convert a pair of vectors to a vector of pairs of independent type
		template<class T, class U> static vector<pair<T,U> > reparam( const pair<vector<T>,vector<U> > input );

		//	Print to screen in a human readable(ish) format a pair of vectors of independent type
		template<class T, class U> static void print( pair<vector<T>,vector<U> > input );

		//	Make a vector which contains the flat result of a vector of vectors of the same type
		template<class T> static vector<T> concatonnate( const vector<vector<T> >& input );

		//	Convert a vector of type U to a vector of type T as long as each object can be directly recast
		template<class T, class U> static vector<T> convert_vector( const vector<U>& input );

		//	Convert a vector of vectors of U to a vector of vector of U (as ling as T can be ditectly recast to U)
		template<class T, class U> static vector<vector<T> > convert_nested_vector( const vector<vector<U> >& input );

		//	Unofotunatley I'm quite proud of this one for the fact the level of abstraction makes my head hurt but was easy to think up
		//	This 'rotates' a vector from being col x row to being row x col
		template<class T> static vector<vector<T> > rotate( const vector<vector<T> >& input );

		//	Return the maxima from a vector of type T (as long as T knows of '<' operator)
		template<class T> static T get_maximum( const vector<T>& input );

		//	Return the minima from a vector of type T (as long as T knows of '>' operator)
		template<class T> static T get_minimum( const vector<T>& input );

		//	Get the 'statistically optimal' number of bins for a vector of data if it is to be put in a histogram
		template<class T> static int get_optimal_histo_bins( const vector<T>& input );

		template<class T> static vector<T> reverse( const vector<T>& input );

		template<class T> static vector<T> request_element( const vector<vector<T> >& input, unsigned int index );

		static void PrintPrimatives( TList* thisList );

		static void Print();

		Template_Functions() {};
		virtual ~Template_Functions() {};
		ClassDef( Template_Functions, 1 );

};

#endif

