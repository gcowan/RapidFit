#ifndef RAPIDFIT_UTILS_TEMPLATE_H
#define RAPIDFIT_UTILS_TEMPLATE_H

#include <vector>
#include <iostream>
#include <cmath>
using namespace::std;

//      Print a vector of variables which are compatible with cout with a single print command
template<class T> void print ( vector<T> output, bool new_line=true )
{
	//cout << "T:" <<output[0]<< endl;
	typename std::vector<T>::iterator output_i;
	for( output_i = output.begin(); output_i != output.end(); ++output_i )
	{
		cout << *output_i ;
		if( new_line ) { cout << endl; }
		else { cout << ",\t" << endl; }
	}
}

//	Print a vector of vectors to screen in some almost human reable format
template<class T> void print ( vector<vector<T> > output )
{
	typename std::vector<vector<T> >::iterator output_i;
	typename std::vector<T>::iterator output_j;
	for( output_i = output.begin(); output_i != output.end(); ++output_i )
	{
		for( output_j = output_i->begin(); output_j != output_i->end(); ++ output_j )
		{
			cout << *output_j << ",\t";
		}
		cout << endl;
	}
}

//      Print a vector of pairs of variables which are compatible with cout with a single command
template<class T, class U> void print ( vector<pair<T,U> > output, bool new_line=true )
{
	typename std::vector<std::pair<T,U> >::iterator output_i;
	for( output_i = output.begin(); output_i != output.end(); ++output_i )
	{
		cout << output_i->first << "  ,  " << output_i->second ;
		if( new_line ) { cout << endl; }
		else { cout << " ::  " << endl; }
	}
}

//      Return a vector of the first objects in a vector of pairs
template<class T, class U > vector<T> return_first( vector<pair<T,U> > input )
{
	typename std::vector<std::pair<T,U> >::iterator input_i;
	vector<T> output;
	for( input_i = input.begin(); input_i != input.end(); ++input_i )
	{
		output.push_back( input_i->first );
	}
	return output;
}

//      return a vector of the second objects in a vector of pairs
template<class T, class U > vector<U> return_second( vector<pair<T,U> > input )
{
	typename std::vector<std::pair<T,U> >::iterator input_i;
	vector<U> output;
	for( input_i = input.begin(); input_i != input.end(); ++input_i )
	{
		output.push_back( input_i->second );
	}
	return output;
}

//      empty a vector
template<class T> void clear ( vector<T>* input )
{
	while( !input->empty() ) { input->pop_back(); }
}

//      empty a vector of pairs
template<class T, class U> void clear ( vector<pair<T,U> >* input )
{
	while( !input->empty() ) { input->pop_back(); }
}

//	Convert a vector of pairs to a pair of vectors of independent type
template<class T, class U > pair<vector<T>,vector<U> > reparam ( vector<pair<T,U> > input )
{
	typename std::vector<std::pair<T,U> >::iterator input_i;
	vector<T> output_first; vector<U> output_second;
	for( input_i=input.begin(); input_i!=input.end(); ++input_i )
	{
		output_first.push_back( input_i->first );
		output_second.push_back( input_i->second );
	}
	pair<vector<T>,vector<U> > output;
	output.first = output_first;
	output.second = output_second;
	return output;
}

//	Convert a pair of vectors to a vector of pairs of independent type
template<class T, class U> vector<pair<T,U> > reparam ( pair<vector<T>,vector<U> > input )
{
	typename std::vector<T>::iterator input_i = input.first.begin();
	typename std::vector<U>::iterator input_j = input.second.begin();
	vector<pair<T,U> > output;
	for(; input_i != input.first.end(); ++input_i, ++input_j )
	{
		pair<T,U> temp;
		temp.first = *input_i;
		temp.second = *input_j;
		output.push_back( temp );
	}
	return output;
}

//	Print to screen in a human readable(ish) format a pair of vectors of independent type
template<class T, class U> void print( pair<vector<T>,vector<U> > input )
{
	print( reparam( input ) );
}

//	Make a vector which contains the flat result of a vector of vectors of the same type
template<class T> vector<T> concatonnate( vector<vector<T> > input )
{
	typename std::vector<vector<T> >::iterator input_i=input.begin();
	vector<T> output;
	for( ; input_i != input.end(); ++input_i )
	{
		typename std::vector<T>::iterator input_j=input_i->begin();
		for( ; input_j != input_i->end(); ++input_j )
		{
			output.push_back( *input_j );
		}
	}
	return output;
}

//	Convert a vector of type U to a vector of type T as long as each object can be directly recast
template<class T, class U> vector<T> convert_vector( vector<U> input )
{
	typename std::vector<U>::iterator input_i = input.begin();
	vector<T> output;
	for( ; input_i != input.end(); ++input_i )
	{
		output.push_back( *input_i );
	}
	return output;
}

//	Convert a vector of vectors of U to a vector of vector of U (as ling as T can be ditectly recast to U)
template<class T, class U> vector<vector<T> > convert_nested_vector( vector<vector<U> > input )
{
	typename std::vector<vector<U> >::iterator input_i = input.begin();
	vector<vector<T> > output;
	for( ; input_i != input.end(); ++input_i )
	{
		typename std::vector<U>::iterator input_j = input_i->begin();
		vector<T> output_temp;
		for( ; input_j != input_i->end(); ++input_j )
		{
			output_temp.push_back( *input_j );
		}
		output.push_back( output_temp );
	}
	return output;
}


//	Unofotunatley I'm quite proud of this one for the fact the level of abstraction makes my head hurt but was easy to think up
//	This 'rotates' a vector from being col x row to being row x col
template<class T> vector<vector<T> > rotate( vector<vector<T> >& input )
{
	if( input.empty() ) return input;
	vector<vector<T> > output;
	typename vector<vector<T> >::iterator input_i = input.begin();
	typename vector<vector<T> >::iterator input_i_begin = input_i;
	typename vector<vector<T> >::iterator input_i_end = input.end();
	typename vector<vector<T> >::iterator output_i;
	typename vector<vector<T> >::iterator output_i_begin;
	typename vector<vector<T> >::iterator output_i_end;
	typename vector<T>::iterator input_ij;
	typename vector<T>::iterator input_ij_end;
	typename vector<T>::iterator input_ij_begin;
	typename vector<T>::iterator output_ij;
	ptrdiff_t output_ij_step;
	size_t i = input[0].size();
	size_t j = input.size();

	output.resize( i, vector<T>() );
	output_i_begin = output.begin();
	output_i_end = output.end();
	for( output_i = output_i_begin; output_i != output_i_end; ++output_i )
	{
		output_i->resize( j );
	}

	input_i = input.begin();

	for( ; input_i != input_i_end; ++input_i )
	{
		input_ij_begin = input_i->begin();
		input_ij_end = input_i->end();
		output_ij_step = (input_i-input_i_begin);
		input_ij = input_ij_begin;
		output_i = output_i_begin;

		for( ; input_ij != input_ij_end; ++input_ij, ++output_i )
		{
			output_ij = output_i->begin() + output_ij_step;
			*output_ij = *input_ij;
		}
	}

	return output;
}

//	Return the minima from a vector of type T (as long as T knows of '<' operator)
template<class T> T get_minimum( vector<T> input )
{
	typename vector<T>::iterator index_i = input.begin();
	typename vector<T>::iterator index_e = input.end();
	T output = *index_i;

	for( ; index_i != index_e; ++index_i )
	{
		if( output < *index_i ) output = *index_i;
	}
	return output;
}

//	Return the maxima from a vector of type T (as long as T knows of '>' operator)
template<class T> T get_maximum( vector<T> input )
{
	typename vector<T>::iterator index_i = input.begin();
	typename vector<T>::iterator index_e = input.end();
	T output = *index_i;

	for( ; index_i != index_e; ++index_i )
	{
		if( output > *index_i ) output = *index_i;
	}
	return output;
}

//	Get the 'statistically optimal' number of bins for a vector of data if it is to be put in a histogram
template<class T> int get_optimal_histo_bins( vector<T> input )
{
	double mean = 0.;
	typename vector<T>::iterator index_i = input.begin();
	typename vector<T>::iterator index_e = input.end();
	for( ; index_i != index_e; ++index_i )
	{
		mean+=*index_i;
	}
	mean/=(double)input.size();
	double variance=0.;
	for( index_i = input.begin(); index_i != index_e; ++index_i )
	{
		variance+=(*index_i)*(*index_i) - mean;
	}
	double width = 3.49 * sqrt( variance ) * pow( (double)input.size(), -(1.0/3.0) );

	double range = get_maximum( input ) - get_minimum( input );
	int wanted_bins = ceil( range / width );

	return wanted_bins;
}

#endif

