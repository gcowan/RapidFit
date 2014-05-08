
#include "TList.h"
#include "TObject.h"
#include "TClass.h"
#include "TString.h"

#include "Template_Functions.h"

#include <iostream>
#include <complex>

using namespace::std;


//      Print a vector of variables which are compatible with cout with a single print command
template<class T> void Template_Functions::print( const vector<T>& output, bool new_line )
{
	//cout << "T:" <<output[0]<< endl;
	typename std::vector<T>::const_iterator output_i;
	for( output_i = output.begin(); output_i != output.end(); ++output_i )
	{
		cout << *output_i ;
		if( new_line ) { cout << endl; }
		else { cout << ",\t" << endl; }
	}
}

template void Template_Functions::print( const vector<int>&, bool new_line );
template void Template_Functions::print( const vector<double>&, bool new_line );
template void Template_Functions::print( const vector<bool>&, bool new_line );
template void Template_Functions::print( const vector<float>&, bool new_line );
template void Template_Functions::print( const vector<string>&, bool new_line );
template void Template_Functions::print( const vector<TString>&, bool new_line );
template void Template_Functions::print( const vector<complex<double> >&, bool new_line );
template void Template_Functions::print( const vector<complex<int> >&, bool new_line );

//      Print a vector of vectors to screen in some almost human reable format
template<class T> void Template_Functions::print( const vector<vector<T> >& output )
{
	typename std::vector<vector<T> >::const_iterator output_i;
	typename std::vector<T>::const_iterator output_j;
	for( output_i = output.begin(); output_i != output.end(); ++output_i )
	{
		for( output_j = output_i->begin(); output_j != output_i->end(); ++ output_j )
		{
			cout << *output_j << ",\t";
		}
		cout << endl;
	}
}

template void Template_Functions::print( const vector<vector<int> >& );
template void Template_Functions::print( const vector<vector<double> >& );
template void Template_Functions::print( const vector<vector<bool> >& );
template void Template_Functions::print( const vector<vector<float> >& );
template void Template_Functions::print( const vector<vector<string> >& );
template void Template_Functions::print( const vector<vector<TString> >& );
template void Template_Functions::print( const vector<vector<complex<double> > >& );
template void Template_Functions::print( const vector<vector<complex<int> > >& );

//      Print a vector of pairs of variables which are compatible with cout with a single command
template<class T, class U> void Template_Functions::print( const vector<pair<T,U> >& output, bool new_line )
{
	typename std::vector<std::pair<T,U> >::const_iterator output_i;
	for( output_i = output.begin(); output_i != output.end(); ++output_i )
	{
		cout << output_i->first << "  ,  " << output_i->second ;
		if( new_line ) { cout << endl; }
		else { cout << " ::  " << endl; }
	}
}

template void Template_Functions::print( const vector< pair<int,int> >&, bool );
template void Template_Functions::print( const vector< pair<int,double> >&, bool );
template void Template_Functions::print( const vector< pair<int,bool> >&, bool );
template void Template_Functions::print( const vector< pair<int,string> >&, bool );
template void Template_Functions::print( const vector< pair<int,TString> >&, bool );
template void Template_Functions::print( const vector< pair<int,float> >&, bool );
template void Template_Functions::print( const vector< pair<string,int> >&, bool );
template void Template_Functions::print( const vector< pair<string,double> >&, bool );
template void Template_Functions::print( const vector< pair<string,bool> >&, bool );
template void Template_Functions::print( const vector< pair<string,float> >&, bool );
template void Template_Functions::print( const vector< pair<string,string> >&, bool );
template void Template_Functions::print( const vector< pair<string,TString> >&, bool );
template void Template_Functions::print( const vector< pair<double,double> >&, bool );
template void Template_Functions::print( const vector< pair<float,float> >&, bool );
template void Template_Functions::print( const vector< pair<bool,bool> >&, bool );

//      Return a vector of the first objects in a vector of pairs
template<class T, class U > vector<T> Template_Functions::return_first( const vector<pair<T,U> >& input )
{
	typename std::vector<std::pair<T,U> >::const_iterator input_i;
	vector<T> output;
	for( input_i = input.begin(); input_i != input.end(); ++input_i )
	{
		output.push_back( input_i->first );
	}
	return output;
}

template vector<double> Template_Functions::return_first( const vector<pair<double,double> >& );
template vector<double> Template_Functions::return_first( const vector<pair<double,TString> >& );
template vector<int> Template_Functions::return_first( const vector<pair<int,int> >& );
template vector<int> Template_Functions::return_first( const vector<pair<int,double> >& );
template vector<int> Template_Functions::return_first( const vector<pair<int,string> >& );
template vector<float> Template_Functions::return_first( const vector<pair<float,float> >& );
template vector<string> Template_Functions::return_first( const vector<pair<string,int> >& );
template vector<string> Template_Functions::return_first( const vector<pair<string,double> >& );

//      Return a vector of the first objects in a vector of pairs
template<class T, class U > vector<vector<T> > Template_Functions::return_first( const vector<pair<vector<T>,U> >& input )
{
	typename std::vector<std::pair<std::vector<T>,U> >::const_iterator input_i;
	vector<vector<T> > output;
	for( input_i = input.begin(); input_i != input.end(); ++input_i )
	{
		output.push_back( input_i->first );
	}
	return output;
}

template vector<vector<double> > Template_Functions::return_first( const vector<pair<vector<double>,double> >& );
template vector<vector<double> > Template_Functions::return_first( const vector<pair<vector<double>,TString> >& );
template vector<vector<int> > Template_Functions::return_first( const vector<pair<vector<int>,int> >& );
template vector<vector<int> > Template_Functions::return_first( const vector<pair<vector<int>,double> >& );
template vector<vector<int> > Template_Functions::return_first( const vector<pair<vector<int>,string> >& );
template vector<vector<float> > Template_Functions::return_first( const vector<pair<vector<float>,float> >& );
template vector<vector<string> > Template_Functions::return_first( const vector<pair<vector<string>,int> >& );
template vector<vector<string> > Template_Functions::return_first( const vector<pair<vector<string>,double> >& );

//      return a vector of the second objects in a vector of pairs
template<class T, class U > vector<U> Template_Functions::return_second( const vector<pair<T,U> >& input )
{
	typename std::vector<std::pair<T,U> >::const_iterator input_i;
	vector<U> output;
	for( input_i = input.begin(); input_i != input.end(); ++input_i )
	{
		output.push_back( input_i->second );
	}
	return output;
}

template vector<double> Template_Functions::return_second( const vector<pair<double,double> >& );
template vector<int> Template_Functions::return_second( const vector<pair<int,int> >& );
template vector<double> Template_Functions::return_second( const vector<pair<int,double> >& );
template vector<string> Template_Functions::return_second( const vector<pair<int,string> >& );
template vector<float> Template_Functions::return_second( const vector<pair<float,float> >& );
template vector<int> Template_Functions::return_second( const vector<pair<string,int> >& );
template vector<double> Template_Functions::return_second( const vector<pair<string,double> >& );
template vector<string> Template_Functions::return_second( const vector<pair<string,string> >& );
template vector<string> Template_Functions::return_second( const vector<pair<double,string> >& );
template vector<TString> Template_Functions::return_second( const vector<pair<double,TString> >& );

//      return a vector of the second objects in a vector of pairs
template<class T, class U > vector<vector<U> > Template_Functions::return_second( const vector<pair<T,vector<U> > >& input )
{
	typename std::vector<std::pair<T,std::vector<U> > >::const_iterator input_i;
	vector<vector<U> > output;
	for( input_i = input.begin(); input_i != input.end(); ++input_i )
	{
		output.push_back( input_i->second );
	}
	return output;
}

template vector<vector<double> > Template_Functions::return_second( const vector<pair<double,vector<double> > >& );
template vector<vector<int> > Template_Functions::return_second( const vector<pair<int,vector<int> > >& );
template vector<vector<int> > Template_Functions::return_second( const vector<pair<double,vector<int> > >& );
template vector<vector<int> > Template_Functions::return_second( const vector<pair<string,vector<int> > >& );
template vector<vector<float> > Template_Functions::return_second( const vector<pair<float,vector<float> > >& );
template vector<vector<string> > Template_Functions::return_second( const vector<pair<int,vector<string> > >& );
template vector<vector<string> > Template_Functions::return_second( const vector<pair<double,vector<string> > >& );

//      empty a vector
template<class T> void Template_Functions::clear( vector<T>* input )
{
	while( !input->empty() ) { input->pop_back(); }
}

template void Template_Functions::clear( vector<double>* );
template void Template_Functions::clear( vector<int>* );
template void Template_Functions::clear( vector<bool>* );
template void Template_Functions::clear( vector<string>* );
template void Template_Functions::clear( vector<TString>* );
template void Template_Functions::clear( vector<complex<int> >* );
template void Template_Functions::clear( vector<complex<double> >* );

//      empty a vector of pairs
template<class T, class U> void Template_Functions::clear( vector<pair<T,U> >* input )
{
	while( !input->empty() ) { input->pop_back(); }
}

template void Template_Functions::clear( vector<pair<double,double> >* );
template void Template_Functions::clear( vector<pair<int,int> >* );
template void Template_Functions::clear( vector<pair<bool,bool> >* );
template void Template_Functions::clear( vector<pair<string,string> >* );
template void Template_Functions::clear( vector<pair<TString,TString> >* );

//      Convert a vector of pairs to a pair of vectors of independent type
template<class T, class U > pair<vector<T>,vector<U> > Template_Functions::reparam ( const vector<pair<T,U> > input )
{
	typename std::vector<std::pair<T,U> >::const_iterator input_i;
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

template pair<vector<int>,vector<int> > Template_Functions::reparam ( const vector<pair<int,int> > );
template pair<vector<int>,vector<double> > Template_Functions::reparam ( const vector<pair<int,double> > );
template pair<vector<int>,vector<float> > Template_Functions::reparam ( const vector<pair<int,float> > );
template pair<vector<double>,vector<double> > Template_Functions::reparam ( const vector<pair<double,double> > );
template pair<vector<float>,vector<float> > Template_Functions::reparam ( const vector<pair<float,float> > );

//      Convert a pair of vectors to a vector of pairs of independent type
template<class T, class U> vector<pair<T,U> > Template_Functions::reparam( const pair<vector<T>,vector<U> > input )
{
	typename std::vector<T>::const_iterator input_i = input.first.begin();
	typename std::vector<U>::const_iterator input_j = input.second.begin();
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

template vector<pair<int,int> > Template_Functions::reparam( const pair<vector<int>,vector<int> > input );
template vector<pair<int,double> > Template_Functions::reparam( const pair<vector<int>,vector<double> > input );
template vector<pair<int,float> > Template_Functions::reparam( const pair<vector<int>,vector<float> > input );
template vector<pair<double,double> > Template_Functions::reparam( const pair<vector<double>,vector<double> > input );
template vector<pair<float,float> > Template_Functions::reparam( const pair<vector<float>,vector<float> > input );

//      Print to screen in a human readable(ish) format a pair of vectors of independent type
template<class T, class U> void Template_Functions::print( pair<vector<T>,vector<U> > input )
{
	print( reparam( input ) );
}

template void Template_Functions::print( pair<vector<int>,vector<int> > input );
template void Template_Functions::print( pair<vector<double>,vector<double> > input );
template void Template_Functions::print( pair<vector<float>,vector<float> > input );
template void Template_Functions::print( pair<vector<string>,vector<string> > input );
template void Template_Functions::print( pair<vector<TString>,vector<TString> > input );

//      Make a vector which contains the flat result of a vector of vectors of the same type
template<class T> vector<T> Template_Functions::concatonnate( const vector<vector<T> >& input )
{
	typename std::vector<vector<T> >::const_iterator input_i=input.begin();
	vector<T> output;
	for( ; input_i != input.end(); ++input_i )
	{
		typename std::vector<T>::const_iterator input_j=input_i->begin();
		for( ; input_j != input_i->end(); ++input_j )
		{
			output.push_back( *input_j );
		}
	}
	return output;
}

template vector<int> Template_Functions::concatonnate( const vector<vector<int> >& input );
template vector<bool> Template_Functions::concatonnate( const vector<vector<bool> >& input );
template vector<float> Template_Functions::concatonnate( const vector<vector<float> >& input );
template vector<double> Template_Functions::concatonnate( const vector<vector<double> >& input );
template vector<string> Template_Functions::concatonnate( const vector<vector<string> >& input );
template vector<TString> Template_Functions::concatonnate( const vector<vector<TString> >& input );
template vector<complex<int> > Template_Functions::concatonnate( const vector<vector<complex<int> > >& input );
template vector<complex<double> > Template_Functions::concatonnate( const vector<vector<complex<double> > >& input );

//      Convert a vector of type U to a vector of type T as long as each object can be directly recast
template<class T, class U> vector<T> Template_Functions::convert_vector( const vector<U>& input )
{
	typename std::vector<U>::const_iterator input_i = input.begin();
	vector<T> output;
	for( ; input_i != input.end(); ++input_i )
	{
		output.push_back( T(*input_i) );
	}
	return output;
}

template vector<int> Template_Functions::convert_vector( const vector<double>& input );
template vector<int> Template_Functions::convert_vector( const vector<float>& input );
template vector<float> Template_Functions::convert_vector( const vector<int>& input );
template vector<double> Template_Functions::convert_vector( const vector<int>& input );

//      Convert a vector of vectors of U to a vector of vector of U (as ling as T can be ditectly recast to U)
template<class T, class U> vector<vector<T> > Template_Functions::convert_nested_vector( const vector<vector<U> >& input )
{
	typename std::vector<vector<U> >::const_iterator input_i = input.begin();
	vector<vector<T> > output;
	for( ; input_i != input.end(); ++input_i )
	{
		typename std::vector<U>::const_iterator input_j = input_i->begin();
		vector<T> output_temp;
		for( ; input_j != input_i->end(); ++input_j )
		{
			output_temp.push_back( T(*input_j) );
		}
		output.push_back( output_temp );
	}
	return output;
}

template vector<vector<int> > Template_Functions::convert_nested_vector( const vector<vector<float> >& input );
template vector<vector<int> > Template_Functions::convert_nested_vector( const vector<vector<double> >& input );
template vector<vector<float> > Template_Functions::convert_nested_vector( const vector<vector<int> >& input );
template vector<vector<double> > Template_Functions::convert_nested_vector( const vector<vector<int> >& input );

//      Unofotunatley I'm quite proud of this one for the fact the level of abstraction makes my head hurt but was easy to think up
//      This 'rotates' a vector from being col x row to being row x col
template<class T> vector<vector<T> > Template_Functions::rotate( const vector<vector<T> >& input )
{
	if( input.empty() ) return input;
	vector<vector<T> > output;
	typename vector<vector<T> >::const_iterator input_i = input.begin();
	typename vector<vector<T> >::const_iterator input_i_begin = input_i;
	typename vector<vector<T> >::const_iterator input_i_end = input.end();
	typename vector<vector<T> >::iterator output_i;
	typename vector<vector<T> >::iterator output_i_begin;
	typename vector<vector<T> >::iterator output_i_end;
	typename vector<T>::const_iterator input_ij;
	typename vector<T>::const_iterator input_ij_end;
	typename vector<T>::const_iterator input_ij_begin;
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

template vector<vector<int> > Template_Functions::rotate( const vector<vector<int> >& input );
template vector<vector<double> > Template_Functions::rotate( const vector<vector<double> >& input );
template vector<vector<float> > Template_Functions::rotate( const vector<vector<float> >& input );
template vector<vector<bool> > Template_Functions::rotate( const vector<vector<bool> >& input );
template vector<vector<string> > Template_Functions::rotate( const vector<vector<string> >& input );
template vector<vector<TString> > Template_Functions::rotate( const vector<vector<TString> >& input );

//      Return the maxima from a vector of type T (as long as T knows of '<' operator)
template<class T> T Template_Functions::get_maximum( const vector<T>& input )
{
	typename vector<T>::const_iterator index_i = input.begin();
	typename vector<T>::const_iterator index_e = input.end();
	T output = *index_i;

	for( ; index_i != index_e; ++index_i )
	{
		if( output < *index_i ) output = *index_i;
	}
	return output;
}

template int Template_Functions::get_maximum( const vector<int>& input );
template float Template_Functions::get_maximum( const vector<float>& input );
template double Template_Functions::get_maximum( const vector<double>& input );

//      Return the minima from a vector of type T (as long as T knows of '>' operator)
template<class T> T Template_Functions::get_minimum( const vector<T>& input )
{
	typename vector<T>::const_iterator index_i = input.begin();
	typename vector<T>::const_iterator index_e = input.end();
	T output = *index_i;

	for( ; index_i != index_e; ++index_i )
	{
		if( output > *index_i ) output = *index_i;
	}
	return output;
}

template int Template_Functions::get_minimum( const vector<int>& input );
template double Template_Functions::get_minimum( const vector<double>& input );
template float Template_Functions::get_minimum( const vector<float>& input );

//      Get the 'statistically optimal' number of bins for a vector of data if it is to be put in a histogram
template<class T> int Template_Functions::get_optimal_histo_bins( const vector<T>& input )
{
	double mean = 0.;
	typename vector<T>::const_iterator index_i = input.begin();
	typename vector<T>::const_iterator index_e = input.end();
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

template int Template_Functions::get_optimal_histo_bins( const vector<int>& input );
template int Template_Functions::get_optimal_histo_bins( const vector<double>& input );
template int Template_Functions::get_optimal_histo_bins( const vector<float>& input );

template<class T> vector<T> Template_Functions::reverse( const vector<T>& input )
{
	typename vector<T>::const_iterator index_b = input.begin()-1;
	typename vector<T>::const_iterator index_i = input.end();
	vector<T> output;
	for( --index_i; index_i != index_b; --index_i )
	{
		output.push_back( *index_i );
	}
	return output;
}

template vector<int> Template_Functions::reverse( const vector<int>& input );
template vector<bool> Template_Functions::reverse( const vector<bool>& input );
template vector<double> Template_Functions::reverse( const vector<double>& input );
template vector<float> Template_Functions::reverse( const vector<float>& input );
template vector<string> Template_Functions::reverse( const vector<string>& input );
template vector<TString> Template_Functions::reverse( const vector<TString>& input );
template vector<complex<int> > Template_Functions::reverse( const vector<complex<int> >& input );
template vector<complex<double> > Template_Functions::reverse( const vector<complex<double> >& input );

template<class T> vector<T> request_element( const vector<vector<T> >& input, unsigned int index )
{
	typename vector<vector<T> >::const_iterator index_i = input.begin();
	vector<T> output;
	for( ; index_i != input.end(); ++index_i )
	{
		output.push_back( (*index_i)[index] );
	}
	return output;
}

template vector<int> request_element( const vector<vector<int> >& input, unsigned int index );
template vector<double> request_element( const vector<vector<double> >& input, unsigned int index );
template vector<float> request_element( const vector<vector<float> >& input, unsigned int index );
template vector<bool> request_element( const vector<vector<bool> >& input, unsigned int index );
template vector<string> request_element( const vector<vector<string> >& input, unsigned int index );
template vector<TString> request_element( const vector<vector<TString> >& input, unsigned int index );
template vector<complex<double> > request_element( const vector<vector<complex<double> > >& input, unsigned int index );
template vector<complex<int> > request_element( const vector<vector<complex<int> > >& input, unsigned int index );

void Template_Functions::PrintPrimatives( TList* thisList )
{
	//TList* thisList = thisObject->GetListOfPrimitives();

	//cout << "TObject: " << thisObject->GetName() << "\t of type: " << thisObject->GetTypeName() << "\tat: " << thisObject << endl;
	cout << "I have knowledge of: " << endl;

	for( unsigned int i=0; i< (unsigned) thisList->GetSize(); ++i )
	{
		TObject* objPointer = thisList->At( i );
		TString Name = objPointer->GetName();
		TString Type = objPointer->IsA()->GetImplFileName();
		cout << "Name: " << Name << "\tType: " << Type << "\tAt: " << objPointer << endl;
	}

	cout << endl;
}

void Template_Functions::Print()
{
	cout << "Hello From Template_Functions!" << endl;
}

