
#include "TF1.h"
#include "TPolyMarker3D.h"
#include "TString.h"

#include "Mathematics.h"

#include <vector>
#include <cmath>
#include <fstream>

using namespace::std;

bool Mathematics::double_equals_test(Double_t first, Double_t second)
{
	return fabs(first-second) < 0.0001;
}

//      This will find all of the unique coordinates stored in a 1D vector of doubles
vector<vector<Double_t> > Mathematics::Unique_Coords( vector<Double_t> input )
{
	vector<Double_t> output_temp = input;
	sort( output_temp.begin(), output_temp.end() );
	vector<Double_t>::iterator len;
	len = unique( output_temp.begin(), output_temp.end(), double_equals_test );
	output_temp.resize( len - output_temp.begin() );
	vector<vector<Double_t> > output( 1, output_temp );
	return output;
}
//	This will find all of the unique coordinates stored in a 2D vector of pairs of doubles
vector<vector<Double_t> > Mathematics::Unique_Coords( vector<pair<Double_t,Double_t> > input )
{
	//      Tolerance of the DOUBLE
	double DT = 1E-5;

	//      Get the number of coordinates in the TPolyMarker3D
	int number_of_points = (int)input.size();

	//      Will populate and return to the user
	vector<vector<Double_t> > Returnable_Coord_Data;

	//      Used within a for loop in logic, externally created/destroyed
	bool add_point = true;
	bool temp_decision_1 = true;
	bool temp_decision_2 = true;

	//      Run over all of the points contained within the TPolyMarker3D object
	for( int i=0; i< number_of_points; ++i )
	{
		//      Assume we haven't seen this point yet
		add_point = true;
		//      Check the anzats
		for( unsigned int j=0; j< Returnable_Coord_Data.size(); ++j )
		{
			temp_decision_1 = true;
			temp_decision_2 = true;
			if( fabs( input[i].first - Returnable_Coord_Data[j][0] ) < DT ) temp_decision_1 = false;
			if( fabs( input[i].second - Returnable_Coord_Data[j][1] ) < DT ) temp_decision_2 = false;

			//      If all of the information is the same do NOT add the point to the new vector of data
			if( ( temp_decision_1 == temp_decision_2 ) && ( temp_decision_1 == false ) )
				add_point = false;
		}
		//      If we haven't seen this point yet add it to the array of points
		if( add_point )
		{
			//cout << Coord_Data_pointer[i].first << "\t" << Coord_Data_pointer[i].second << endl;
			vector<Double_t> temp_vector;
			temp_vector.push_back( input[i].first );
			temp_vector.push_back( input[i].second );
			Returnable_Coord_Data.push_back( temp_vector );
		}
	}

	//      Return a 2D vector of unique corrdinates
	return Returnable_Coord_Data;
}

//      This will find all of the unique coordinates stored in a n-D vector of doubles
vector<vector<Double_t> > Mathematics::Unique_Coords( vector<vector<Double_t> > input )
{
	//      Tolerance of the DOUBLE
	double DT = 1E-5;

	//      Get the number of coordinates in the TPolyMarker3D
	int number_of_points = (int)input.size();

	//      Will populate and return to the user
	vector<vector<Double_t> > Returnable_Coord_Data;

	//      Used within a for loop in logic, externally created/destroyed
	bool add_point = true;

	//      Run over all of the points contained within the TPolyMarker3D object
	for( int i=0; i< number_of_points; ++i )
	{
		//      Assume we haven't seen this point yet
		add_point = true;
		//      Check the anzats
		for( vector<vector<Double_t> >::iterator ret_i = Returnable_Coord_Data.begin(); ret_i != Returnable_Coord_Data.end(); ++ret_i )
		{
			vector<bool> decisions;
			vector<Double_t>::iterator coord_i = input[i].begin();
			vector<Double_t>::iterator coord_j = ret_i->begin();
			for( ; coord_i != input[i].end(); ++coord_i, ++coord_j )
			{
				if( fabs( *coord_i - *coord_j ) < DT ) decisions.push_back( false );
				else decisions.push_back( true );
			}
			for( vector<bool>::iterator dec_i = decisions.begin(); dec_i != decisions.end(); ++dec_i )
			{
				add_point = add_point && *dec_i;
			}
		}
		//      If we haven't seen this point yet add it to the array of points
		if( add_point )
		{
			//cout << Coord_Data_pointer[i].first << "\t" << Coord_Data_pointer[i].second << endl;
			vector<Double_t> temp_vector;
			for( vector<Double_t>::iterator wanted_i = input[i].begin(); wanted_i != input[i].end(); ++wanted_i )
			{
				temp_vector.push_back( *wanted_i );
			}
			Returnable_Coord_Data.push_back( temp_vector );
		}
	}

	//      Return a 2D vector of unique corrdinates
	return Returnable_Coord_Data;
}

//      This will find all of the unique coordinates stored in a TPolyMarker3D object
//      NB this was written due to the fact that ROOT thows away the contents of {3/4}D TTree->Draw() objects... (God only knows why)
vector<vector<Double_t> > Mathematics::Unique_Coords( TPolyMarker3D *pm )
{
	//      Tolerance of the DOUBLE
	double DT = 1E-5;

	//      Get the number of coordinates in the TPolyMarker3D
	int number_of_points = pm->GetN();
	//      Get the data contained in the TPolyMarker3D
	Float_t* Coord_Data_pointer = pm->GetP();

	//      Will populate and return to the user
	vector<vector<Double_t> > Returnable_Coord_Data;

	//      Used within a for loop in logic, externally created/destroyed
	bool add_point = true;
	bool temp_decision_1 = true;
	bool temp_decision_2 = true;
	bool temp_decision_3 = true;

	//      Run over all of the points contained within the TPolyMarker3D object
	for( int i=0; i< number_of_points; ++i )
	{
		//      Assume we haven't seen this point yet
		add_point = true;
		//      Check the anzats
		for( unsigned int j=0; j< Returnable_Coord_Data.size(); ++j )
		{
			temp_decision_1 = true;
			temp_decision_2 = true;
			temp_decision_3 = true;
			if( fabs( Coord_Data_pointer[i*3] - Returnable_Coord_Data[j][0] ) < DT ) temp_decision_1 = false;
			if( fabs( Coord_Data_pointer[i*3+1] - Returnable_Coord_Data[j][1] ) < DT ) temp_decision_2 = false;
			if( fabs( Coord_Data_pointer[i*3+2] - Returnable_Coord_Data[j][2] ) < DT ) temp_decision_3 = false;

			//      If all of the information is the same do NOT add the point to the new vector of data
			if( ( ( temp_decision_1 == temp_decision_2 ) && ( temp_decision_2 == temp_decision_3 ) ) && ( temp_decision_1 == false ) )
				add_point = false;
		}
		//      If we haven't seen this point yet add it to the array of points
		if( add_point )
		{
			//cout << Coord_Data_pointer[i*3] << "\t" << Coord_Data_pointer[i*3+1] << "\t" << Coord_Data_pointer[i*3+2] << endl;
			vector<Double_t> temp_vector;
			temp_vector.push_back( Coord_Data_pointer[i*3] );
			temp_vector.push_back( Coord_Data_pointer[i*3+1] );
			temp_vector.push_back( Coord_Data_pointer[i*3+2] );
			Returnable_Coord_Data.push_back( temp_vector );
		}
	}

	//      Return a 2D vector of unique 3D corrdinates of size npoints*3
	return Returnable_Coord_Data;
}

//  The gamma distribution coded up in root is the more general form of that found on wikipedia (there's a surprise)
//
//  Using the root definition:                                  wiki:
//                              gamma = mean^2 / sigma^2                k     = mu^2 / sigma^2
//                              beta  = sigma^2 / mean                  theta = sigma^2 / mu
//
//              For:            mu == 0                         The 2 conditions above are ONLY valid for this condition

TF1* Mathematics::gamma_func( int OutputLevel )
{
	//      the gamma function in ROOT has issues with verbosity, let's silence it and then return the verbosity back at the end

	streambuf *nullbuf=NULL, *cout_bak=NULL, *cerr_bak=NULL, *clog_bak=NULL;
	ofstream filestr;
	filestr.open ("/dev/null");
	//      If the user wanted silence we point the Std Output Streams to /dev/null
	if( OutputLevel <= -1 )
	{
		cout_bak = cout.rdbuf();
		cerr_bak = cerr.rdbuf();
		clog_bak = clog.rdbuf();
		nullbuf = filestr.rdbuf();
		cout.rdbuf(nullbuf);
		cerr.rdbuf(nullbuf);
		clog.rdbuf(nullbuf);
	}

	TF1* output = new TF1( "gammaf", "[0]*TMath::GammaDist( x, ([1]*[1])/([2]*[2]), 0, ([2]*[2])/[1] )" );

	//  TF1 inherits from TFormula so we can use it's functions to rename the parameters to have consistancy with gaus / landau functions
	output->SetParName( 0, "Constant" );
	output->SetParName( 1, "Mean" );
	output->SetParName( 2, "Sigma" );

	//      Reset Std Output Streams
	if( OutputLevel <= -1 )
	{
		cout.rdbuf(cout_bak);
		cerr.rdbuf(cerr_bak);
		clog.rdbuf(clog_bak);
	}

	return output;
}

TF1* Mathematics::raw_gamma_func( int OutputLevel )
{
	(void) OutputLevel;
	TF1* output = new TF1( "gammaf", "[0]*TMath::GammaDist( x, [1], [2], [3] )" );

	output->SetParName( 0, "Const" );
	output->SetParName( 1, "#gamma" );
	output->SetParName( 2, "#mu" );
	output->SetParName( 3, "#beta" );

	return output;
}

TF1* Mathematics::landau_func()
{
	//  Want plotting consistancy so have defined my own landau function
	TF1 *land = new TF1( "mylandau", "[0]*TMath::Landau( x, [1], [2] )" );
	land->SetParName( 0, "Constant" );
	land->SetParName( 1, "Mean" );
	land->SetParName( 2, "Sigma" );
	return land;
}

//	Return if the first of each of these pairs is unique
bool Mathematics::Sort_first_Double( pair<double,double> one_pair, pair<double,double> two_pair )
{
	bool x_larger = one_pair.first < two_pair.first ;
	if( x_larger ) return true;
	else return false;
}

//	Return if the second of each of these pairs is unique
bool Mathematics::Sort_second_Double( pair<double,double> one_pair, pair<double,double> two_pair )
{
	bool y_larger = one_pair.second < two_pair.second ;
	if( y_larger ) return true;
	else return false;
}

//	Return if these are a unique pair of doubles or not
bool Mathematics::Unique_2D_Double( pair<double,double> one_pair, pair<double,double> two_pair )
{
	bool x_same = fabs(one_pair.first - two_pair.first) < 1E-3 ;
	bool y_same = fabs(one_pair.second - two_pair.second) < 1E-3 ;
	return x_same && y_same;
}
