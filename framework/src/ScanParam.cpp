/**
  @class ScanParam
      an object to hold the 

  @author Rob Currie
  @date 2011-02
  */

#include "ScanParam.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cstdio>

using namespace::std;

ScanParam::ScanParam()
{
}

ScanParam::ScanParam( vector<string> nName, vector<string> nType, vector<double> nMinimum, vector<double> nMaximum, vector<int> nSigma, vector<int> nPoints) :
	name(nName), type(nType), minimum(nMinimum), maximum(nMaximum), sigma(nSigma), points(nPoints)
{
	if( !minimum.empty() && !maximum.empty() )
	{
		if( minimum[0] < maximum[0] )
		{
			cerr << "Scan Limit \"" << nName[0] << "\" has maximum less than minimum: values swapped" << endl;
			minimum.swap( maximum );
		}
		if( (minimum[0] - maximum[0]) < 1E-6 )
		{
			cerr << "Scan Limit \"" << nName[0] << "\" has maximum = minimum, assuming you forgot a sign in the XML" << endl;
			minimum[0]=-minimum[0];
		}
	}
}

ScanParam::ScanParam( string nName, string nType, double nMaximum, double nMinimum, int nPoints )
{
	name.push_back( nName );
	type.push_back( nType );
	points.push_back( nPoints );
	if ( maximum < minimum )
	{
		cerr << "Scan Limit \"" << nName << "\" has maximum less than minimum: values swapped" << endl;
		minimum.push_back( nMaximum );
		maximum.push_back( nMinimum );
	} else {
		minimum.push_back( nMinimum );
		maximum.push_back( nMaximum );
	}
}

ScanParam::ScanParam( string nName, string nType, int nSigma, int nPoints )
{
	name.push_back( nName );
	type.push_back( nType );
	points.push_back( nPoints );
	sigma.push_back( nSigma );
}

ScanParam::ScanParam( string nName, string nType, int nPoints)
{
	name.push_back( nName );
	type.push_back( nType );
	points.push_back( nPoints );
}

ScanParam::ScanParam( string nName )
{
	name.push_back( nName );
}

//Destructor
ScanParam::~ScanParam()
{
}

//General print
void ScanParam::print()
{
	if( !name.empty() )     cout << "   name        " << name[0] << endl ;
	if( !type.empty() )     cout << "   type        " << type[0] << endl;
	if( !minimum.empty() )  cout << "   minimum     " << minimum[0] << endl ;
	if( !maximum.empty() )  cout << "   maximum     " << maximum[0] << endl ;
	if( !sigma.empty() )    cout << "   sigma       " << sigma[0] << endl ;
	if( !points.empty() )   cout << "   points      " << points[0] << endl ;
}
