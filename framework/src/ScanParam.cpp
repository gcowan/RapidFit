/**
  @class ScanParam
      an object to hold the 

  @author Rob Currie
  @date 2011-02
  */

//	RapidFit Headers
#include "ScanParam.h"
//	System Headers
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cstdio>

using namespace::std;

/*ScanParam::ScanParam() : name(), minimum(), maximum(), sigma(), points()
{
}*/

ScanParam::ScanParam( vector<string> nName, vector<double> nMinimum, vector<double> nMaximum, vector<int> nSigma, vector<int> nPoints) :
	name(nName), minimum(nMinimum), maximum(nMaximum), sigma(nSigma), points(nPoints)
{
	if( !minimum.empty() && !maximum.empty() )
	{
		if( minimum[0] < maximum[0] )
		{
			cerr << "Scan Limit \"" << nName[0] << "\" has maximum less than minimum: values swapped" << endl;
			minimum.swap( maximum );
		}
	}
}


bool ScanParam::HasName() {  return !name.empty();  }
string ScanParam::GetName() {
	if( !name.empty() ){  return name[0];  }
	else return "";
}

void ScanParam::SetName( string new_val ) {
	while( !name.empty() ){   name.pop_back(); }
	name.push_back(    new_val );
}

bool ScanParam::HasMax() { return !maximum.empty(); }

double ScanParam::GetMax() {
	if( !maximum.empty() ){ return maximum[0]; }
	else return 0;
}

void ScanParam::SetMax(double new_val) {
	while( !maximum.empty() ){
		maximum.pop_back(); };
		maximum.push_back( new_val );
}

bool ScanParam::HasMin() {  return !minimum.empty();  }

double ScanParam::GetMin() {
	if( !minimum.empty() ){ return minimum[0]; }
	else return 0;
}

void ScanParam::SetMin(double new_val) {
	while( !minimum.empty() ){ minimum.pop_back(); };
	minimum.push_back( new_val );
}

bool ScanParam::HasSigma(){  return !sigma.empty();  }

int ScanParam::GetSigma() {
	if( !sigma.empty() ) {   return sigma[0]; }
	else return 5;
}

void ScanParam::SetSigma(int new_val) {
	while( !sigma.empty() ){    sigma.pop_back(); };
	sigma.push_back( new_val );
}

bool ScanParam::HasPoints(){ return !points.empty(); }

int ScanParam::GetPoints() {
	if( !points.empty() ){  return points[0];  }
	else return 10;
}

void ScanParam::SetPoints(int new_val) {
	while( !points.empty() ){  points.pop_back(); };
	points.push_back(  new_val );
}


ScanParam::ScanParam( string nName, double nMaximum, double nMinimum, int nPoints ) : name(), minimum(), maximum(), sigma(), points()
{
	name.push_back( nName );
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

ScanParam::ScanParam( string nName, int nSigma, int nPoints ) : name(), minimum(), maximum(), sigma(), points()
{
	name.push_back( nName );
	points.push_back( nPoints );
	sigma.push_back( nSigma );
}

ScanParam::ScanParam( string nName, int nPoints) : name(), minimum(), maximum(), sigma(), points()
{
	name.push_back( nName );
	points.push_back( nPoints );
}

ScanParam::ScanParam( string nName ) : name(), minimum(), maximum(), sigma(), points()
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
	if( !minimum.empty() )  cout << "   minimum     " << minimum[0] << endl ;
	if( !maximum.empty() )  cout << "   maximum     " << maximum[0] << endl ;
	if( !sigma.empty() )    cout << "   sigma       " << sigma[0] << endl ;
	if( !points.empty() )   cout << "   points      " << points[0] << endl ;
}
