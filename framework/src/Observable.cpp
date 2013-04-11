/**
  @class Observable

  The measured value of a particular variable for an event

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
 */

//	RapidFit Headers
#include "Observable.h"
//	System Headers
#include <iostream>
#include <vector>
#include <string>

using namespace::std;

//Constructor with correct argument
Observable::Observable( const string Name, const double NewValue, const string NewUnit )
	: name(Name), value(NewValue), unit(NewUnit), bin_num(-1), acceptance(-1.), bkg_bin_num(-1), bkg_acceptance(-1.), offset(-999.)
{
	if (unit == "")
	{
		cerr << "Observable \"" << Name << "\" has no unit! What kind of physicist are you?" << endl;
	}
}

Observable::Observable( const string Name ) : name(Name), value(0.), unit("Uninitialised"), bin_num(-1), acceptance(-1.), bkg_bin_num(-1), bkg_acceptance(-1.), offset(-999.)
{
}

Observable::Observable( const Observable& input ) :
	name(input.name), value(input.value), unit(input.unit), bin_num(input.bin_num),
	acceptance(input.acceptance), bkg_bin_num(input.bkg_bin_num), bkg_acceptance(input.bkg_acceptance), offset(input.offset)
{
}

//Destructor
Observable::~Observable()
{
}

//Get and set the value
double Observable::GetValue() const
{
	return value;
}

//Get the unit
string Observable::GetUnit() const
{
	return unit;
}

void Observable::SetBinNumber( const int input ) const
{
	bin_num = input;
}

void Observable::SetBkgBinNumber( const int input ) const
{
	bkg_bin_num = input;
}

int Observable::GetBinNumber() const
{
	return bin_num;
}

int Observable::GetBkgBinNumber() const
{
	return bkg_bin_num;
}

void Observable::SetAcceptance( const double input ) const
{
	acceptance = input;
}

void Observable::SetBkgAcceptance( const double input ) const
{
	bkg_acceptance = input;
}

double Observable::GetAcceptance() const
{
	return acceptance;
}

double Observable::GetBkgAcceptance() const
{
	return bkg_acceptance;
}

string Observable::GetName() const
{
	return name;
}

void Observable::Print() const
{
	cout << "Name: " << name;
	cout << "\tValue: " << value;
	cout << "\tUnit: " << unit;
	if( bin_num != -1 ) cout << "\tBinNum: " << bin_num << "\tAcceptance: " << acceptance;
	if( bkg_bin_num != -1 ) cout << "\tBkgBinNum: " << bkg_bin_num << "\tBkgAcceptance: " << bkg_acceptance;
	cout << endl;
}

void Observable::SetObservable( const Observable* input )
{
	name = input->GetName();
	value = input->GetValue();
	unit = input->GetUnit();
	bin_num = input->GetBinNumber();
	acceptance = input->GetAcceptance();
	bkg_bin_num = input->GetBkgBinNumber();
	bkg_acceptance = input->GetBkgAcceptance();
	offset = input->GetOffSet();
}

void Observable::ExternallySetValue( const double input )
{
	value = input;
}

void Observable::SetOffSet( const double input ) const
{
	offset = input;
}

double Observable::GetOffSet() const
{
	return offset;
}

