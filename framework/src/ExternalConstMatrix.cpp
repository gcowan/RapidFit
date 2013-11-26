/**
  @class ExternalConstMatrix

  A class that holds experimentally dervied constraints on fit parameters

  @author Benjamin M Wynne bwynne@cern.ch
  @date 21-01-10
  */

#include "TMatrixT.h"
///	RapidFit Headers
#include "ExternalConstMatrix.h"
#include "StringProcessing.h"
///	System Headers
#include <string>
#include <sstream>
#include <iostream>

using namespace::std;

//Constructor with correct arguments
ExternalConstMatrix::ExternalConstMatrix( string NewName, string NewValue, string NewError, string NewCorrelations ) :
	names(NewName), values(NewValue), errors(NewError), correlations(NewCorrelations), internalParameterSet(NULL), wantedParameters(),
	names_val(), values_val(), errors_val(), corr_matrix()
{
	vector<string> allNames = StringProcessing::SplitString( names, ':' ); wantedParameters = allNames;
	vector<string> allValues = StringProcessing::SplitString( values, ':' );
	vector<string> allErrors = StringProcessing::SplitString( errors, ':' );
	vector<string> allCorrs = StringProcessing::SplitString( correlations, ':' );

	if( ( allNames.size() != allValues.size() ) || ( allNames.size() != allErrors.size() ) ||
			( allNames.size()*allNames.size() != allCorrs.size() ) )
	{
		cerr << "Wrong number of arguments seperated by ':' provided for ExternalCorrMatrix" << endl;
		exit(38740912);
	}

	for( unsigned int i=0; i< allValues.size(); ++i )
	{
		values_val.push_back( strtod( allValues[i].c_str(), NULL ) );
		errors_val.push_back( strtod( allErrors[i].c_str(), NULL ) );
	}

	vector<double>flatCorrs;
	for( unsigned int i=0; i< allCorrs.size(); ++i )
	{
		flatCorrs.push_back( strtod( allCorrs[i].c_str(), NULL ) );
	}

	for( unsigned int i=0; i< allNames.size(); ++i )
	{
		vector<double> thisRow;
		for( unsigned int j=0; i< allNames.size(); ++j )
		{
			thisRow.push_back( flatCorrs[i*allNames.size()+j] );
		}
		corr_matrix.push_back( thisRow );
	}

	internalParameterSet = new ParameterSet( allNames );
}

ExternalConstMatrix::ExternalConstMatrix( const ExternalConstMatrix& input ) :
	names(input.names), values(input.values), errors(input.errors), correlations(input.correlations), internalParameterSet(NULL), wantedParameters(input.wantedParameters),
	names_val(input.names_val), values_val(input.values_val), errors_val(input.errors_val), corr_matrix(input.corr_matrix)
{
	if( input.internalParameterSet != NULL )
	{
		internalParameterSet = new ParameterSet( *(input.internalParameterSet) );
	}
}
//Destructor
ExternalConstMatrix::~ExternalConstMatrix()
{
	if( internalParameterSet != NULL ) delete internalParameterSet;
}

//Return the information held
string ExternalConstMatrix::GetName() const
{
	return names;
}


void ExternalConstMatrix::SetPhysicsParameters( const ParameterSet* input )
{
	internalParameterSet->SetPhysicsParameters( input );
}

bool ExternalConstMatrix::CanApply( const ParameterSet* input ) const
{
	vector<string> input_names = input->GetAllNames();
	vector<int> wanted_locations;
	for( unsigned int i=0; i< wantedParameters.size(); ++i )
	{
		wanted_locations.push_back( StringProcessing::VectorContains( &input_names, &(wantedParameters[i]) ) );
	}
	bool decision = true;
	for( unsigned int i=0; i< wanted_locations.size(); ++i )
	{
		bool thisdecision = wanted_locations[i] != -1;
		decision = decision && thisdecision;
	}
	if( !decision )
	{
		cerr << "Cannot Apply " << names << " constraint" << endl;
	}
	return decision;
}

double ExternalConstMatrix::GetChi2() const
{
	double returnable=0.;

	int dim = (int)names.size();

	TMatrixT<double> diff( dim, 1 );
	TMatrixT<double> diff_T( 1, dim );

	for( unsigned int i=0; i< (unsigned)dim; ++i )
	{
		diff( i, 0 ) =  internalParameterSet->GetPhysicsParameter( wantedParameters[i] )->GetValue() - values_val[i];
		diff_T( 0, i ) = diff( i, 0 );
	}

	TMatrixT<double> cov_matrix( dim, dim );
	for( unsigned int i=0; i< (unsigned)dim; ++i )
	{
		for( unsigned int j=0; j< (unsigned)dim; ++j )
		{
			cov_matrix( i, j ) =  corr_matrix[i][j]*errors_val[i]*errors_val[j];
		}
	}

	TMatrixT<double> finalMatrix( 1, 1 );

	finalMatrix = diff * cov_matrix * diff_T;

	returnable = finalMatrix( 0, 0 );

	return returnable;
}

void ExternalConstMatrix::Print() const
{
	cout << "External Constrint Matrix: " << names << "\tValues: " << values << "\tErrors: " << errors << endl;
}

string ExternalConstMatrix::XML() const
{
	stringstream xml;

	xml << "<ToFit>" << endl;
	xml << "\t<ConstraintFunction>" << endl;
	xml << "\t\t<ExternalConstMatrix>" << endl;
	xml << "\t\t\t<Names>" << names << "</Names>" << endl;
	xml << "\t\t\t<Values>" << values << "</Values>" << endl;
	xml << "\t\t\t<Errors>" << errors << "</Errors>" << endl;
	xml << "\t\t\t<Correlations>" << correlations << "</Correlations>" << endl;
	xml << "\t\t</ExternalConstMatrix>" << endl;
	xml << "\t</ConstraintFunction>" << endl;
	xml << "</ToFit>" << endl;

	return xml.str();
}

bool ExternalConstMatrix::isExternalConst() const
{
	return false;
}

bool ExternalConstMatrix::isExternalMatrix() const
{
	return true;
}
