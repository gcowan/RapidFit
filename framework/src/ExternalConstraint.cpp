/**
  @class ExternalConstraint

  A class that holds experimentally dervied constraints on fit parameters

  @author Benjamin M Wynne bwynne@cern.ch
  @date 21-01-10
  */

///	RapidFit Headers
#include "ExternalConstraint.h"
#include "StringProcessing.h"
///	System Headers
#include <string>
#include <sstream>
#include <iostream>

using namespace::std;

//Constructor with correct arguments
ExternalConstraint::ExternalConstraint( string NewName, double NewValue, double NewError ) : name(NewName), value(NewValue), error(NewError), internalParameterSet(NULL), wantedParameters()
{
	if( (name == "GammaL") || (name == "GammaH") || (name == "GammaObs") || (name == "GLandGH") ||
			(name == "SpecialForGreig") || (name == "Gammas-DeltaGamma-1fbPaper") || (name == "Gammas-DeltaGamma-1fbPaper-combined") )
	{
		wantedParameters.push_back("gamma");
		wantedParameters.push_back("deltaGamma");
	}
	else if( name == "CosSqPlusSinSq" )
	{
		wantedParameters.push_back( "cosphis" );
		wantedParameters.push_back( "sinphis" );
	}
	else if( name == "ATOTAL" )
	{
		wantedParameters.push_back( "Aperp_sq" );
		wantedParameters.push_back( "Azero_sq" );
	}
	else
	{
		wantedParameters.push_back( name );
	}


	internalParameterSet = new ParameterSet( wantedParameters );
}

string ExternalConstraint::GetValueStr() const
{
	stringstream this_stream;
	this_stream << value;
	return this_stream.str();
}

string ExternalConstraint::GetErrorStr() const
{
	stringstream this_stream;
	this_stream << error;
	return this_stream.str();
}

ExternalConstraint::ExternalConstraint( const ExternalConstraint& input ) :
	name( input.name ), value( input.value ), error( input.error ), internalParameterSet(NULL), wantedParameters(input.wantedParameters)
{
	if( input.internalParameterSet != NULL )
	{
		internalParameterSet = new ParameterSet( *(input.internalParameterSet) );
	}
}

//Destructor
ExternalConstraint::~ExternalConstraint()
{
	if( internalParameterSet != NULL ) delete internalParameterSet;
}

//Return the information held
string ExternalConstraint::GetName() const
{
	return name;
}


void ExternalConstraint::SetPhysicsParameters( const ParameterSet* input )
{
	internalParameterSet->SetPhysicsParameters( input );
}

bool ExternalConstraint::CanApply( const ParameterSet* input ) const
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
		cerr << "Cannot Apply " << name << " constraint" << endl;
	}
	return decision;
}

double ExternalConstraint::GetChi2() const
{
	double returnable=0.;
	if( (name == "GammaL") || (name == "GammaH") || (name == "GammaObs") || (name == "GLandGH") ||
			(name == "SpecialForGreig") || (name == "Gammas-DeltaGamma-1fbPaper") || (name == "Gammas-DeltaGamma-1fbPaper-combined") )
	{
		double gammaVal = internalParameterSet->GetPhysicsParameter( "gamma" )->GetValue();
		double deltaGammaVal = internalParameterSet->GetPhysicsParameter( "deltaGamma" )->GetValue();
		if( name == "GammaL" )
		{
			double gaml_fit = gammaVal + (deltaGammaVal * 0.5);
			double gaussSqrt = (gaml_fit - value) / error;
			returnable += gaussSqrt * gaussSqrt;
		}
		else if( name == "GammaH" )
		{
			double gaml_fit = gammaVal - (deltaGammaVal * 0.5);
			double gaussSqrt = (gaml_fit - value) / error;
			returnable += gaussSqrt * gaussSqrt;
		}
		else if( name == "GLandGH" )
		{
			double G1 = gammaVal - (deltaGammaVal * 0.5);
			double G2 = gammaVal + (deltaGammaVal * 0.5);

			if( G1 < 0.) returnable += (G1 * G1) / (error*error);
			if( G2 < 0. ) returnable += (G2 * G2) / (error*error);
		}
		else if( name == "GammaObs" )
		{
			double gamobs_fit = gammaVal * (1. - (deltaGammaVal / 2. / gammaVal) * (deltaGammaVal / 2. / gammaVal))
				/ (1. + (deltaGammaVal / 2. / gammaVal) * (deltaGammaVal / 2. / gammaVal));
			double gaussSqrt = (gamobs_fit - value) / error;
			returnable += gaussSqrt * gaussSqrt;
		}
		else if( name == "SpecialForGreig" )
		{
			double gamma_con = 0.657;
			double dgam_con = 0.123;
			double v[2];
			v[0] = gammaVal - gamma_con;
			v[1] = deltaGammaVal- dgam_con;
			double E[2][2];
			E[0][0] = 7213.6; // 7213.6;  //13566.7
			E[1][1] = 1087.3; // 1087.3;  //1221.
			E[1][0] = 587.1;  // 587.1;   //1221.
			E[0][1] = E[1][0];
			double chisq=0.;
			for (int ix = 0; ix < 2; ++ix) {
				for (int iy = 0; iy < 2; ++iy) {
					chisq += v[ix] * E[ix][iy] * v[iy];
				}
			}
			returnable += chisq;
		}
		else if( name == "Gammas-DeltaGamma-1fbPaper" )
		{
			double gamma_con = 0.6631;
			double dgam_con = 0.100;
			double v[2];
			v[0] = gammaVal - gamma_con;
			v[1] = deltaGammaVal - dgam_con;
			double E[2][2];
			E[0][0] = 16850.7;
			E[1][1] = 3988.85;
			E[1][0] = 1904.58;
			E[0][1] = E[1][0];
			double chisq = 0.;
			for (int ix = 0; ix < 2; ++ix) {
				for (int iy = 0; iy < 2; ++iy) {
					chisq += v[ix] * E[ix][iy] * v[iy];
				}
			}
			returnable += chisq;
		}
		else if( name == "Gammas-DeltaGamma-1fbPaper-combined" )
		{
			double gamma_con = 0.661;
			double dgam_con = 0.106;
			double v[2];
			v[0] = gammaVal - gamma_con;
			v[1] = deltaGammaVal - dgam_con;
			double E[2][2];
			E[0][0] = 19273.;
			E[1][1] = 5895.26;
			E[1][0] = -498.83;
			E[0][1] = E[1][0];
			double chisq = 0.;
			for (int ix = 0; ix < 2; ++ix) {
				for (int iy = 0; iy < 2; ++iy) {
					chisq += v[ix] * E[ix][iy] * v[iy];
				}
			}
			returnable += chisq;
		}
	}
	else if( name == "ATOTAL" )
	{
            double Aperpsq = internalParameterSet->GetPhysicsParameter( "Aperp_sq" )->GetValue();
            double Azerosq = internalParameterSet->GetPhysicsParameter( "Azero_sq" )->GetValue();
            // double Assq =  internalParameterSet->GetPhysicsParameter(AsName)->GetValue();
            double excess = ( value - Aperpsq - Azerosq); //-Assq) ;
            double penalty=0.;
            if (excess >= 0) penalty = 0;
            else penalty = ( excess * excess ) / ( error*error );
            returnable += penalty;
	}
	else
	{
		//if( name == "mistagP1_OS" )
		//{
//			cout << endl;
//			cout << name << ": " << internalParameterSet->GetPhysicsParameter(name)->GetValue() << endl;
//			double parameterValuea = internalParameterSet->GetPhysicsParameter(name)->GetValue();
//			cout << "value: " << value << endl;
//			cout << "error: " << error << endl;
//			cout << "Chi2:  " << ((parameterValuea - value ) / error)*((parameterValuea - value ) / error) << endl;
//			cout << endl;
		//}
		double parameterValue = internalParameterSet->GetPhysicsParameter(name)->GetValue();
		double gaussSqrt = (parameterValue - value ) / error;
		returnable += gaussSqrt * gaussSqrt;
	}
	return returnable;
}

double ExternalConstraint::GetValue() const
{
	return value;
}

double ExternalConstraint::GetError() const
{
	return error;
}

void ExternalConstraint::Print() const
{
	cout << "External Constrint: " << name << "\tValue: " << value << "\tError: " << error << endl;
}

string ExternalConstraint::XML() const
{
	stringstream xml;

	xml << "<ToFit>" << endl;
	xml << "\t<ConstraintFunction>" << endl;
	xml << "\t\t<ExternalConstraint>" << endl;
	xml << "\t\t\t<Name>" << name << "</Name>" << endl;
	xml << "\t\t\t<Value>" << value << "</Value>" << endl;
	xml << "\t\t\t<Error>" << error << "</Error>" << endl;
	xml << "\t\t</ExternalConstraint>" << endl;
	xml << "\t</ConstraintFunction>" << endl;
	xml << "</ToFit>" << endl;

	return xml.str();
}

bool ExternalConstraint::isExternalConst() const
{
	return true;
}

bool ExternalConstraint::isExternalMatrix() const
{
	return false;
}
