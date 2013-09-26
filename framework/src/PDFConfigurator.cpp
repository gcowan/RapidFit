/**
  @class PDFConfigurator

  This is a class to contain configuration information for a PDF
  The information is taken from the XML  under the <PDF> tag
  Objects of this class are passed to PDF constructors to configure them

  @author Pete Clarke
  @date 2011-05-01
 */

#include "PDFConfigurator.h"
#include "StringProcessing.h"
#include "ClassLookUp.h"
#include "IPDF.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

using namespace::std;

///..........................................
//Default constructor & destructor
PDFConfigurator::PDFConfigurator() : defaultNames(), replacementNames(), configParameters(), configValues(), daughterPDFs(), thisPDFsBoundary(NULL), fractionName("UnsetFractionName"), PDFLabel("")
{
}

PDFConfigurator::PDFConfigurator( const PDFConfigurator& input ) :
	defaultNames( input.defaultNames ), replacementNames( input.replacementNames ), configParameters( input.configParameters ), configValues( input.configValues ),
	daughterPDFs(), thisPDFsBoundary(NULL), fractionName( input.fractionName ), PDFLabel(input.PDFLabel)
{
	for( unsigned int i=0; i< input.daughterPDFs.size(); ++i )
	{
		daughterPDFs.push_back( ClassLookUp::CopyPDF( input.daughterPDFs[i] ) );
	}
	if( input.thisPDFsBoundary != NULL ) thisPDFsBoundary = new PhaseSpaceBoundary( *(input.thisPDFsBoundary) );
}

PDFConfigurator::~PDFConfigurator()
{
	if( thisPDFsBoundary != NULL ) delete thisPDFsBoundary;
	while( !daughterPDFs.empty() )
	{
		if( daughterPDFs.back() != NULL ) delete daughterPDFs.back();
		daughterPDFs.pop_back();
	}
}

void PDFConfigurator::AddDaughterPDF( const IPDF* input )
{
	daughterPDFs.push_back( ClassLookUp::CopyPDF( input ) );
}

vector<IPDF*> PDFConfigurator::GetDaughterPDFs() const
{
	return daughterPDFs;
}

void PDFConfigurator::SetPhaseSpaceBoundary( PhaseSpaceBoundary* input )
{
	if( thisPDFsBoundary != NULL ) delete thisPDFsBoundary;
	thisPDFsBoundary = new PhaseSpaceBoundary( *input );
}

PhaseSpaceBoundary* PDFConfigurator::GetPhaseSpaceBoundary() const
{
	return thisPDFsBoundary;
}

void PDFConfigurator::appendParameterNames( string names_list )
{
	vector<string> all_params = StringProcessing::SplitString( names_list, ':' );

	if( all_params.size() < 2 )
	{
		cerr << "Cannot Understand AppendParameterNames:\t" << names_list << endl;
		exit(-5434);
	}

	string suffix = all_params.back();

	for( unsigned int i=0; i< (all_params.size()-1); ++i )
	{
		string old_name = all_params[i];
		string new_name = old_name;
		new_name.append( suffix );
		string substitution = old_name;
		substitution.append(":");
		substitution.append( new_name );
		this->addParameterSubstitution( substitution );
	}
}

//...........................................
// Parameter substitution methods
// These allow you to configure a PDF with a substituttion of a default parameter name for another name
// Example :  substitute mass -> massBd

// Method to store a substitution
void PDFConfigurator::addParameterSubstitution( string substitution ) {

	// The string substitution is of the form defaultName:replacementName

	// Find position of ":" character
	size_t pos = substitution.find_first_of(":") ;
	if( pos == string::npos )
	{
		cout << "  In PDFConfigurator::addParameterSubstitution : No separator found  in string : " << substitution << endl ;
		exit(1);
	}

	// Break into 2 strings
	string defaultName     =  substitution.substr( 0, pos ) ;
	string replacementName =  substitution.substr( pos+1, string::npos ) ;
	//cout << " Sub[" << defaultName << "->" << replacementName << "]" << endl;//setw(3) << " ";
	defaultNames.push_back( defaultName );
	replacementNames.push_back( replacementName );

	return ;
}

// Method to return a substitute name against an input default name.
// If no substitution then the input name is returned unchanged
string PDFConfigurator::getName( string defaultName )
{
	for( unsigned int ii=0; ii<defaultNames.size() ; ++ii )
	{
		//cout << " Testing " << defaultName <<  "   against   " << defaultNames[ii] << endl ;
		if( defaultName == defaultNames[ii] )
		{
			//cout << " PDFConfig::sub[" << defaultName << "->" << replacementNames[ii] << "]" << endl;
			return replacementNames[ii] ;
		}
	}
	return defaultName ;
}


//............................................
// Configuration parameters

// Method to store a configuration parmeter
// The configuration parameter is of the form parameter:value
void PDFConfigurator::addConfigurationParameter( string configString )
{
	// Find position of ":" character
	size_t pos = configString.find_first_of(":") ;
	if( pos == string::npos )
	{
		cout << "  In PDFConfigurator::addConfigurationParameter : No separator found  in string : " << configString << endl ;
		exit(1);
	}

	// Break into 2 strings
	string configParameter  =  configString.substr( 0, pos ) ;
	string configValue      =  configString.substr( pos+1, string::npos ) ;
	//cout << " ConfigParam[" << configParameter << " : " << configValue << "]" << endl;//setw(3) << " ";
	configParameters.push_back( configParameter ) ;
	configValues.push_back( configValue ) ;

	return ;
}


// Method to return a configuration parameter value
// If not found then returns " "
string PDFConfigurator::getConfigurationValue( string configParam )
{
	for( unsigned int ii=0; ii<configParameters.size() ; ++ii )
	{
		//cout << " Testing " << configParam <<  "   against   " << configParameters[ii] << endl ;
		if( configParam == configParameters[ii] )
		{
			//cout << " PDFConfiguratorgetConfigurationValue setting [" << configParam << "--->" << configValues[ii] << "]" << endl ;
			return configValues[ii] ;
		}
	}
	return string("") ;
}

// Method to check for a configuration parameter value
bool PDFConfigurator::hasConfigurationValue( string configParam, string paramValue )
{
	for( unsigned int ii=0; ii<configParameters.size() ; ++ii )
	{
		//cout << " Testing " << configParam <<  "   against   " << configParameters[ii] << endl ;
		if( configParam == configParameters[ii] )
		{
			//cout << " PDFConfigurator::hasConfigurationValue found [" << configParam << " == " << configValues[ii] << "]" << endl ;
			return ( configValues[ii] == paramValue ) ;
		}
	}
	return false ;
}

// Method to check for a configuration parameter value boolean
bool PDFConfigurator::isTrue( string configParam )
{
	for( unsigned int ii=0; ii<configParameters.size() ; ++ii )
	{
		//cout << " Testing " << configParam <<  "   against   " << configParameters[ii] << endl ;
		if( configParam == configParameters[ii] )
		{
			bool returnValue =  ( strcasecmp( configValues[ii].c_str() , "TRUE" ) == 0 )  ;
			//cout << " PDFConfigurator::isTrue found [" << configParam << " == " << returnValue << " ]" << endl ;
			return returnValue ;
		}
	}
	return false ;
}

string PDFConfigurator::XML() const
{
	stringstream xml;
	for( unsigned int i=0; i< configParameters.size(); ++i )
	{
		xml << "\t<ConfigurationParameter>" << configParameters[i] << ":" << configValues[i] << "</ConfigurationParameter>" << endl;
	}
	for( unsigned int i=0; i< defaultNames.size(); ++i )
	{
		xml << "\t<ParameterSubstitution>" << defaultNames[i] << ":" << replacementNames[i] << "</ParameterSubstitution>" << endl;
	}
	return xml.str();
}

void PDFConfigurator::Print() const
{
	for( unsigned int i=0; i< configParameters.size() ; ++i )
	{
		cout << " PDFConfigurator::Param   [ " << configParameters[i] << " : " << configValues[i] << " ]" << endl;
	}
	for( unsigned int i=0; i< defaultNames.size(); ++i )
	{
		cout << " PDFConfigurator::Sub     [ " << defaultNames[i] << " ---> " << replacementNames[i] << " ]" << endl;
	}
	return;
}

bool PDFConfigurator::empty() const
{
	return configParameters.empty() && defaultNames.empty();
}

void PDFConfigurator::SetFractionName( const string thisName )
{
	fractionName = thisName;
}

string PDFConfigurator::GetFractionName() const
{
	return fractionName;
}

void PDFConfigurator::SetPDFLabel( const string input )
{
	PDFLabel = input;
}

string PDFConfigurator::GetPDFLabel() const
{
	return PDFLabel;
}

