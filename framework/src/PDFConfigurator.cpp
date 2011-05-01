/**
 @class PDFConfigurator

 This is a class to contain configuration information for a PDF
 The information is taken from the XML  under the <PDF> tag
 Objects of this class are passed to PDF constructors to configure them
 
 @author Pete Clarke
 @date 2011-05-01
 */

#include <iostream>
#include <cstdlib>
#include "PDFConfigurator.h"

///..........................................
//Default constructor & destructor
PDFConfigurator::PDFConfigurator() { }
PDFConfigurator::~PDFConfigurator() { }

//...........................................
// Parameter substitution methods
// These allow you to configure a PDF with a substituttion of a default parameter name for another name
// Example :  substitute mass -> massBd 

// Method to store a substitution
void PDFConfigurator::addParameterSubstitution( string substitution ) {

	// The string substitution is of the form defaultName:replacementName

	// Find position of ":" character
	size_t pos = substitution.find_first_of(":") ;
	if( pos == string::npos ) {
		cout << "  In PDFConfigurator::addParameterSubstitution : No separator found  in string : " << substitution << endl ;
		exit(1);
	}
	
	// Break into 2 strings
	string defaultName     =  substitution.substr( 0, pos ) ;
	string replacementName =  substitution.substr( pos+1, string::npos ) ;
	cout << " In PDFConfigurator::addParameterSubstitution : storing substitution [" << defaultName << "--->" << replacementName << "]" << endl ;
	defaultNames.push_back( defaultName ) ;
	replacementNames.push_back( replacementName ) ;
	
	return ;
}

// Method to return a substitute name against an input default name.
// If no substitution then the input name is returned unchanged
string PDFConfigurator::getName( string defaultName ) {
	for( unsigned int ii=0; ii<defaultNames.size() ; ii++ ) {
		//cout << " Testing " << defaultName <<  "   against   " << defaultNames[ii] << endl ;
		if( defaultName == defaultNames[ii] ) {
			cout << " PDFConfigurator substituting [" << defaultName << "--->" << replacementNames[ii] << "]" << endl ;
			return replacementNames[ii] ;
		}
	}
	return defaultName ;
};

		   
	
