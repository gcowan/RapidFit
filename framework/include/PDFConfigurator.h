/**
 @class PDFConfigurator

 This is a class to contain configuration information for a PDF
 The information is taken from the XML  under the <PDF> tag
 Objects of this class are passed to PDF constructors to configure them

 @author Pete Clarke
 @date 2011-05-01
*/
#ifndef PDFCONFIGURATOR_H
#define PDFCONFIGURATOR_H

#include <string>
#include <vector>

using namespace std;

class PDFConfigurator
{
	public:
		PDFConfigurator();
		~PDFConfigurator();

		// Methods for use when substituting a parameter name for a different name
		// E.g. to substitute tau_LL  ->  tauLL2
		void addParameterSubstitution( string substitution );
		string getName( string defaultName ) ; 
	
	    // Set configuration parameter
		void addConfigurationParameter( string configString ) ;
		string getConfigurationValue( string configParam );
		bool hasConfigurationValue( string configParam, string val );
		bool isTrue( string configParam );
	    

	private:
		vector<string> defaultNames;
		vector<string> replacementNames;

		vector<string> configParameters;
		vector<string> configValues;
};

#endif
