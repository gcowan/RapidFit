/**
  @class PDFConfigurator

  This is a class to contain configuration information for a PDF
  The information is taken from the XML  under the <PDF> tag
  Objects of this class are passed to PDF constructors to configure them

  @author Pete Clarke
  @date 2011-05-01
 */

#pragma once
#ifndef PDFCONFIGURATOR_H
#define PDFCONFIGURATOR_H

//	RapidFit Headers
#include "IDataSet.h"
#include "PhaseSpaceBoundary.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class IDataSet;
class IPDF;

class PDFConfigurator
{
	public:
		/*!
		 * @brief Constructor
		 */
		PDFConfigurator();

		/*!
		 * @brief Copy Constructor
		 */
		PDFConfigurator( const PDFConfigurator& );

		/*!
		 * @brief Destructor
		 */
		~PDFConfigurator();

		// Methods for use when substituting a parameter name for a different name
		// E.g. to substitute tau_LL  ->  tauLL2
		void appendParameterNames( string names_list );
		void addParameterSubstitution( string substitution );
		string getName( string defaultName ) ; 

		// Set configuration parameter
		/*!
		 * @brief
		 *
		 * @return void
		 */
		void addConfigurationParameter( string configString ) ;

		/*!
		 * @brief
		 *
		 * @return
		 */
		string getConfigurationValue( string configParam );

		/*!
		 * @brief
		 *
		 * @return
		 */
		bool hasConfigurationValue( string configParam, string val );

		/*!
		 * @brief
		 *
		 * @return
		 */
		bool isTrue( string configParam );

		string XML() const;

		void Print() const;

		bool empty() const;

		void AddDaughterPDF( const IPDF* input );

		vector<IPDF*> GetDaughterPDFs() const;

		void SetPhaseSpaceBoundary( PhaseSpaceBoundary* input );

		PhaseSpaceBoundary* GetPhaseSpaceBoundary() const;

		void SetFractionName( const string thisName );

		string GetFractionName() const;

		void SetPDFLabel( const string input );

		string GetPDFLabel() const;

	private:

		//      Uncopyable This Way!
		PDFConfigurator& operator = ( const PDFConfigurator& );

		vector<string> defaultNames;
		vector<string> replacementNames;

		vector<string> configParameters;
		vector<string> configValues;

		vector<IPDF*> daughterPDFs;
		PhaseSpaceBoundary* thisPDFsBoundary;

		string fractionName;

		string PDFLabel;
};

#endif

