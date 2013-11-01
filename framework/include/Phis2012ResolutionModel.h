/**
  @class ResolutionModel

  A class for implementing a decay time resolution model

  This one is a single Gaussian event-by-event model, with a single scale factor

  @author Dianne Ferguson
  @data 2013-10-17
  */

#pragma once
#ifndef Phis2012Resolution_Model_H
#define Phis2012Resolution_Model_H

//	RapidFir Headers
#include "ParameterSet.h"
#include "DataPoint.h"
#include "PDFConfigurator.h"
#include "Observable.h"
#include "IResolutionModel.h"
//	System Headers
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

using namespace::std;

//=======================================

class Phis2012ResolutionModel : public IResolutionModel
{
	public:

		Phis2012ResolutionModel( PDFConfigurator* configurator, bool quiet=false ) ;

		void addParameters( vector<string> & parameterNames ) ;
		void setParameters( ParameterSet & parameters ) ;

		void addObservables( vector<string> & observableNames ) ;
		void setObservables( DataPoint * measurement ) ;

		double Exp( double time, double gamma ) ;
		double ExpInt( double tlow, double thigh, double gamma ) ;

		double ExpSin( double time, double gamma, double dms ) ;
		double ExpSinInt( double tlow, double thigh, double gamma, double dms ) ;

		double ExpCos( double time, double gamma, double dms ) ;
		double ExpCosInt( double tlow, double thigh, double gamma, double dms ) ;

		bool isPerEvent() ;

	protected:

		unsigned int numComponents();
		void requestComponent( unsigned int );

		double GetFraction( unsigned int );

	private:

		double GetThisScale();

		ObservableRef sfBarSlopeName;			// Scale to multiply e-by-e resolution
		ObservableRef sfBarOffsetName;			// Scale to multiply e-by-e resolution
		ObservableRef sfSigmaSlopeName;			// Scale to multiply e-by-e resolution
		ObservableRef sfSigmaOffsetName;			// Scale to multiply e-by-e resolution
		ObservableRef sigmaBarName;			// Scale to multiply e-by-e resolution
		double sfBarSlope;
		double sfBarOffset;
		double sfSigmaSlope;
		double sfSigmaOffset;
		double sfbar;
		double sfsigma;
		double sigmaBar;
		double resScale;
		double resScale2;

		ObservableRef eventResolutionName;  // Event-by-event resolution observable
		double eventResolution;

		unsigned int numberComponents;
		unsigned int wantedComponent;

		ObservableRef timeResFracName;
		double resFrac;
};

#endif

