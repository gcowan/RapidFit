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
#include "SlicedAcceptance.h"
//	System Headers
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

using namespace::std;

//=======================================

class TimeAccRes : public IResolutionModel
{
	public:

		TimeAccRes( PDFConfigurator* configurator, bool quiet=false ) ;
		TimeAccRes( const TimeAccRes& );
		~TimeAccRes();

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

		bool isPerEvent();

		bool CacheValid() const;

	protected:

		unsigned int numComponents() { return 0; };
		void requestComponent( unsigned int input ) { (void) input; };

		double GetFraction( unsigned int input ) { (void) input; return 0.; };

	private:
		double GetThisScale() { return 0.; };

		IResolutionModel* resolutionModel;
		SlicedAcceptance* timeAcc;

		PDFConfigurator* _config;

		void ConfigTimeAcc( PDFConfigurator* configurator, bool quiet );
		void ConfigTimeRes( PDFConfigurator* configurator, bool quiet );
};

#endif

