/**
  @class IResolutionModel

  Interface class for Resolution Models

  @author Pete Clarke
  @data 2013-06-03
  */

#pragma once
#ifndef _I_Resolution_Model_H
#define _I_Resolution_Model_H

//	RapidFir Headers
#include "ParameterSet.h"

#include "PDFConfigurator.h"
#include "Observable.h"
//	System Headers
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

using namespace::std;

//=======================================

class IResolutionModel
{
	public:

		virtual void addParameters( vector<string> & parameterNames ) = 0;
		virtual void setParameters( ParameterSet & parameters ) = 0;

		virtual void addObservables( vector<string> & observableNames ) = 0;
		virtual void setObservables( DataPoint * measurement ) = 0;

		virtual double Exp( double time, double gamma ) = 0;
		virtual double ExpInt( double tlow, double thigh, double gamma ) = 0;

		virtual double ExpSin( double time, double gamma, double dms ) = 0;
		virtual double ExpSinInt( double tlow, double thigh, double gamma, double dms ) = 0;

		virtual double ExpCos( double time, double gamma, double dms ) = 0;
		virtual double ExpCosInt( double tlow, double thigh, double gamma, double dms ) = 0;

		virtual bool isPerEvent() = 0;

		virtual ~IResolutionModel() {};
	protected:

		virtual unsigned int numComponents() = 0;
		virtual void requestComponent( unsigned int ) = 0;

		virtual double GetFraction( unsigned int ) = 0;

		IResolutionModel() {};
};

#endif

