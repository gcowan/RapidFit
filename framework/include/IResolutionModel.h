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

class IResolutionModel;
//=======================================
/*!
 * @brief typedef for the class-factory objects which actually create the new class instances in memory
 */
typedef IResolutionModel* CreateResModel_t( PDFConfigurator*, bool );

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

		virtual pair<double,double> ExpCosSin( double time, double gamma, double dms ) { return make_pair(0.,0.); };
		virtual pair<double,double> ExpCosSinInt( double tlow, double thigh, double gamma, double dms ) { return make_pair(0.,0.); };

		virtual bool isPerEvent() = 0;

		virtual ~IResolutionModel() {};

		virtual bool CacheValid() const = 0;

	protected:

		virtual unsigned int numComponents() = 0;
		virtual void requestComponent( unsigned int ) = 0;

		virtual double GetFraction( unsigned int ) = 0;

		IResolutionModel() {};
};

#ifndef __CINT__

/*!
 * @brief Macro for adding a class lookup and copy function instance to the main function index
 *
 * This allows for any standard Configuration functions inherited from each PDF that must be called after the object has been initialized to be called here
 *
 * This allows for objects to be correctly configured without the PDF developer having to care about initializing the PDFs
 * @return This Returns the PDF that has been constructed
 */
#define RESMODEL_CREATOR( X ) \
        extern "C" IResolutionModel* CreateResModel_##X( PDFConfigurator* config, bool quiet ) { \
                IResolutionModel* thisObject = (IResolutionModel*) new X( config, quiet ); \
                return thisObject; \
}

#endif

#endif

