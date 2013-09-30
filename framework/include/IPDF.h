/*!
 * @interface IPDF
 *
 * @brief Common interface for all PDFs
 *
 * For more details and the implementation of most of these functions see BasePDF
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef RAPIDFIT_IPDF_H
#define RAPIDFIT_IPDF_H

///	ROOT Headers
#include "TString.h"
#include "TRandom3.h"
///	RapidFit Headers
#include "IPDF_MCCaching.h"
#include "IPDF_Framework.h"
#include "IPDF_NormalisationCaching.h"
#include "ParameterSet.h"
#include "ComponentRef.h"
///	System Headers
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

using namespace::std;

class IPDF : public virtual IPDF_NormalisationCaching, public virtual IPDF_MCCaching, public virtual IPDF_Framework
{
	public:
		/*!
		 * Virtual Destructor
		 */
		virtual ~IPDF()
		{
		};

		/*!
		 * Interface Function:
		 * Externally update the PDFs in the PDF
		 */
		virtual void UpdatePhysicsParameters( ParameterSet* ) = 0;

		/*!
		 * Interface Function:
		 * Return the integral of the function over the given boundary
		 */
		virtual double Integral( DataPoint*, PhaseSpaceBoundary* ) = 0;

		/*!
		 * Interface Function:
		 * Return the function value at the given point
		 */
		virtual double Evaluate( DataPoint* ) = 0;

		/*!
		 * Interface Function:
		 * Return the function value at the given point for generation
		 */
		virtual double EvaluateForNumericGeneration( DataPoint* ) = 0;

		/*!
		 * Interface Function
		 * Return the function value at the given point for use by numeric integral
		 */
		virtual double EvaluateForNumericIntegral( DataPoint* ) = 0;

		virtual double EvaluateTimeOnly( DataPoint* ) = 0;

		/*!
		 * Interface Function:
		 * Get the function parameters
		 */
		virtual ParameterSet* GetPhysicsParameters() = 0;

		/*!
		 * Interface Function:
		 * Return a prototype data point
		 */
		virtual vector<string> GetPrototypeDataPoint() = 0;

		/*!
		 * Interface Function:
		 * Return a prototype set of physics parameters
		 */
		virtual vector<string> GetPrototypeParameterSet() = 0;

		/*!
		 * Interface Function:
		 * Return a list of parameters not to be integrated
		 */
		virtual vector<string> GetDoNotIntegrateList() = 0;

		/*!
		 * Interface Function:
		 * Return a list of PDF components addresses in string format
		 */
		virtual vector<string> PDFComponents() = 0;

		virtual string GetComponentName( ComponentRef* ) = 0;

		/*!
		 * Interface Function:
		 * Return the function value at the given point
		 */
		virtual double EvaluateComponent( DataPoint*, ComponentRef* ) = 0;

	protected:

		/*!
		 * Default Constructor
		 */
		IPDF() {};

		virtual bool SetPhysicsParameters( ParameterSet* Input ) = 0;

		/*!
		 * Interface Function:
		 * Return the Integral over the whole PhaseSpace
		 */
		virtual double Normalisation( PhaseSpaceBoundary* ) = 0;

		/*!
		 * Interface Function:
		 * Return the Integral of this DataPoint in this PhaseSpace
		 */
		virtual double Normalisation( DataPoint*, PhaseSpaceBoundary* ) = 0;
};

#endif

