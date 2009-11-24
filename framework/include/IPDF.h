/**
        @interface IPDF

        Common interface for all PDFs

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef IPDF_H
#define IPDF_H

#include "DataPoint.h"
#include "PhaseSpaceBoundary.h"
#include "ParameterSet.h"
#include <vector>
#include <string>

using namespace std;

class IPDF
{
	public:
		//Indicate whether the function has been set up correctly
		virtual bool IsValid() = 0;

		//Set the function parameters
		virtual bool SetPhysicsParameters( ParameterSet* ) = 0;

		//Return the integral of the function over the given boundary
		virtual double Integral( DataPoint*, PhaseSpaceBoundary* ) = 0;

		//Return the function value at the given point
		virtual double Evaluate( DataPoint* ) = 0;

		//Return a prototype data point
		virtual vector<string> GetPrototypeDataPoint() = 0;

		//Return a prototype set of physics parameters
		virtual vector<string> GetPrototypeParameterSet() = 0;

		//Return a list of parameters not to be integrated
		virtual vector<string> GetDoNotIntegrateList() = 0;
};

#endif
