#ifndef FOAM_INTEGRATOR_H
#define FOAM_INTEGRATOR_H

//	RapidFit Headers
#include "IPDF.h"
#include "IDataSet.h"
#include "MakeFoam.h"
//	System Headers
#include <vector>

class FoamIntegrator
{
	public:
		FoamIntegrator();
		FoamIntegrator( IPDF*, IDataSet* );
		~FoamIntegrator();

		double Integral( DataPoint*, PhaseSpaceBoundary* );

	private:
		vector<MakeFoam*> allIntegrators;
		vector<string> discreteNames;
		vector< vector<double> > discreteValues;
};

#endif
