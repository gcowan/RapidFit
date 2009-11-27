#ifndef MAKE_FOAM_H
#define MAKE_FOAM_H

#include "IPDF.h"
#include "PhaseSpaceBoundary.h"

class MakeFoam
{
	public:
		MakeFoam();
		MakeFoam( IPDF*, PhaseSpaceBoundary*, DataPoint* );
		~MakeFoam();

		void Debug();
		double Integral();

	private:
		vector<PhaseSpaceBoundary> finishedCells;
		vector<DataPoint> centerPoints;
		vector<double> centerValues, cellIntegrals;
		IPDF * integratePDF;
};

#endif
