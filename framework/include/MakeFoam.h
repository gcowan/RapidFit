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

	private:
		vector<PhaseSpaceBoundary> finishedCells;
		vector<DataPoint> centerPoints;
};

#endif
