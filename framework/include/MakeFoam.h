
#pragma once
#ifndef MAKE_FOAM_H
#define MAKE_FOAM_H

//	RapidFit Headers
#include "IPDF.h"
#include "PhaseSpaceBoundary.h"

class MakeFoam
{
	public:
		MakeFoam( IPDF*, PhaseSpaceBoundary*, DataPoint* );
		MakeFoam( const MakeFoam& );
		~MakeFoam();

		void Debug();
		double Integral();

	private:
		//	Uncopyable!
		MakeFoam& operator = ( const MakeFoam& );

		vector<PhaseSpaceBoundary*> finishedCells;
		vector<DataPoint*> centerPoints;
		vector<double> centerValues, cellIntegrals;
		IPDF * integratePDF;
};

#endif

