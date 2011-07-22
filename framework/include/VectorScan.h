//	Class for Constructing Vector paths for the scanning code in RapidFit to scan over

#ifndef VECTORSCAN_H
#define VECTORSCAN_H

//	RapidFit Headers
#include "FitResultVector.h"
#include "MinimiserConfiguration.h"
#include "ScanParam.h"
//	System Headers
#include <vector>

class VectorScan{

	public:

		//	Construct a 'Classic' Scan in 1D
		static FitResultVector* Scan1D( MinimiserConfiguration*, ScanParam* );

		//	Construct a 'Classic' Scan in 2D
		static FitResultVector* Scan2D( MinimiserConfiguration*, pair<ScanParam*, ScanParam*> );		

	private:

		//	Internal Logic for simple fitting with a single straight line vector of points
		static FitResultVector* VectorScanEngine( MinimiserConfiguration*, vector<vector<double> >, vector<string> );

};

#endif

