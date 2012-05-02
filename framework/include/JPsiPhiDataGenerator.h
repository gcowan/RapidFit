/*!
 * @class JPsiPhiDataGenerator
 *
 * A data generator implementing preselection for JPsiPhi compatible events
 *
 * @author Benjamin M Wynne bwynne@cern.ch
*/

#pragma once
#ifndef JPSIPHI_DATA_GENERATOR_H
#define JPSIPHI_DATA_GENERATOR_H

//	RapidFit Headers
#include "AcceptReject.h"

class JPsiPhiDataGenerator : public AcceptReject
{
	public:
		JPsiPhiDataGenerator( PhaseSpaceBoundary*, IPDF* );
		~JPsiPhiDataGenerator();

	protected:
		virtual bool Preselection( DataPoint*, double );
};

#endif

