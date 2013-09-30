
/*!
 * @brief Undocumented!
 *
 * @author Unknown
 */

#pragma once
#ifndef FOAM_INTEGRATOR_H
#define FOAM_INTEGRATOR_H

///	RapidFit Headers
#include "IDataSet.h"
#include "MakeFoam.h"
///	System Headers
#include <string>
#include <vector>

class IPDF;
class IDataSet;
class MakeFoam;

using namespace::std;

class FoamIntegrator
{
	public:
		/*!
		 * Correct Copy Constructor
		 */
		FoamIntegrator( const FoamIntegrator& );

		/*!
		 * Constructor with correct number of Arguments
		 */
		FoamIntegrator( IPDF*, IDataSet* );

		/*!
		 * Destructor
		 */
		~FoamIntegrator();

		/*!
		 * Internal Function for Calculating the Integral
		 * The behaviour of this object needs to be explored/documented!
		 */
		double Integral( DataPoint*, PhaseSpaceBoundary* );

	private:
		/*!
		 * Don't Copy the class this way!
		 */
		FoamIntegrator& operator = ( const FoamIntegrator& );

		vector<MakeFoam*> allIntegrators;		/*!	Undocumented	*/
		vector<string> discreteNames;			/*!	Undocumented	*/
		vector< vector<double> > discreteValues;	/*!	Undocumented	*/
};

#endif

