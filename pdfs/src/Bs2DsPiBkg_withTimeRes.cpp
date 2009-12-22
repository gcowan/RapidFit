// $Id: Bs2DsPiBkg_withTimeRes.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class Bs2DsPiBkg_withTimeRes Bs2DsPiBkg_withTimeRes.cpp
 *
 *  PDF for Bs2DsPi  background from Bd->Dpi
 *
 *  @author Greig A Cowan, Peter Clarke
 *  @date 2009-11-13
 */

#include "Bs2DsPiBkg_withTimeRes.h"
#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"
#include "Mathematics.h"

//Constructor
Bs2DsPiBkg_withTimeRes::Bs2DsPiBkg_withTimeRes() : 
	// Observables
	  timeName	( "time" )

	// Physics parameters
	, lifetimeBdName	( "lifetimeBd" )
	, timeResName	( "timeRes_DsPi_background" )
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bs2DsPiBkg_withTimeRes::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( lifetimeBdName );
	parameterNames.push_back( timeResName );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bs2DsPiBkg_withTimeRes::~Bs2DsPiBkg_withTimeRes()
{
}

//..................................
//Calculate the PDF value

double Bs2DsPiBkg_withTimeRes::Evaluate(DataPoint * measurement)
{
	// Observable
	double time = measurement->GetObservable( timeName )->GetValue();

	// Parameters
  	double timeRes = allParameters.GetPhysicsParameter( timeResName )->GetValue();
  	double lifetimeBd = allParameters.GetPhysicsParameter( lifetimeBdName )->GetValue();
  
	double gammaBd = 1. / lifetimeBd ;
	
	double val           = Mathematics::Exp( time, gammaBd, timeRes ) ;

	return val;
}


//.................................
// Calculate the PDF normalisation 
double Bs2DsPiBkg_withTimeRes::Normalisation(PhaseSpaceBoundary * boundary)

{
	double tmin = 0.;
	double tmax = 0.;
	IConstraint * timeBound = boundary->GetConstraint("time");
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
			cerr << "Bound on time not provided" << endl;
			return -1.;
	}
	else
	{
			tmin = timeBound->GetMinimum();
			tmax = timeBound->GetMaximum();
	}
	
	// Parameters
  	double timeRes = allParameters.GetPhysicsParameter( timeResName )->GetValue();
  	double lifetimeBd = allParameters.GetPhysicsParameter( lifetimeBdName )->GetValue();
	
	double gammaBd = 1. / lifetimeBd ;
	
	double val           = Mathematics::ExpInt( tmin, tmax, gammaBd, timeRes ) ;
	
 
	return val;
}

