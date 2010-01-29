// $Id: Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.cpp,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.cpp
 *
 *  RapidFit PDF for Bd2JpsiPhi with average angular acceptance input as a paramter
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2010-01-12
 */

#include "Bd2JpsiKstar_withTimeRes_withAverageAngAcc.h"
#include "Mathematics.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

//Constructor
Bd2JpsiKstar_withTimeRes_withAverageAngAcc::Bd2JpsiKstar_withTimeRes_withAverageAngAcc() : 
	// Physics parameters
	gammaName     ( "gamma" )    
	, deltaMName    ( "deltaM")
	, delta_paraName( "delta_para" )
	, delta_perpName( "delta_perp" )
	, Azero_sqName  ( "Azero_sq" )
	, Apara_sqName  ( "Apara_sq" )
	, angAccI1Name	( "angAccI1" )
	, angAccI2Name	( "angAccI2" )
	, angAccI3Name	( "angAccI3" )
	, angAccI4Name	( "angAccI4" )
	, angAccI5Name	( "angAccI5" )
	, angAccI6Name	( "angAccI6" )
	, timeRes1Name	( "timeResolution1" )
	, timeRes2Name	( "timeResolution2" )
	, timeRes1FractionName	( "timeResolution1Fraction" )

	// Observables (What we want to gain from the pdf after inserting physics parameter values)
	, timeName	( "time" )
	, cosThetaName	( "cosTheta" )
	, phiName	( "phi" )
	, cosPsiName	( "cosPsi" )
	//, timeres	( "resolution" )
	, normalisationCacheValid(false)
, evaluationCacheValid(false)
{
	MakePrototypes();
}

//Make the data point and parameter set
void Bd2JpsiKstar_withTimeRes_withAverageAngAcc::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( timeName );
	allObservables.push_back( cosThetaName );
	allObservables.push_back( phiName );
	allObservables.push_back( cosPsiName );

	// Need to think about additional parameters like
	// event-by-event propertime resolution and acceptance.
	// This will require event-by-event PDF normalisation,
	// but we are already doing this for tagging.

	//Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( gammaName );
	parameterNames.push_back( Apara_sqName );
	parameterNames.push_back( Azero_sqName );
	parameterNames.push_back( delta_paraName );
	parameterNames.push_back( delta_perpName );
	parameterNames.push_back( deltaMName );
	parameterNames.push_back( timeRes1Name );
	parameterNames.push_back( timeRes2Name );
	parameterNames.push_back( timeRes1FractionName );
	parameterNames.push_back( angAccI1Name );
	parameterNames.push_back( angAccI2Name );
	parameterNames.push_back( angAccI3Name );
	parameterNames.push_back( angAccI4Name );
	parameterNames.push_back( angAccI5Name );
	parameterNames.push_back( angAccI6Name );
	allParameters = *( new ParameterSet(parameterNames) );

	valid = true;
}

//Destructor
Bd2JpsiKstar_withTimeRes_withAverageAngAcc::~Bd2JpsiKstar_withTimeRes_withAverageAngAcc()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bd2JpsiKstar_withTimeRes_withAverageAngAcc::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	normalisationCacheValid = false;
	evaluationCacheValid = false;
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
        // Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
        gamma      = allParameters.GetPhysicsParameter( gammaName )->GetValue();
        deltaMs    = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
        Azero_sq   = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
        Apara_sq   = allParameters.GetPhysicsParameter( Apara_sqName )->GetValue();
		delta_para = allParameters.GetPhysicsParameter( delta_paraName )->GetValue();
		delta_perp = allParameters.GetPhysicsParameter( delta_perpName )->GetValue();
        timeRes1 = allParameters.GetPhysicsParameter( timeRes1Name )->GetValue();
        timeRes2 = allParameters.GetPhysicsParameter( timeRes2Name )->GetValue();
        timeRes1Frac = allParameters.GetPhysicsParameter( timeRes1FractionName )->GetValue();
        angAccI1 = allParameters.GetPhysicsParameter( angAccI1Name )->GetValue();
        angAccI2 = allParameters.GetPhysicsParameter( angAccI2Name )->GetValue();
        angAccI3 = allParameters.GetPhysicsParameter( angAccI3Name )->GetValue();
        angAccI4 = allParameters.GetPhysicsParameter( angAccI4Name )->GetValue();
        angAccI5 = allParameters.GetPhysicsParameter( angAccI5Name )->GetValue();
        angAccI6 = allParameters.GetPhysicsParameter( angAccI6Name )->GetValue();

        Aperp_sq = 1 - Azero_sq - Apara_sq;
		AparaAperp = sqrt(Apara_sq)*sqrt(Aperp_sq);
        AzeroApara = sqrt(Azero_sq)*sqrt(Apara_sq);
        AzeroAperp = sqrt(Azero_sq)*sqrt(Aperp_sq);

	return isOK;
}

//Return a list of parameters not to be integrated
vector<string> Bd2JpsiKstar_withTimeRes_withAverageAngAcc::GetDoNotIntegrateList()
{
	vector<string> doNotIntList;
	return doNotIntList;
}

//Calculate the function value
double Bd2JpsiKstar_withTimeRes_withAverageAngAcc::Evaluate(DataPoint * measurement)
{
	time = measurement->GetObservable( timeName )->GetValue();	
        cosTheta = measurement->GetObservable( cosThetaName )->GetValue();
        phi      = measurement->GetObservable( phiName )->GetValue();
        cosPsi   = measurement->GetObservable( cosPsiName )->GetValue();

	//cout << gamma << " " << Aperp_sq << " " << Azero_sq << endl;
	
        if(timeRes1Frac >= 0.9999)
	{
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                return buildPDFnumerator();
        }
        else
	{
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                double val1 = buildPDFnumerator();
                // Set the member variable for time resolution to the second value and calculate
                timeRes = timeRes2;
                double val2 = buildPDFnumerator();
                return timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
        }
}

double Bd2JpsiKstar_withTimeRes_withAverageAngAcc::buildPDFnumerator()
{
	// The angular functions f1->f6 as defined in roadmap Table 1.(same for Kstar)
	double f1, f2, f3, f4, f5, f6;
	Mathematics::getBs2JpsiPhiAngularFunctions( f1, f2, f3, f4, f5, f6, cosTheta, phi, cosPsi );	
	
	// The time dependent amplitudes as defined in roadmap Eqns 48 -> 59  //No tagging so only need 2 (h± pg 72)
	// First for the B
	double AzeroAzeroB, AparaAparaB, AperpAperpB;
	double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB;
	
	getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB
			, ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
							   );
	
	//W+
	double v1 = f1 * AzeroAzeroB
		+ f2 * AparaAparaB
		+ f3 * AperpAperpB
		+ f4 * ImAparaAperpB
		+ f5 * ReAzeroAparaB
		+ f6 * ImAzeroAperpB; 

	/*
	if (isnan(v1))
	{
		cout << f1 << " " << f2 << " " << f3 << " " <<f4 << " " << f5 << " " << f6 << endl;
			cout << AzeroAzeroB << " " << AperpAperpB << " " << AparaAparaB << " " << ImAparaAperpB << " " << ReAzeroAparaB << " " << ImAzeroAperpB << endl;
	}
	 */
		return v1;
}


double Bd2JpsiKstar_withTimeRes_withAverageAngAcc::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
		
	IConstraint * timeBound = boundary->GetConstraint("time");
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		return -1.;
	}
	else
	{
		tlow = timeBound->GetMinimum();
		thigh = timeBound->GetMaximum();
	}


        if(timeRes1Frac >= 0.9999)
        {
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                return buildPDFdenominator();
        }
        else
        {
                // Set the member variable for time resolution to the first value and calculate
                timeRes = timeRes1;
                double val1 = buildPDFdenominator();
                // Set the member variable for time resolution to the second value and calculate
                timeRes = timeRes2;
                double val2 = buildPDFdenominator();
                return timeRes1Frac*val1 + (1. - timeRes1Frac)*val2;
        }
}
	
double Bd2JpsiKstar_withTimeRes_withAverageAngAcc::buildPDFdenominator()
{

	if (!normalisationCacheValid)
	{
		// The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
		getTimeAmplitudeIntegrals(  cachedAzeroAzeroIntB
				, cachedAparaAparaIntB
				, cachedAperpAperpIntB
				, cachedAparaAperpIntB
				, cachedAzeroAparaIntB
								  , cachedAzeroAperpIntB);

		
		normalisationCacheValid = true;
	}

	double v1 = cachedAzeroAzeroIntB * angAccI1
		+ cachedAparaAparaIntB * angAccI2
		+ cachedAperpAperpIntB * angAccI3 
		+ cachedAparaAperpIntB * angAccI4 
		+ cachedAzeroAparaIntB * angAccI5 
		+ cachedAzeroAperpIntB * angAccI6; 

	return v1;
}

void Bd2JpsiKstar_withTimeRes_withAverageAngAcc::getTimeDependentAmplitudes(  double & AzeroAzero
		, double & AparaApara
		, double & AperpAperp
		, double & ImAparaAperp
		, double & ReAzeroApara	
		, double & ImAzeroAperp		)
{
	// Quantities depending only on physics parameters can be cached
	if ( !evaluationCacheValid )
	{
		cachedAzero = sqrt( Azero_sq );
		cachedApara = sqrt( Apara_sq );
		cachedAperp = sqrt( Aperp_sq );
		
		cachedSinDeltaPerpPara	= sin( delta_perp - delta_para );
		cachedCosDeltaPara	= cos( delta_para );
		cachedSinDeltaPerp	= sin( delta_perp );

		
		evaluationCacheValid = true;
	}

      
	//cout << gamma << " " << deltaGamma << " " << Azero_sq << " " << Aperp_sq << endl;
	//cout << cachedExpCosh << " " << cachedExpSinh << " " << cachedExpCos << " " << cachedExpSin << endl;
	// Now calculate the amplitudes
	
	double Exp = Mathematics::Exp(time, gamma, timeRes);
	
		AzeroAzero = Azero_sq * Exp;  // changed- see note 2009-015 eq 11-13
        AparaApara = Apara_sq * Exp;  //
        AperpAperp = Aperp_sq * Exp;  //

        ImAparaAperp = cachedApara*cachedAperp * cachedSinDeltaPerpPara * Exp;    //See http://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=33933 page14

		ReAzeroApara = cachedAzero*cachedApara * cachedCosDeltaPara * Exp;  
	
        ImAzeroAperp = cachedAzero*cachedAperp * cachedSinDeltaPerp * Exp;
	
	//if ( isnan(ImAparaAperp)) cout << Azero_sq << " " << Apara_sq << " " << Aperp_sq << " " << Exp << endl;
	
	return;
}

void Bd2JpsiKstar_withTimeRes_withAverageAngAcc::getTimeAmplitudeIntegrals( double & AzeroAzeroInt
		, double & AparaAparaInt
		, double & AperpAperpInt
		, double & AparaAperpInt
		, double & AzeroAparaInt
		, double & AzeroAperpInt	)
{

	double ExpInt = Mathematics::ExpInt(tlow, thigh, gamma, timeRes);
	
		AzeroAzeroInt = Azero_sq * ExpInt;     
        AparaAparaInt = Apara_sq  * ExpInt; 
        AperpAperpInt = Aperp_sq * ExpInt;  
	
		AparaAperpInt = AparaAperp * cachedSinDeltaPerpPara * ExpInt;    
															
		AzeroAparaInt = AzeroApara * cachedCosDeltaPara * ExpInt;	
															
		AzeroAperpInt = AzeroAperp * cachedSinDeltaPerp * ExpInt;	
															
															
	return;
}