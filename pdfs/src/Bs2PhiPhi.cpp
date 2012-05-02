/** @class Bs2PhiPhi Bs2PhiPhi.cpp
 *
 *  RapidFit PDF for Bs2PhiPhi
 *
 *  @author Young Min Kim
 *  @date 12 Nov 2009
 */

#include "Bs2PhiPhi.h"
#include <iostream>
#include "math.h"
#include "TMath.h"

PDF_CREATOR( Bs2PhiPhi );

//Constructor
Bs2PhiPhi::Bs2PhiPhi( PDFConfigurator* configurator ) : cach_v1(), cach_v2(), cach_Azero(), cach_Apara(), cach_Aperp(), cach_sinPhis(), cach_cosPhis(), cach_sinDelta_(), cach_cosDelta_(), cach_cosDeltaSubtr(),
      normalisationCacheValid(false)
    , evaluationCacheValid(false)
    // Physics parameters
    , gamma_sName     ( configurator->getName("gamma_s") )
    , gamma_lName     ( configurator->getName("gamma_l") )
    , gamma_hName     ( configurator->getName("gamma_h") )
    , deltaMName      ( configurator->getName("deltaM") )
    , Phi_sName       ( configurator->getName("Phi_s") )
    , Azero_sqName    ( configurator->getName("Azero_sq") )
    , Apara_sqName    ( configurator->getName("Apara_sq") )
    , Aperp_sqName    ( configurator->getName("Aperp_sq") )
    , delta_1Name     ( configurator->getName("delta_1") )
    , delta_2Name     ( configurator->getName("delta_2") )

    // Observables
    , timeName   ( configurator->getName("time") )
    , theta_1Name( configurator->getName("theta_1") )
    , theta_2Name( configurator->getName("theta_2") )
    , phiName    ( configurator->getName("phi") )
    , tagName    ( configurator->getName("tag") )
    , mistagName ( configurator->getName("mistag") )
    , timeconstraintName( configurator->getName("time") )
    //, timeres    ( "resolution" )
{
    MakePrototypes();
    cout << "Making PhiPhi" << endl;
}

//Make the data point and parameter set
void Bs2PhiPhi::MakePrototypes()
{
    //Make the DataPoint prototype
    allObservables.push_back( timeName );
    allObservables.push_back( theta_1Name );
    allObservables.push_back( theta_2Name );
    allObservables.push_back( phiName );
    allObservables.push_back( tagName );
    allObservables.push_back( mistagName );

    // Need to think about additional parameters like
    // event-by-event propertime resolution and acceptance.
    // This will require event-by-event PDF normalisation,
    // but we are already doing this for tagging.

    //Make the parameter set
    vector<string> parameterNames;
    parameterNames.push_back( gamma_sName );
    parameterNames.push_back( gamma_lName );
    parameterNames.push_back( gamma_hName );
    parameterNames.push_back( Aperp_sqName );
    parameterNames.push_back( Azero_sqName );
    //parameterNames.push_back( Apara_sqName );
    parameterNames.push_back( delta_1Name );
    parameterNames.push_back( delta_2Name );
    parameterNames.push_back( deltaMName );
    parameterNames.push_back( Phi_sName );
    allParameters = ParameterSet(parameterNames);
}

//Destructor
Bs2PhiPhi::~Bs2PhiPhi()
{
}

//Not only set the physics parameters, but indicate that the cache is no longer valid
bool Bs2PhiPhi::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
    normalisationCacheValid = false;
    evaluationCacheValid = false;
    return allParameters.SetPhysicsParameters(NewParameterSet);
}

//Return a list of parameters not to be integrated
vector<string> Bs2PhiPhi::GetDoNotIntegrateList()
{
    vector<string> mistagList;
    mistagList.push_back(mistagName);
    return mistagList;
}

//Calculate the function value
double Bs2PhiPhi::Evaluate(DataPoint * measurement)
{
    // The angular functions f1->f6 as defined in roadmap Table 1.
    double f1, f2, f3, f4, f5, f6;
    getAngularFunctions( f1, f2, f3, f4, f5, f6, measurement );    

    // The time dependent amplitudes as defined in roadmap Eqns 48 -> 59
    // First for the B
    double AzeroAzeroB, AparaAparaB, AperpAperpB;
    double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB;
    getTimeDependentAmplitudes( AzeroAzeroB, AparaAparaB, AperpAperpB
                              , ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB
                              , measurement, 1);

    // Now for the Bbar
    double AzeroAzeroBbar, AparaAparaBbar, AperpAperpBbar;
    double ImAparaAperpBbar, ReAzeroAparaBbar, ImAzeroAperpBbar;
    getTimeDependentAmplitudes( AzeroAzeroBbar, AparaAparaBbar, AperpAperpBbar
                              , ImAparaAperpBbar, ReAzeroAparaBbar, ImAzeroAperpBbar
                              , measurement, -1);

    // Now need to know the tag and the mistag
    int q = (int)measurement->GetObservable( tagName )->GetValue(); //-1, 0 or +1
    double omega = measurement->GetObservable( mistagName )->GetValue();    
      
    double epsilon[3];
    epsilon[0] = omega;
    epsilon[1] = 0.5;
    epsilon[2] = (1.0 - omega);
    double w1  = epsilon[q + 1];
    double w2  = 1.0 - w1;

    //W+
    double v1 = f1 * AzeroAzeroB
              + f2 * AparaAparaB
              + f3 * AperpAperpB
              + f4 * ImAparaAperpB
              + f5 * ReAzeroAparaB
              + f6 * ImAzeroAperpB; 
  
    //W-
    double v2 = f1 * AzeroAzeroBbar
              + f2 * AparaAparaBbar
              + f3 * AperpAperpBbar
              + f4 * ImAparaAperpBbar
              + f5 * ReAzeroAparaBbar
              + f6 * ImAzeroAperpBbar;  
  
    return ( w1*v1 + w2*v2 );
}


double Bs2PhiPhi::Normalisation(DataPoint * measurement, PhaseSpaceBoundary * boundary)
{
    cout << "Normalisation" << endl;
    // Now need to know the tag and the mistag
    int q = (int)measurement->GetObservable( tagName )->GetValue(); //-1, 0 or +1
    double omega = measurement->GetObservable( mistagName )->GetValue();    

    double epsilon[3];
    epsilon[0] = omega;
    epsilon[1] = 0.5;
    epsilon[2] = (1.0 - omega);
    double w1  = epsilon[q + 1];
    double w2  = 1.0 - w1;

    if (!normalisationCacheValid)
    {
        // The integrals of the time dependent amplitudes as defined in roadmap Eqns 48 -> 59
        double AzeroAzeroIntB, AparaAparaIntB, AperpAperpIntB;
        getTimeAmplitudeIntegrals( AzeroAzeroIntB
                                 , AparaAparaIntB
                                 , AperpAperpIntB
                                 , boundary
                                 , 1); //+1 for B

        double AzeroAzeroIntBbar, AparaAparaIntBbar, AperpAperpIntBbar;
        getTimeAmplitudeIntegrals( AzeroAzeroIntBbar
                                 , AparaAparaIntBbar
                                 , AperpAperpIntBbar
                                 , boundary
                                 , -1); //-1 for Bbar
    
        cach_v1 = AzeroAzeroIntB + AparaAparaIntB + AperpAperpIntB;
        cach_v2 = AzeroAzeroIntBbar + AparaAparaIntBbar + AperpAperpIntBbar;
        normalisationCacheValid = true;
    }
    
    return ( w1*cach_v1 + w2*cach_v2 );
}

void Bs2PhiPhi::getAngularFunctions( double & f1
                                         , double & f2
                                         , double & f3
                                         , double & f4
                                         , double & f5
                                         , double & f6
                                         , DataPoint * measurement)
{
    // Observables (the stuff your experiment measures)
    double theta_1 = measurement->GetObservable( theta_1Name )->GetValue();
    double theta_2 = measurement->GetObservable( theta_2Name )->GetValue();
    double phi     = measurement->GetObservable( phiName )->GetValue();
    
    // Intermediate values for theta1,2
    double costheta_1 = cos(theta_1);
    double costheta_2 = cos(theta_2);
    double sintheta_1 = sin(theta_1);
    double sintheta_2 = sin(theta_2);
    
    double Ct1Ct1Ct2Ct2 = costheta_1*costheta_1*costheta_2*costheta_2;           //cos^2(theta_1) * cos^2(theta_2)
    double St1St1St2St2 = sintheta_1*sintheta_1*sintheta_2*sintheta_2;           //sin^2(theta_1) * sin^2(theta_2)
    double S2t1S2t2     = 2.*sintheta_1*costheta_1  *  2.*sintheta_2*costheta_2; //sin(2*theta_1) * sin(2*theta_2), sin2A = 2sinAcosA
    
    // Intermediate values for phi
    double cosPhi  = cos(phi);
    double sinPhi  = sin(phi);
    double cos2Phi = 2.*cosPhi - 1.; //cos(2A) = cos^2(A) - sin^2(A) = 2cos^2(A) - 1
    double sin2Phi = 2.*sinPhi*cosPhi;
    
    double Pi = TMath::Pi();
    double norm = Pi*Pi*Pi;  // Pi^3, the factor you get when integrating f1->6 over all angles (theta1,2 0->2Pi, phi 0->Pi)
                             // Same factor for f1->3 (after removing f1's "4" coefficient) and f4->6 give 0
    
    //Expressions from Equation 6.21 in Nick's Thesis
    f1 = 4.* Ct1Ct1Ct2Ct2 * norm;
    f2 = St1St1St2St2*(1. + cos2Phi) * norm;
    f3 = St1St1St2St2*(1. - cos2Phi) * norm;
    f4 = -2.*St1St1St2St2*sin2Phi * norm;
    f5 = sqrt(2.)*S2t1S2t2*cosPhi * norm;
    f6 = -1.*sqrt(2.)*S2t1S2t2*sinPhi * norm;
    
    return;
}

void Bs2PhiPhi::getTimeDependentAmplitudes(  double & AzeroAzero //(These 6 are output parameters)
                                                 , double & AparaApara
                                                 , double & AperpAperp
                                                 , double & ImAparaAperp
                                                 , double & ReAzeroApara
                                                 , double & ImAzeroAperp
                                                 , DataPoint * measurement
                                                 , int Btype )
{
    // Observable
    double time = measurement->GetObservable( timeName )->GetValue();
        
    // Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
    double gamma_s, gamma_l, gamma_h;
    double deltaMs, Phi_s;
    double Azero_sq, Apara_sq, Aperp_sq;
    double delta_1, delta_2;
    getPhysicsParameters( gamma_s, gamma_l, gamma_h, deltaMs, Phi_s
                        , Azero_sq, Apara_sq, Aperp_sq, delta_1, delta_2);

    // Quantities depending only on physics parameters can be cached
    if ( !evaluationCacheValid )
    {
        cach_Azero = sqrt( Azero_sq );
        cach_Apara = sqrt( Apara_sq );
        cach_Aperp = sqrt( Aperp_sq );
        
        cach_sinPhis = sin( Phi_s );
        cach_cosPhis = cos( Phi_s );
        
        cach_sinDelta_[0]  = sin( delta_1 );
        cach_sinDelta_[1]  = sin( delta_2 );
        cach_cosDelta_[0]  = cos( delta_1 );
        cach_cosDelta_[1]  = cos( delta_2 );
        
        cach_cosDeltaSubtr = cos( delta_2 - delta_1 );
        evaluationCacheValid = true;
    }

    // Quantities depending on time cannot be cached
    double expTG_L = exp( -gamma_l*time );
    double expTG_H = exp( -gamma_h*time );
    double expTG_s = exp( -gamma_s*time );
    
    double sinDeltaMsT = sin( deltaMs*time );
    double cosDeltaMsT = cos( deltaMs*time );
    
    // Now calculate the amplitudes
    // as defined in Equation 6.24 from Nick's Thesis
    
    int Swap = 1; //Sign Swap, for LongExpression[1]'s use
    double LongExpression[2];  //K1, K2 and K5 use the same Long Expression inside square brackets
                               //K3 is the same Long Expression with sign changes for all the Phi terms
                               //All 4 expressions get multiplied by 1/2, which is included in LongExpression[]
                              
    double LongerExpression[2]; //K4 and K6 use the same LongerExpression inside square brackets aside from delta_1,2
    
    for(int c = 0; c < 2; c++) //This for loop is to defend against typo propagation through duplicate code lines
    {
        //LongExpression[0] for K1, K2, K5; [1] for K3
        LongExpression[c] = (
                                (1. + Swap * cach_cosPhis) * expTG_L                     //the Gamma_L term
                              + (1. - Swap * cach_cosPhis) * expTG_H                     //the Gamma_H term
                              + Btype * Swap * 2. * expTG_s * sinDeltaMsT * cach_sinPhis //the Gamma_s term
                            ) / 2.;
        
        //LongerExpression[0] for K4; [1] for K6
        LongerExpression[c] = Btype * expTG_s * (   cach_sinDelta_[c] * cosDeltaMsT 
                                                  - cach_cosDelta_[c] * sinDeltaMsT * cach_cosPhis )
                            - 0.5 * ( expTG_H - expTG_L ) * cach_cosDelta_[c] * cach_sinPhis;
        Swap = -1;
    }
    
    AzeroAzero = Azero_sq * LongExpression[0];
    AparaApara = Apara_sq * LongExpression[0];
    AperpAperp = Aperp_sq * LongExpression[1];
    ImAparaAperp = cach_Apara * cach_Aperp * LongerExpression[0];
    ReAzeroApara = cach_Azero * cach_Apara * cach_cosDeltaSubtr * LongExpression[0];
    ImAzeroAperp = cach_Azero * cach_Aperp * LongerExpression[1];

    return;
}

void Bs2PhiPhi::getTimeAmplitudeIntegrals( double & AzeroAzeroInt
                                               , double & AparaAparaInt
                                               , double & AperpAperpInt
                                               , PhaseSpaceBoundary * boundary
                                               , int Btype)
{
    double tlow = 0.;
    double thigh = 0.;
    IConstraint * timeBound = boundary->GetConstraint( timeconstraintName );
    if ( timeBound->GetUnit() == "NameNotFoundError" )
    {
        cerr << "Bound on time not provided" << endl;
        AzeroAzeroInt = -999.;
        AparaAparaInt = -999.;
        AperpAperpInt = -999.;
        return;
    }
    else
    {
        tlow = timeBound->GetMinimum();
        thigh = timeBound->GetMaximum();
    }
    
    // Physics parameters
    double gamma_s, gamma_l, gamma_h;
    double deltaMs, Phi_s;
    double Azero_sq, Apara_sq, Aperp_sq;
    double delta_1, delta_2;
    getPhysicsParameters( gamma_s, gamma_l, gamma_h, deltaMs, Phi_s
                        , Azero_sq, Apara_sq, Aperp_sq, delta_1, delta_2);
    
    //!! Ok, need to redo this part then
//    int Swap = 1; //Sign Swap, for K3 use -1
    
    AzeroAzeroInt = 0.5 * Azero_sq * getMainIntAnswer(tlow, thigh, gamma_s, gamma_l, gamma_h, deltaMs, Phi_s, Btype,  1);
    AparaAparaInt = 0.5 * Apara_sq * getMainIntAnswer(tlow, thigh, gamma_s, gamma_l, gamma_h, deltaMs, Phi_s, Btype,  1);
    AperpAperpInt = 0.5 * Aperp_sq * getMainIntAnswer(tlow, thigh, gamma_s, gamma_l, gamma_h, deltaMs, Phi_s, Btype, -1);
    /*
    AzeroAzeroInt = getAzeroAzeroInt( tlow, thigh, Azero_sq, tauL, tauH, tauBar, deltaMs, Phi_s, Btype);
    AparaAparaInt = getAparaAparaInt( tlow, thigh, Apara_sq, tauL, tauH, tauBar, deltaMs, Phi_s, Btype);
    AperpAperpInt = getAperpAperpInt( tlow, thigh, Aperp_sq, tauL, tauH, tauBar, deltaMs, Phi_s, Btype);
    // No contribution from interference terms here since they drop out when the angular integration is done.
    */
    return;
}

//inline double Bs2PhiPhi::getAzeroAzeroInt(double tmin, double tmax, 
//                               double k0, double tauL, double tauH, double tauBar, double Dms, double phis, int Btype)
inline double Bs2PhiPhi::getMainIntAnswer(double tmin, double tmax, double gamma_s, double gamma_l, double gamma_h, 
                                                double Dms, double phis, int Btype, int Swap)
{
    double gamma_sDms_sq = ( gamma_s*gamma_s + Dms*Dms );
    double cosphis = cos(phis);
    double sinphis = sin(phis);
    
    double answer[2]; //solutions to the integral with t1 and t2 respectively
    double t[2]; t[0] = tmin; t[1] = tmax;
    
    for (int c = 0; c < 2; c++) //This for loop is to defend against typo propagation through duplicate code lines
    {
        //Note: the - sign present in all 3 terms is accounted for at the end,
        //      as we do ans(t1)-ans(t2) not ans(t2)-ans(t1)
        answer[c] = ( (1 + Swap * cosphis) * exp(-gamma_l * t[c]) / gamma_l
                    + (1 - Swap * cosphis) * exp(-gamma_h * t[c]) / gamma_h
                    + Btype * Swap *  2.0  * exp(-gamma_s * t[c]) * sinphis
                            * (gamma_s * sin(Dms * t[c]) + Dms * cos(Dms * t[c]))
                            / gamma_sDms_sq
                    );
    }
    
    return (answer[0]-answer[1]);
}

void Bs2PhiPhi::getPhysicsParameters( double & gamma_s
                                          , double & gamma_l
                                          , double & gamma_h
                                          , double & deltaM
                                          , double & Phi_s
                                          , double & Azero_sq
                                          , double & Apara_sq
                                          , double & Aperp_sq
                                          , double & delta_1
                                          , double & delta_2)
{
    // Physics parameters (the stuff you want to extract from the physics model by plugging in the experimental measurements)
    gamma_s      = allParameters.GetPhysicsParameter( gamma_sName )->GetValue();
    gamma_l      = allParameters.GetPhysicsParameter( gamma_lName )->GetValue();
    gamma_h      = allParameters.GetPhysicsParameter( gamma_hName )->GetValue();
    deltaM       = allParameters.GetPhysicsParameter( deltaMName )->GetValue();
    Phi_s        = allParameters.GetPhysicsParameter( Phi_sName )->GetValue();
    Azero_sq     = allParameters.GetPhysicsParameter( Azero_sqName )->GetValue();
    //Apara_sq     = allParameters.GetPhysicsParameter( Apara_sqName )->GetValue();
    Aperp_sq     = allParameters.GetPhysicsParameter( Aperp_sqName )->GetValue();
    delta_1 = allParameters.GetPhysicsParameter( delta_1Name )->GetValue();
    delta_2 = allParameters.GetPhysicsParameter( delta_2Name )->GetValue();

    Apara_sq = 1 - Azero_sq - Aperp_sq;

    return;
}
