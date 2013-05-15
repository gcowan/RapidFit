// $Id: DPBackground.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class DPBackground DPBackground.cpp
 *
 *  PDF for Bd2JpsiKpi spin-0
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#include "DPBackground.h"

#include <iostream>
#include <gsl/gsl_sf_legendre.h>

PDF_CREATOR( DPBackground );

//Constructor
DPBackground::DPBackground( PDFConfigurator* configurator) :
    // Observables
	  m23Name	( configurator->getName("m23") )
	, cosTheta1Name	( configurator->getName("cosTheta1") )
	, cosTheta2Name	( configurator->getName("cosTheta2") )
	, phiName	( configurator->getName("phi") )

   	, m23(), cosTheta1(), cosTheta2(), phi()
{
	MakePrototypes();

	this->SetNumericalNormalisation( true );
	this->TurnCachingOff();
    for ( int l = 0; l < l_max + 1; l++ )
    {
        for ( int i = 0; i < i_max + 1; i++ )
        {
            for ( int k = 0; k < k_max + 1; k++ )
            {
                for ( int j = 0; j < j_max + 1; j++ )
                {
                    c[l][i][k][j] = 0.;
                }
            }
        }
    }
// These are terms > 3sigma significance
c[0][0][0][0] = 0.070524;// +- 0.000000
c[0][0][0][2] = -0.005394;// +- 0.001064
c[0][0][1][2] = 0.005547;// +- 0.001096
c[0][1][0][0] = -0.006344;// +- 0.001872
c[1][1][0][0] = 0.013542;// +- 0.002817
c[1][1][1][2] = -0.009051;// +- 0.002820
c[2][0][0][0] = -0.048220;// +- 0.001889
c[2][0][0][2] = 0.010622;// +- 0.001962
c[2][1][0][0] = -0.017056;// +- 0.003475
c[2][2][0][0] = 0.016528;// +- 0.004483
c[3][0][0][0] = 0.019681;// +- 0.002697
c[3][1][0][0] = -0.016729;// +- 0.004664
c[4][0][0][0] = -0.029545;// +- 0.002847
c[4][1][0][0] = 0.030923;// +- 0.004890
c[5][0][0][0] = -0.026361;// +- 0.003057
c[5][4][2][2] = 0.028289;// +- 0.009332
c[6][0][0][0] = 0.039573;// +- 0.003638
}

DPBackground::DPBackground( const DPBackground &copy ) :
	BasePDF( (BasePDF)copy )
	,m23Name(copy.m23Name)
	,cosTheta1Name(copy.cosTheta1Name)
	,cosTheta2Name(copy.cosTheta2Name)
	,phiName(copy.phiName)
    ,m23(copy.m23)
	,cosTheta1(copy.cosTheta1)
	,cosTheta2(copy.cosTheta2)
	,phi(copy.phi)
{
	this->SetNumericalNormalisation(true);
	this->TurnCachingOff();

    for ( int l = 0; l < l_max + 1; l++ )
    {
        for ( int i = 0; i < i_max + 1; i++ )
        {
            for ( int k = 0; k < k_max + 1; k++ )
            {
                for ( int j = k; j < j_max + 1; j++ )
                {
                    c[l][i][k][j] = copy.c[l][i][k][j];
                }
            }
        }
    }
}

//Make the data point and parameter set
void DPBackground::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( m23Name );
	allObservables.push_back( cosTheta1Name );
	allObservables.push_back( cosTheta2Name );
	allObservables.push_back( phiName );

    //Make the parameter set
	vector<string> parameterNames;
	allParameters = ParameterSet(parameterNames);
}

//Destructor
DPBackground::~DPBackground()
{
    /*
	if ( useAngularAcceptance ) {
        for ( int k = 0; k < k_max; k+=1 ) // must have l >= k
        {
            for ( int j = 0; j < j_max; j+=1)
            {
                delete[] c[k][j];
            }
            delete[] c[k];
        }
        delete[] c;
    }
    */
}

bool DPBackground::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	return isOK;
}

//Calculate the function value
double DPBackground::Evaluate(DataPoint * measurement)
{
	// Observables
	m23       = measurement->GetObservable( m23Name )->GetValue();
	cosTheta1 = measurement->GetObservable( cosTheta1Name )->GetValue();
	cosTheta2 = measurement->GetObservable( cosTheta2Name )->GetValue();
	phi       = measurement->GetObservable( phiName )->GetValue();
    double m23_mapped = (m23 - 0.64)/(1.59 - 0.64)*2. + (-1); // should really do this in a generic way

	double returnable_value(0.);
    double Q_l(0.);
    double P_i(0.);
    double Y_jk(0.);
    for ( int l = 0; l < l_max+1; l++ )
    {
        for ( int i = 0; i < i_max+1; i++ )
        {
            for ( int k = 0; k < 3; k++)
            {
                for ( int j = 0; j < 3; j+=2 ) // limiting the loop here to only look at terms we need
                {
                    //cout << "likj" << l << " " << i <<  " " << k << " " << j << endl;
                    if (j < k) continue; // must have l >= k
                    Q_l  = gsl_sf_legendre_Pl     (l,    m23_mapped);
                    P_i  = gsl_sf_legendre_Pl     (i,    cosTheta2);
                    // only consider case where k >= 0
                    // these are the real valued spherical harmonics
                    if ( k == 0 ) Y_jk =           gsl_sf_legendre_sphPlm (j, k, cosTheta1);
                    else          Y_jk = sqrt(2) * gsl_sf_legendre_sphPlm (j, k, cosTheta1) * cos(k*phi);
                    returnable_value += c[l][i][k][j]*(Q_l * P_i * Y_jk);
                }
            }
        }
    }
    if( std::isnan(returnable_value) || returnable_value < 0 ) return 0.;
	else return returnable_value;
}

double DPBackground::Normalisation(PhaseSpaceBoundary * boundary)
{
        (void) boundary;
	return -1.;
}

