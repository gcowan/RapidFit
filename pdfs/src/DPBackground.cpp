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
#ifdef __RAPIDFIT_USE_GSL
#include <gsl/gsl_sf_legendre.h>
#endif

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
	//this->TurnCachingOff();
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

  c[0][0][0][0] = 0.070524;// +- 0.000000
  c[0][0][0][2] = -0.006140;// +- 0.001129
  c[0][0][1][2] = 0.005151;// +- 0.001165
  c[1][1][0][0] = 0.013604;// +- 0.003015
  c[1][1][1][2] = -0.009820;// +- 0.003020
  c[2][0][0][0] = -0.047773;// +- 0.002055
  c[2][0][0][2] = 0.011558;// +- 0.002108
  c[2][1][0][0] = -0.019409;// +- 0.003751
  c[2][2][0][0] = 0.021200;// +- 0.004837
  c[3][0][0][0] = 0.015202;// +- 0.002853
  c[3][1][0][0] = -0.017231;// +- 0.004905
  c[4][0][0][0] = -0.024147;// +- 0.003050
  c[4][1][0][0] = 0.029003;// +- 0.005208
  c[5][0][0][0] = -0.024514;// +- 0.003292
  c[5][4][2][2] = 0.030472;// +- 0.010019
  c[6][0][0][0] = 0.029808;// +- 0.003891
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
	//this->TurnCachingOff();

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

#ifdef __RAPIDFIT_USE_GSL

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
#endif
    return 0 ;
}

double DPBackground::Normalisation(PhaseSpaceBoundary * boundary)
{
        (void) boundary;
	return -1.;
}

