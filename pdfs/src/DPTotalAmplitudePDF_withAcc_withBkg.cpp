// $Id: DPTotalAmplitudePDF_withAcc_withBkg.cpp,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class DPTotalAmplitudePDF_withAcc_withBkg DPTotalAmplitudePDF_withAcc_withBkg.cpp
 *
 *  PDF for Bd2JpsiKpi spin-0
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-11-13
 */

#include "TMath.h"

#include "DPTotalAmplitudePDF_withAcc_withBkg.h"
#include "DPJpsiKaon.hh"
#include "DPZplusK.hh"
#include "DPHelpers.hh"
#include "CalculateAngles.hh"
#include "DPComponent.hh"

#include <iostream>
#include <cmath>
#include "TComplex.h"
#include "RooMath.h"
#ifdef __RAPIDFIT_USE_GSL
#include <gsl/gsl_sf_legendre.h>
#endif

PDF_CREATOR( DPTotalAmplitudePDF_withAcc_withBkg );

//Constructor
DPTotalAmplitudePDF_withAcc_withBkg::DPTotalAmplitudePDF_withAcc_withBkg( PDFConfigurator* configurator) :

	// Physics parameters
	  fractionName		( configurator->getName("BkgFraction") )
	, magA0ZplusName		( configurator->getName("magA0Zplus") )
	, magApZplusName		( configurator->getName("magApZplus") )
	, magAmZplusName		( configurator->getName("magAmZplus") )
	, phaseA0ZplusName	( configurator->getName("phaseA0Zplus") )
	, phaseApZplusName	( configurator->getName("phaseApZplus") )
	, phaseAmZplusName	( configurator->getName("phaseAmZplus") )

	, magA0Kst892Name	( configurator->getName("magA0Kst892") )
	, magApKst892Name	( configurator->getName("magApKst892") )
	, magAmKst892Name	( configurator->getName("magAmKst892") )
	, phaseA0Kst892Name	( configurator->getName("phaseA0Kst892") )
	, phaseApKst892Name	( configurator->getName("phaseApKst892") )
	, phaseAmKst892Name	( configurator->getName("phaseAmKst892") )

	, magA0Kst1410Name	( configurator->getName("magA0Kst1410") )
	, magApKst1410Name	( configurator->getName("magApKst1410") )
	, magAmKst1410Name	( configurator->getName("magAmKst1410") )
	, phaseA0Kst1410Name	( configurator->getName("phaseA0Kst1410") )
	, phaseApKst1410Name	( configurator->getName("phaseApKst1410") )
	, phaseAmKst1410Name	( configurator->getName("phaseAmKst1410") )

	, magA0Kst1680Name	( configurator->getName("magA0Kst1680") )
	, magApKst1680Name	( configurator->getName("magApKst1680") )
	, magAmKst1680Name	( configurator->getName("magAmKst1680") )
	, phaseA0Kst1680Name	( configurator->getName("phaseA0Kst1680") )
	, phaseApKst1680Name	( configurator->getName("phaseApKst1680") )
	, phaseAmKst1680Name	( configurator->getName("phaseAmKst1680") )

	, magA0K01430Name	( configurator->getName("magA0K01430") )
	, phaseA0K01430Name	( configurator->getName("phaseA0K01430") )

	, magA0K21430Name	( configurator->getName("magA0K21430") )
	, magApK21430Name	( configurator->getName("magApK21430") )
	, magAmK21430Name	( configurator->getName("magAmK21430") )
	, phaseA0K21430Name	( configurator->getName("phaseA0K21430") )
	, phaseApK21430Name	( configurator->getName("phaseApK21430") )
	, phaseAmK21430Name	( configurator->getName("phaseAmK21430") )

	, magA0K31780Name	( configurator->getName("magA0K31780") )
	, magApK31780Name	( configurator->getName("magApK31780") )
	, magAmK31780Name	( configurator->getName("magAmK31780") )
	, phaseA0K31780Name	( configurator->getName("phaseA0K31780") )
	, phaseApK31780Name	( configurator->getName("phaseApK31780") )
	, phaseAmK31780Name	( configurator->getName("phaseAmK31780") )

	, magA0K42045Name	( configurator->getName("magA0K42045") )
	, magApK42045Name	( configurator->getName("magApK42045") )
	, magAmK42045Name	( configurator->getName("magAmK42045") )
	, phaseA0K42045Name	( configurator->getName("phaseA0K42045") )
	, phaseApK42045Name	( configurator->getName("phaseApK42045") )
	, phaseAmK42045Name	( configurator->getName("phaseAmK42045") )

	, magA0K52380Name	( configurator->getName("magA0K52380") )
	, magApK52380Name	( configurator->getName("magApK52380") )
	, magAmK52380Name	( configurator->getName("magAmK52380") )
	, phaseA0K52380Name	( configurator->getName("phaseA0K52380") )
	, phaseApK52380Name	( configurator->getName("phaseApK52380") )
	, phaseAmK52380Name	( configurator->getName("phaseAmK52380") )

	, magA0K800Name		( configurator->getName("magA0K800") )
	, phaseA0K800Name	( configurator->getName("phaseA0K800") )

	, magA0NRName		( configurator->getName("magA0NR") )
	, phaseA0NRName  	( configurator->getName("phaseA0NR") )

	, massZplusName		( configurator->getName("massZplus") )
	, widthZplusName	( configurator->getName("widthZplus") )
	, massKst892Name	( configurator->getName("massKst892") )
	, widthKst892Name	( configurator->getName("widthKst892") )
	, massKst1410Name	( configurator->getName("massKst1410") )
	, widthKst1410Name	( configurator->getName("widthKst1410") )
	, massKst1680Name	( configurator->getName("massKst1680") )
	, widthKst1680Name	( configurator->getName("widthKst1680") )
	, massK01430Name	( configurator->getName("massK01430") )
	, widthK01430Name	( configurator->getName("widthK01430") )
	, massK21430Name	( configurator->getName("massK21430") )
	, widthK21430Name	( configurator->getName("widthK21430") )
	, massK31780Name	( configurator->getName("massK31780") )
	, widthK31780Name	( configurator->getName("widthK31780") )
	, massK42045Name	( configurator->getName("massK42045") )
	, widthK42045Name	( configurator->getName("widthK42045") )
	, massK52380Name	( configurator->getName("massK52380") )
	, widthK52380Name	( configurator->getName("widthK52380") )
	, massK800Name	( configurator->getName("massK800") )
	, widthK800Name	( configurator->getName("widthK800") )
	// Observables
	, m23Name	( configurator->getName("m23") )
	, cosTheta1Name	( configurator->getName("cosTheta1") )
	, cosTheta2Name	( configurator->getName("cosTheta2") )
	, phiName	( configurator->getName("phi") )
	, m13Name	( configurator->getName("m13") )
	, cosZName	( configurator->getName("cosZ") )
	, cosPsiZName	( configurator->getName("cosPsiZ") )
	, phiZName	( configurator->getName("phiZ") )
	, alphaName	( configurator->getName("alpha") )
    , pionIDName( configurator->getName("pionID") )

    , mag_LASSName	( configurator->getName("mag_LASS") )
  	, phase_LASSName	( configurator->getName("phase_LASS") )
	, a_LASSName	( configurator->getName("a_LASS") )
	, r_LASSName	( configurator->getName("r_LASS") )
	// The actual values of the parameters and observables
	, fraction()
    , magA0Zplus(),  magApZplus(),  magAmZplus(),  phaseA0Zplus(),  phaseApZplus(),  phaseAmZplus()
	, magA0Kst892(),  magApKst892(), magAmKst892(), phaseA0Kst892(),  phaseApKst892(),  phaseAmKst892()
	, magA0Kst1410(), magApKst1410(), magAmKst1410(), phaseA0Kst1410(), phaseApKst1410(), phaseAmKst1410()
	, magA0Kst1680(), magApKst1680(), magAmKst1680(), phaseA0Kst1680(), phaseApKst1680(), phaseAmKst1680()
	, magA0K01430(),  				  phaseA0K01430()
	, magA0K21430(), magApK21430(), magAmK21430(), phaseA0K21430(), phaseApK21430(), phaseAmK21430()
	, magA0K31780(), magApK31780(), magAmK31780(), phaseA0K31780(), phaseApK31780(), phaseAmK31780()
	, magA0K42045(), magApK42045(), magAmK42045(), phaseA0K42045(), phaseApK42045(), phaseAmK42045()
	, magA0K52380(), magApK52380(), magAmK52380(), phaseA0K52380(), phaseApK52380(), phaseAmK52380()
	, magA0K800(),  				  phaseA0K800()
	, magA0NR(),    				  phaseA0NR()
	, massZplus(), widthZplus()
	, massKst892(), widthKst892()
	, massKst1410(), widthKst1410()
	, massKst1680(), widthKst1680()
	, massK01430(), widthK01430()
	, massK21430(), widthK21430()
	, massK31780(), widthK31780()
	, massK42045(), widthK42045()
	, massK52380(), widthK52380()
	, massK800(), widthK800()
	, m23(), cosTheta1(), cosTheta2(), phi(), pionID()
	, m13(), cosZ(), cosPsiZ(), phiZ(), alpha()
	//LASS parameters
	, a_LASS(), r_LASS(), mag_LASS(), phase_LASS()
// 	, massPsi(3.096916) // Jpsi
    , massPsi(3.686093) // psi(2S)
    , massB(5.2794) // B0
//    , massB(5.36677) // B_s0
	, pMuPlus(0., 0., 0., 0.), pMuMinus(0., 0., 0., 0.), pPi(0., 0., 0., 0.), pK(0., 0., 0., 0.), pB(0., 0., 0., massB)
	, cosARefs(0.)
{
	MakePrototypes();

	componentIndex = 0;

	// Construct all components we need
	DPComponent * tmp;
	// B0 --> Z+ K-
 	//tmp=new DPZplusK(1,0,massB,4.430,0.100,0.493677,
 	//		 0.13957018, 1.6, 1.6, massPsi, 1, 23); // spin 1 Z, for MC testing

	tmp=new DPZplusK(0,1,massB,1.430,0.0100,0.493677,
			0.13957018, 1.6, 1.6, massPsi, 0, 23); // spin 0 Z for datafit
	ZComponents.push_back(tmp);
	// B0 --> J/psi K*
	tmp=new DPJpsiKaon(0, 1, massB, 0.89594, 0.0487, 0.493677,
			0.13957018, 1.6, 1.6, massPsi, 1);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K*(1410)
	tmp=new DPJpsiKaon(0, 1, massB, 1.414, 0.232, 0.493677,
			0.13957018, 1.6, 1.6, massPsi, 1);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K*(1680)
	tmp=new DPJpsiKaon(0, 1, massB, 1.717, 0.322, 0.493677,
			0.13957018, 1.6, 1.6, massPsi, 1);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K0(1430)
	tmp=new DPJpsiKaon(1, 0, massB, 1.425, 0.270, 0.493677,
			0.13957018, 1.6, 1.6, massPsi, 0);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K2(1430)
	tmp=new DPJpsiKaon(1, 2, massB, 1.4324, 0.109, 0.493677,
			0.13957018, 1.6, 1.6, massPsi, 2);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K3(1780)
	tmp=new DPJpsiKaon(2, 3, massB, 1.4324, 0.109, 0.493677,
			0.13957018, 1.6, 1.6, massPsi, 3);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K4(2045)
	tmp=new DPJpsiKaon(3, 4, massB, 1.4324, 0.109, 0.493677,
			0.13957018, 1.6, 1.6, massPsi, 4);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K5(2380)
	tmp=new DPJpsiKaon(4, 5, massB, 1.4324, 0.109, 0.493677,
			0.13957018, 1.6, 1.6, massPsi, 5);
	KpiComponents.push_back(tmp);
	// B0 --> J/psi K(800)
	tmp=new DPJpsiKaon(1, 0, massB, 0.682, 0.574, 0.493677,
			0.13957018, 1.6, 1.6, massPsi, 0);
	KpiComponents.push_back(tmp);

	// Kpi s-wave using LASS
  	tmp=new DPJpsiKaon(1, 0, massB, 1.425, 0.270, 0.493677,
                     0.13957018, 1.6, 1.6, massPsi, 0,
                     "LASS", 1.94, 1.76);

	KpiComponents.push_back(tmp);

	// Kpi s-wave using Non Resonnant
  	tmp=new DPJpsiKaon(1, 0, massB, 1., 0.270, 0.493677,
                     0.13957018, 1.6, 1.6, massPsi, 0,
                     "NR", 1.94, 1.76);

	KpiComponents.push_back(tmp);

	this->SetNumericalNormalisation( true );
	//this->TurnCachingOff();
    useAngularAcceptance = false;
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
    for ( int l = 0; l < l_max_b + 1; l++ )
    {
        for ( int i = 0; i < i_max_b + 1; i++ )
        {
            for ( int k = 0; k < k_max_b + 1; k++ )
            {
                for ( int j = 0; j < j_max_b + 1; j++ )
                {
                    b[l][i][k][j] = 0.;
                }
            }
        }
    }
	useAngularAcceptance = configurator->isTrue( "UseAngularAcceptance" );

    // From this PDF
    // This is for the JpsiK* analysis, using 70MeV window around K*.
    // Acceptance is flat in Kpi mass.
    /*
    c[0][0][0] = 0.063365;// +- 0.000111
    c[0][0][2] = 0.006549;// +- 0.000168
    c[0][0][4] = 0.000538;// +- 0.000168
    c[0][1][2] = -0.001093;// +- 0.000155
    c[0][2][2] = 0.000539;// +- 0.000140
    c[1][0][0] = -0.014598;// +- 0.000225
    c[1][0][2] = -0.001741;// +- 0.000296
    c[1][1][2] = -0.001824;// +- 0.000234
    c[2][0][0] = -0.017138;// +- 0.000321
    c[2][0][2] = -0.001622;// +- 0.000397
    c[3][0][0] = -0.005164;// +- 0.000370
    c[4][0][0] = -0.003148;// +- 0.000424
    */

    // This is for the psi(2S)Kpi analysis.
    // First dimension is mKpi
c[0][0][0][0] = 0.070537;// +- 0.000042
c[0][0][0][2] = 0.003374;// +- 0.000108
c[0][0][1][2] = -0.000793;// +- 0.000108
c[0][0][2][2] = 0.004150;// +- 0.000104
c[0][1][0][0] = -0.009559;// +- 0.000175
c[0][1][0][2] = -0.001790;// +- 0.000180
c[0][1][1][2] = 0.001620;// +- 0.000180
c[0][1][2][2] = 0.000894;// +- 0.000172
c[0][2][0][0] = -0.014521;// +- 0.000227
c[0][2][1][2] = 0.001962;// +- 0.000233
c[0][2][2][2] = -0.002894;// +- 0.000224
c[0][4][0][0] = -0.002925;// +- 0.000307
c[1][0][0][0] = 0.014735;// +- 0.000218
c[1][0][0][2] = -0.004731;// +- 0.000224
c[1][0][1][2] = 0.001244;// +- 0.000223
c[1][0][2][2] = 0.004222;// +- 0.000215
c[1][1][0][0] = 0.013788;// +- 0.000366
c[1][1][1][2] = 0.008654;// +- 0.000375
c[1][2][0][0] = 0.013075;// +- 0.000476
c[1][2][0][2] = 0.008135;// +- 0.000492
c[1][2][2][2] = -0.003336;// +- 0.000467
c[1][3][1][2] = -0.004594;// +- 0.000575
c[2][0][0][0] = 0.002701;// +- 0.000294
c[2][1][1][2] = 0.004370;// +- 0.000507
c[2][2][0][0] = 0.008449;// +- 0.000640
c[2][2][0][2] = 0.006358;// +- 0.000663
c[2][3][0][0] = -0.004795;// +- 0.000759
c[2][3][1][2] = -0.004152;// +- 0.000771
c[3][1][0][0] = -0.007938;// +- 0.000598
c[3][2][0][0] = -0.006623;// +- 0.000772
c[4][1][0][0] = 0.004372;// +- 0.000681
c[4][3][0][0] = 0.006850;// +- 0.001044
c[5][0][0][0] = -0.003132;// +- 0.000446
c[6][0][0][0] = 0.004619;// +- 0.000485
c[6][1][0][0] = -0.005224;// +- 0.000822

    // Background - with my phi-angle definition or Belle's
b[0][0][0][0] = 0.070524;// +- 0.000000
b[0][0][0][2] = -0.007277;// +- 0.001296
b[0][0][1][2] = 0.005762;// +- 0.001339
b[1][0][0][0] = 0.006331;// +- 0.001999
b[1][1][0][0] = 0.012922;// +- 0.003488
b[2][0][0][0] = -0.047201;// +- 0.002412
b[2][0][0][2] = 0.013230;// +- 0.002452
b[2][1][0][0] = -0.023220;// +- 0.004386
b[2][2][0][0] = 0.018912;// +- 0.005656
b[3][0][0][0] = 0.011316;// +- 0.003258
b[3][1][0][0] = -0.017484;// +- 0.005588
b[4][0][0][0] = -0.020099;// +- 0.003545
b[4][1][0][0] = 0.029596;// +- 0.006075
b[5][0][0][0] = -0.024536;// +- 0.003844
b[6][0][0][0] = 0.019579;// +- 0.004502
}

DPTotalAmplitudePDF_withAcc_withBkg::DPTotalAmplitudePDF_withAcc_withBkg( const DPTotalAmplitudePDF_withAcc_withBkg &copy ) :
	BasePDF( (BasePDF)copy )
	,m23Name(copy.m23Name)
	,cosTheta1Name(copy.cosTheta1Name)
	,cosTheta2Name(copy.cosTheta2Name)
	,phiName(copy.phiName)
	,m13Name(copy.m13Name)
	,cosZName(copy.cosZName)
	,cosPsiZName(copy.cosPsiZName)
	,phiZName(copy.phiZName)
	,alphaName(copy.alphaName)
	,pionIDName(copy.pionIDName)
    ,m23(copy.m23)
	,cosTheta1(copy.cosTheta1)
	,cosTheta2(copy.cosTheta2)
	,phi(copy.phi)
    ,m13(copy.m13)
	,cosZ(copy.cosZ)
	,cosPsiZ(copy.cosPsiZ)
	,phiZ(copy.phiZ)
	,alpha(copy.alpha)
	,pionID(copy.pionID)
	,fractionName(copy.fractionName)
	,magA0ZplusName(copy.magA0ZplusName)
	,magApZplusName(copy.magApZplusName)
	,magAmZplusName(copy.magAmZplusName)
	,phaseA0ZplusName(copy.phaseA0ZplusName)
	,phaseApZplusName(copy.phaseApZplusName)
	,phaseAmZplusName(copy.phaseAmZplusName)

	,magA0Kst892Name(copy.magA0Kst892Name)
	,magApKst892Name(copy.magApKst892Name)
	,magAmKst892Name(copy.magAmKst892Name)
	,phaseA0Kst892Name(copy.phaseA0Kst892Name)
	,phaseApKst892Name(copy.phaseApKst892Name)
	,phaseAmKst892Name(copy.phaseAmKst892Name)

	,magA0Kst1410Name(copy.magA0Kst1410Name)
	,magApKst1410Name(copy.magApKst1410Name)
	,magAmKst1410Name(copy.magAmKst1410Name)
	,phaseA0Kst1410Name(copy.phaseA0Kst1410Name)
	,phaseApKst1410Name(copy.phaseApKst1410Name)
	,phaseAmKst1410Name(copy.phaseAmKst1410Name)

	,magA0Kst1680Name(copy.magA0Kst1680Name)
	,magApKst1680Name(copy.magApKst1680Name)
	,magAmKst1680Name(copy.magAmKst1680Name)
	,phaseA0Kst1680Name(copy.phaseA0Kst1680Name)
	,phaseApKst1680Name(copy.phaseApKst1680Name)
	,phaseAmKst1680Name(copy.phaseAmKst1680Name)

	,magA0K01430Name(copy.magA0K01430Name)
	,phaseA0K01430Name(copy.phaseA0K01430Name)

	,magA0K21430Name(copy.magA0K21430Name)
	,magApK21430Name(copy.magApK21430Name)
	,magAmK21430Name(copy.magAmK21430Name)
	,phaseA0K21430Name(copy.phaseA0K21430Name)
	,phaseApK21430Name(copy.phaseApK21430Name)
	,phaseAmK21430Name(copy.phaseAmK21430Name)

	,magA0K31780Name(copy.magA0K31780Name)
	,magApK31780Name(copy.magApK31780Name)
	,magAmK31780Name(copy.magAmK31780Name)
	,phaseA0K31780Name(copy.phaseA0K31780Name)
	,phaseApK31780Name(copy.phaseApK31780Name)
	,phaseAmK31780Name(copy.phaseAmK31780Name)

	,magA0K42045Name(copy.magA0K42045Name)
	,magApK42045Name(copy.magApK42045Name)
	,magAmK42045Name(copy.magAmK42045Name)
	,phaseA0K42045Name(copy.phaseA0K42045Name)
	,phaseApK42045Name(copy.phaseApK42045Name)
	,phaseAmK42045Name(copy.phaseAmK42045Name)

	,magA0K52380Name(copy.magA0K52380Name)
	,magApK52380Name(copy.magApK52380Name)
	,magAmK52380Name(copy.magAmK52380Name)
	,phaseA0K52380Name(copy.phaseA0K52380Name)
	,phaseApK52380Name(copy.phaseApK52380Name)
	,phaseAmK52380Name(copy.phaseAmK52380Name)

	,magA0K800Name(copy.magA0K800Name)
	,phaseA0K800Name(copy.phaseA0K800Name)

	,magA0NRName(copy.magA0NRName)
	,phaseA0NRName(copy.phaseA0NRName)

	,massZplusName(copy.massZplusName)
	,widthZplusName(copy.widthZplusName)
	,massKst892Name(copy.massKst892Name)
	,widthKst892Name(copy.widthKst892Name)
	,massKst1410Name(copy.massKst1410Name)
	,widthKst1410Name(copy.widthKst1410Name)
	,massKst1680Name(copy.massKst1680Name)
	,widthKst1680Name(copy.widthKst1680Name)
	,massK01430Name(copy.massK01430Name)
	,widthK01430Name(copy.widthK01430Name)
	,massK21430Name(copy.massK21430Name)
	,widthK21430Name(copy.widthK21430Name)
	,massK31780Name(copy.massK31780Name)
	,widthK31780Name(copy.widthK31780Name)
	,massK42045Name(copy.massK42045Name)
	,widthK42045Name(copy.widthK42045Name)
	,massK52380Name(copy.massK52380Name)
	,widthK52380Name(copy.widthK52380Name)
	,massK800Name(copy.massK800Name)
	,widthK800Name(copy.widthK800Name)

  	,mag_LASSName(copy.mag_LASSName)
  	,phase_LASSName(copy.phase_LASSName)
	,a_LASSName(copy.a_LASSName)
	,r_LASSName(copy.r_LASSName)

	,fraction(copy.fraction)
	,magA0Zplus(copy.magA0Zplus)
	,magApZplus(copy.magApZplus)
	,magAmZplus(copy.magAmZplus)
	,phaseA0Zplus(copy.phaseA0Zplus)
	,phaseApZplus(copy.phaseApZplus)
	,phaseAmZplus(copy.phaseAmZplus)

	,magA0Kst892(copy.magA0Kst892)
	,magApKst892(copy.magApKst892)
	,magAmKst892(copy.magAmKst892)
	,phaseA0Kst892(copy.phaseA0Kst892)
	,phaseApKst892(copy.phaseApKst892)
	,phaseAmKst892(copy.phaseAmKst892)

	,magA0Kst1410(copy.magA0Kst1410)
	,magApKst1410(copy.magApKst1410)
	,magAmKst1410(copy.magAmKst1410)
	,phaseA0Kst1410(copy.phaseA0Kst1410)
	,phaseApKst1410(copy.phaseApKst1410)
	,phaseAmKst1410(copy.phaseAmKst1410)

	,magA0Kst1680(copy.magA0Kst1680)
	,magApKst1680(copy.magApKst1680)
	,magAmKst1680(copy.magAmKst1680)
	,phaseA0Kst1680(copy.phaseA0Kst1680)
	,phaseApKst1680(copy.phaseApKst1680)
	,phaseAmKst1680(copy.phaseAmKst1680)

	,magA0K01430(copy.magA0K01430)
	,phaseA0K01430(copy.phaseA0K01430)

	,magA0K21430(copy.magA0K21430)
	,magApK21430(copy.magApK21430)
	,magAmK21430(copy.magAmK21430)
	,phaseA0K21430(copy.phaseA0K21430)
	,phaseApK21430(copy.phaseApK21430)
	,phaseAmK21430(copy.phaseAmK21430)

	,magA0K31780(copy.magA0K31780)
	,magApK31780(copy.magApK31780)
	,magAmK31780(copy.magAmK31780)
	,phaseA0K31780(copy.phaseA0K31780)
	,phaseApK31780(copy.phaseApK31780)
	,phaseAmK31780(copy.phaseAmK31780)

	,magA0K42045(copy.magA0K42045)
	,magApK42045(copy.magApK42045)
	,magAmK42045(copy.magAmK42045)
	,phaseA0K42045(copy.phaseA0K42045)
	,phaseApK42045(copy.phaseApK42045)
	,phaseAmK42045(copy.phaseAmK42045)

	,magA0K52380(copy.magA0K52380)
	,magApK52380(copy.magApK52380)
	,magAmK52380(copy.magAmK52380)
	,phaseA0K52380(copy.phaseA0K52380)
	,phaseApK52380(copy.phaseApK52380)
	,phaseAmK52380(copy.phaseAmK52380)

	,magA0K800(copy.magA0K800)
	,phaseA0K800(copy.phaseA0K800)

	,magA0NR(copy.magA0NR)
	,phaseA0NR(copy.phaseA0NR)


	,massZplus(copy.massZplus)
	,widthZplus(copy.widthZplus)
	,massKst892(copy.massKst892)
	,widthKst892(copy.widthKst892)
	,massKst1410(copy.massKst1410)
	,widthKst1410(copy.widthKst1410)
	,massKst1680(copy.massKst1680)
	,widthKst1680(copy.widthKst1680)
	,massK01430(copy.massK01430)
	,widthK01430(copy.widthK01430)
	,massK21430(copy.massK21430)
	,widthK21430(copy.widthK21430)
	,massK31780(copy.massK31780)
	,widthK31780(copy.widthK31780)
	,massK42045(copy.massK42045)
	,widthK42045(copy.widthK42045)
	,massK52380(copy.massK52380)
	,widthK52380(copy.widthK52380)
	,massK800(copy.massK800)
	,widthK800(copy.widthK800)

	,massPsi(copy.massPsi)
	,massB(copy.massB)

	,pMuPlus(copy.pMuPlus)
	,pMuMinus(copy.pMuMinus)
	,pPi(copy.pPi)
	,pK(copy.pK)
	,KpiComponents()//copy.KpiComponents)
	,ZComponents()//copy.ZComponents)
	,wigner(copy.wigner)
	,useAngularAcceptance(copy.useAngularAcceptance)
	,pB(copy.pB)
	,cosARefs(copy.cosARefs)
	,mag_LASS(copy.mag_LASS)
	,phase_LASS(copy.phase_LASS)
        ,a_LASS(copy.a_LASS)
        ,r_LASS(copy.r_LASS)
{
	this->SetNumericalNormalisation(true);
	//this->TurnCachingOff();
	componentIndex = 0;

	//cout << "Making copy of DPTotalAmplitudePDF_withAcc_withBkg. Acceptance: " << useAngularAcceptance << endl;
    //    std::cout << "In DPTotal copy" << pMuPlus.X() << " " << pMuPlus.Y() << " " << pMuPlus.Z() << std::endl;

	for( unsigned int i=0; i < copy.KpiComponents.size(); ++i )
	{
		KpiComponents.push_back( new DPJpsiKaon( *((DPJpsiKaon*)copy.KpiComponents[i]) ) );
	}

	for( unsigned int i=0; i < copy.ZComponents.size(); ++i )
	{
		ZComponents.push_back( new DPZplusK( *((DPZplusK*)copy.ZComponents[i]) ) );
	}

	if ( useAngularAcceptance )
    {
    for ( int l = 0; l < l_max + 1; l++ )
    {
        for ( int i = 0; i < i_max + 1; i++ )
        {
            for ( int k = 0; k < k_max + 1; k++ )
            {
                for ( int j = k; j < j_max + 1; j++ )
                {
                    c[l][i][k][j] = copy.c[l][i][k][j];
                    //d[i][k][j] = copy.d[i][k][j];
                    //e[i][k][j] = copy.e[i][k][j];
                }
            }
        }
    }
    }
    for ( int l = 0; l < l_max_b + 1; l++ )
    {
        for ( int i = 0; i < i_max_b + 1; i++ )
        {
            for ( int k = 0; k < k_max_b + 1; k++ )
            {
                for ( int j = k; j < j_max_b + 1; j++ )
                {
                    b[l][i][k][j] = copy.b[l][i][k][j];
                }
            }
        }
    }
}


//Make the data point and parameter set
void DPTotalAmplitudePDF_withAcc_withBkg::MakePrototypes()
{
	//Make the DataPoint prototype
	allObservables.push_back( m23Name );
	allObservables.push_back( cosTheta1Name );
	allObservables.push_back( cosTheta2Name );
	allObservables.push_back( phiName );
	//allObservables.push_back( m13Name );
	//allObservables.push_back( cosZName );
	//allObservables.push_back( cosPsiZName );
	//allObservables.push_back( phiZName );
	//allObservables.push_back( alphaName );
    allObservables.push_back( pionIDName );

    //Make the parameter set
	vector<string> parameterNames;
	parameterNames.push_back( fractionName );
	parameterNames.push_back( magA0ZplusName );
	parameterNames.push_back( magApZplusName );
	parameterNames.push_back( magAmZplusName );
	parameterNames.push_back( phaseA0ZplusName );
	parameterNames.push_back( phaseApZplusName );
	parameterNames.push_back( phaseAmZplusName );

	parameterNames.push_back( magA0Kst892Name );
	parameterNames.push_back( magApKst892Name );
	parameterNames.push_back( magAmKst892Name );
	parameterNames.push_back( phaseA0Kst892Name );
	parameterNames.push_back( phaseApKst892Name );
	parameterNames.push_back( phaseAmKst892Name );

	parameterNames.push_back( magA0Kst1410Name );
	parameterNames.push_back( magApKst1410Name );
	parameterNames.push_back( magAmKst1410Name );
	parameterNames.push_back( phaseA0Kst1410Name );
	parameterNames.push_back( phaseApKst1410Name );
	parameterNames.push_back( phaseAmKst1410Name );

	parameterNames.push_back( magA0Kst1680Name );
	parameterNames.push_back( magApKst1680Name );
	parameterNames.push_back( magAmKst1680Name );
	parameterNames.push_back( phaseA0Kst1680Name );
	parameterNames.push_back( phaseApKst1680Name );
	parameterNames.push_back( phaseAmKst1680Name );

	parameterNames.push_back( magA0K01430Name );
	parameterNames.push_back( phaseA0K01430Name );

	parameterNames.push_back( magA0K21430Name );
	parameterNames.push_back( magApK21430Name );
	parameterNames.push_back( magAmK21430Name );
	parameterNames.push_back( phaseA0K21430Name );
	parameterNames.push_back( phaseApK21430Name );
	parameterNames.push_back( phaseAmK21430Name );

	parameterNames.push_back( magA0K31780Name );
	parameterNames.push_back( magApK31780Name );
	parameterNames.push_back( magAmK31780Name );
	parameterNames.push_back( phaseA0K31780Name );
	parameterNames.push_back( phaseApK31780Name );
	parameterNames.push_back( phaseAmK31780Name );

	parameterNames.push_back( magA0K42045Name );
	parameterNames.push_back( magApK42045Name );
	parameterNames.push_back( magAmK42045Name );
	parameterNames.push_back( phaseA0K42045Name );
	parameterNames.push_back( phaseApK42045Name );
	parameterNames.push_back( phaseAmK42045Name );

	parameterNames.push_back( magA0K52380Name );
	parameterNames.push_back( magApK52380Name );
	parameterNames.push_back( magAmK52380Name );
	parameterNames.push_back( phaseA0K52380Name );
	parameterNames.push_back( phaseApK52380Name );
	parameterNames.push_back( phaseAmK52380Name );

	parameterNames.push_back( magA0K800Name );
	parameterNames.push_back( phaseA0K800Name );

	parameterNames.push_back( magA0NRName );
	parameterNames.push_back( phaseA0NRName );

	parameterNames.push_back( massZplusName );
	parameterNames.push_back( widthZplusName );
	parameterNames.push_back( massKst892Name );
	parameterNames.push_back( widthKst892Name );
	parameterNames.push_back( massKst1410Name );
	parameterNames.push_back( widthKst1410Name );
	parameterNames.push_back( massKst1680Name );
	parameterNames.push_back( widthKst1680Name );
	parameterNames.push_back( massK01430Name );
	parameterNames.push_back( widthK01430Name );
	parameterNames.push_back( massK21430Name );
	parameterNames.push_back( widthK21430Name );
	parameterNames.push_back( massK31780Name );
	parameterNames.push_back( widthK31780Name );
	parameterNames.push_back( massK42045Name );
	parameterNames.push_back( widthK42045Name );
	parameterNames.push_back( massK52380Name );
	parameterNames.push_back( widthK52380Name );
	parameterNames.push_back( massK800Name );
	parameterNames.push_back( widthK800Name );

	parameterNames.push_back( mag_LASSName );
	parameterNames.push_back( phase_LASSName );
	parameterNames.push_back( a_LASSName );
	parameterNames.push_back( r_LASSName );


	allParameters = ParameterSet(parameterNames);
}

//Destructor
DPTotalAmplitudePDF_withAcc_withBkg::~DPTotalAmplitudePDF_withAcc_withBkg()
{
	// destroy components
	for (unsigned int i = 0; i < KpiComponents.size(); ++i)
	{
		delete KpiComponents[i];
		KpiComponents[i]=0;
	}
	for (unsigned int i=0;i<ZComponents.size();++i)
	{
		delete ZComponents[i];
		ZComponents[i]=0;
	}
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

vector<string> DPTotalAmplitudePDF_withAcc_withBkg::GetDoNotIntegrateList()
{
    vector<string> list;
    list.push_back(pionIDName);
    return list;
}

bool DPTotalAmplitudePDF_withAcc_withBkg::SetPhysicsParameters( ParameterSet * NewParameterSet )
{
	bool isOK = allParameters.SetPhysicsParameters(NewParameterSet);
	fraction      = allParameters.GetPhysicsParameter( fractionName )->GetValue();
	magA0Zplus    = allParameters.GetPhysicsParameter( magA0ZplusName )->GetValue();
	magApZplus    = allParameters.GetPhysicsParameter( magApZplusName )->GetValue();
	magAmZplus    = allParameters.GetPhysicsParameter( magAmZplusName )->GetValue();
	phaseA0Zplus   = allParameters.GetPhysicsParameter( phaseA0ZplusName )->GetValue();
	phaseApZplus   = allParameters.GetPhysicsParameter( phaseApZplusName )->GetValue();
	phaseAmZplus   = allParameters.GetPhysicsParameter( phaseAmZplusName )->GetValue();

	magA0Kst892    = allParameters.GetPhysicsParameter( magA0Kst892Name )->GetValue();
	magApKst892    = allParameters.GetPhysicsParameter( magApKst892Name )->GetValue();
	magAmKst892        = allParameters.GetPhysicsParameter( magAmKst892Name )->GetValue();
	phaseA0Kst892  = allParameters.GetPhysicsParameter( phaseA0Kst892Name )->GetValue();
	phaseApKst892  = allParameters.GetPhysicsParameter( phaseApKst892Name )->GetValue();
	phaseAmKst892  = allParameters.GetPhysicsParameter( phaseAmKst892Name )->GetValue();

	magA0Kst1410   = allParameters.GetPhysicsParameter( magA0Kst1410Name )->GetValue();
	magApKst1410   = allParameters.GetPhysicsParameter( magApKst1410Name )->GetValue();
	magAmKst1410   = allParameters.GetPhysicsParameter( magAmKst1410Name )->GetValue();
	phaseA0Kst1410 = allParameters.GetPhysicsParameter( phaseA0Kst1410Name )->GetValue();
	phaseApKst1410 = allParameters.GetPhysicsParameter( phaseApKst1410Name )->GetValue();
	phaseAmKst1410 = allParameters.GetPhysicsParameter( phaseAmKst1410Name )->GetValue();

	magA0Kst1680   = allParameters.GetPhysicsParameter( magA0Kst1680Name )->GetValue();
	magApKst1680   = allParameters.GetPhysicsParameter( magApKst1680Name )->GetValue();
	magAmKst1680   = allParameters.GetPhysicsParameter( magAmKst1680Name )->GetValue();
	phaseA0Kst1680 = allParameters.GetPhysicsParameter( phaseA0Kst1680Name )->GetValue();
	phaseApKst1680 = allParameters.GetPhysicsParameter( phaseApKst1680Name )->GetValue();
	phaseAmKst1680 = allParameters.GetPhysicsParameter( phaseAmKst1680Name )->GetValue();

	magA0K01430    = allParameters.GetPhysicsParameter( magA0K01430Name )->GetValue();
	phaseA0K01430  = allParameters.GetPhysicsParameter( phaseA0K01430Name )->GetValue();

	magA0K21430    = allParameters.GetPhysicsParameter( magA0K21430Name )->GetValue();
	magApK21430   = allParameters.GetPhysicsParameter( magApK21430Name )->GetValue();
	magAmK21430   = allParameters.GetPhysicsParameter( magAmK21430Name )->GetValue();
	phaseA0K21430  = allParameters.GetPhysicsParameter( phaseA0K21430Name )->GetValue();
	phaseApK21430 = allParameters.GetPhysicsParameter( phaseApK21430Name )->GetValue();
	phaseAmK21430 = allParameters.GetPhysicsParameter( phaseAmK21430Name )->GetValue();

	magA0K31780    = allParameters.GetPhysicsParameter( magA0K31780Name )->GetValue();
	magApK31780   = allParameters.GetPhysicsParameter( magApK31780Name )->GetValue();
	magAmK31780   = allParameters.GetPhysicsParameter( magAmK31780Name )->GetValue();
	phaseA0K31780  = allParameters.GetPhysicsParameter( phaseA0K31780Name )->GetValue();
	phaseApK31780 = allParameters.GetPhysicsParameter( phaseApK31780Name )->GetValue();
	phaseAmK31780 = allParameters.GetPhysicsParameter( phaseAmK31780Name )->GetValue();

	magA0K42045    = allParameters.GetPhysicsParameter( magA0K42045Name )->GetValue();
	magApK42045   = allParameters.GetPhysicsParameter( magApK42045Name )->GetValue();
	magAmK42045   = allParameters.GetPhysicsParameter( magAmK42045Name )->GetValue();
	phaseA0K42045  = allParameters.GetPhysicsParameter( phaseA0K42045Name )->GetValue();
	phaseApK42045 = allParameters.GetPhysicsParameter( phaseApK42045Name )->GetValue();
	phaseAmK42045 = allParameters.GetPhysicsParameter( phaseAmK42045Name )->GetValue();

	magA0K52380    = allParameters.GetPhysicsParameter( magA0K52380Name )->GetValue();
	magApK52380   = allParameters.GetPhysicsParameter( magApK52380Name )->GetValue();
	magAmK52380   = allParameters.GetPhysicsParameter( magAmK52380Name )->GetValue();
	phaseA0K52380  = allParameters.GetPhysicsParameter( phaseA0K52380Name )->GetValue();
	phaseApK52380 = allParameters.GetPhysicsParameter( phaseApK52380Name )->GetValue();
	phaseAmK52380 = allParameters.GetPhysicsParameter( phaseAmK52380Name )->GetValue();

	magA0K800    = allParameters.GetPhysicsParameter( magA0K800Name )->GetValue();
	phaseA0K800  = allParameters.GetPhysicsParameter( phaseA0K800Name )->GetValue();

	magA0NR    = allParameters.GetPhysicsParameter( magA0NRName )->GetValue();
	phaseA0NR  = allParameters.GetPhysicsParameter( phaseA0NRName )->GetValue();


	massZplus  = allParameters.GetPhysicsParameter( massZplusName )->GetValue();
	widthZplus = allParameters.GetPhysicsParameter( widthZplusName )->GetValue();
	massKst892  = allParameters.GetPhysicsParameter( massKst892Name )->GetValue();
	widthKst892 = allParameters.GetPhysicsParameter( widthKst892Name )->GetValue();
	massKst1410  = allParameters.GetPhysicsParameter( massKst1410Name )->GetValue();
	widthKst1410 = allParameters.GetPhysicsParameter( widthKst1410Name )->GetValue();
	massKst1680  = allParameters.GetPhysicsParameter( massKst1680Name )->GetValue();
	widthKst1680 = allParameters.GetPhysicsParameter( widthKst1680Name )->GetValue();
	massK01430  = allParameters.GetPhysicsParameter( massK01430Name )->GetValue();
	widthK01430 = allParameters.GetPhysicsParameter( widthK01430Name )->GetValue();
	massK21430  = allParameters.GetPhysicsParameter( massK21430Name )->GetValue();
	widthK21430 = allParameters.GetPhysicsParameter( widthK21430Name )->GetValue();
	massK31780  = allParameters.GetPhysicsParameter( massK31780Name )->GetValue();
	widthK31780 = allParameters.GetPhysicsParameter( widthK31780Name )->GetValue();
	massK42045  = allParameters.GetPhysicsParameter( massK42045Name )->GetValue();
	widthK42045 = allParameters.GetPhysicsParameter( widthK42045Name )->GetValue();
	massK52380  = allParameters.GetPhysicsParameter( massK52380Name )->GetValue();
	widthK52380 = allParameters.GetPhysicsParameter( widthK52380Name )->GetValue();
	massK800  = allParameters.GetPhysicsParameter( massK800Name )->GetValue();
	widthK800 = allParameters.GetPhysicsParameter( widthK800Name )->GetValue();

	mag_LASS = allParameters.GetPhysicsParameter( mag_LASSName )->GetValue();
	phase_LASS = allParameters.GetPhysicsParameter( phase_LASSName )->GetValue();
	a_LASS = allParameters.GetPhysicsParameter( a_LASSName )->GetValue();
	r_LASS = allParameters.GetPhysicsParameter( r_LASSName )->GetValue();

	// No checks performed here to ensure that parameters are set correctly
	ZComponents[0]  ->setResonanceParameters( massZplus, widthZplus );
	KpiComponents[0]->setResonanceParameters( massKst892, widthKst892 );
	KpiComponents[1]->setResonanceParameters( massKst1410, widthKst1410 );
	KpiComponents[2]->setResonanceParameters( massKst1680, widthKst1680 );
	KpiComponents[3]->setResonanceParameters( massK01430, widthK01430 );
	KpiComponents[4]->setResonanceParameters( massK21430, widthK21430 );
	KpiComponents[5]->setResonanceParameters( massK31780, widthK31780 );
	KpiComponents[6]->setResonanceParameters( massK42045, widthK42045 );
	KpiComponents[7]->setResonanceParameters( massK52380, widthK52380 );
	KpiComponents[8]->setResonanceParameters( massK800, widthK800 );
	KpiComponents[9]->setResonanceParameters( a_LASS, r_LASS );
	ZComponents[0]  ->setHelicityAmplitudes(magA0Zplus, magApZplus, magAmZplus, phaseA0Zplus, phaseApZplus, phaseAmZplus);
	//ZComponents[0]  ->setHelicityAmplitudes(magA0Zplus, magA0Zplus, magA0Zplus, phaseA0Zplus, phaseA0Zplus, phaseA0Zplus);
	KpiComponents[0]->setHelicityAmplitudes(magA0Kst892,  magApKst892, magAmKst892, phaseA0Kst892, phaseApKst892, phaseAmKst892);
	KpiComponents[1]->setHelicityAmplitudes(magA0Kst1410, magApKst1410, magAmKst1410, phaseA0Kst1410, phaseApKst1410, phaseAmKst1410);
	KpiComponents[2]->setHelicityAmplitudes(magA0Kst1680, magApKst1680, magAmKst1680, phaseA0Kst1680, phaseApKst1680, phaseAmKst1680);
	KpiComponents[3]->setHelicityAmplitudes(magA0K01430, 0., 0., phaseA0K01430, 0., 0.);
	KpiComponents[4]->setHelicityAmplitudes(magA0K21430, magApK21430, magAmK21430, phaseA0K21430, phaseApK21430, phaseAmK21430);
	KpiComponents[5]->setHelicityAmplitudes(magA0K31780, magApK31780, magAmK31780, phaseA0K31780, phaseApK31780, phaseAmK31780);
	KpiComponents[6]->setHelicityAmplitudes(magA0K42045, magApK42045, magAmK42045, phaseA0K42045, phaseApK42045, phaseAmK42045);
	KpiComponents[7]->setHelicityAmplitudes(magA0K52380, magApK52380, magAmK52380, phaseA0K52380, phaseApK52380, phaseAmK52380);
	KpiComponents[8]->setHelicityAmplitudes(magA0K800, 0., 0., phaseA0K800, 0., 0.);
	KpiComponents[9]->setHelicityAmplitudes(mag_LASS, 0., 0., phase_LASS, 0., 0.);
	KpiComponents[10]->setHelicityAmplitudes(magA0NR, 0., 0., phaseA0NR, 0., 0.);

	return isOK;
}

//Calculate the function value
double DPTotalAmplitudePDF_withAcc_withBkg::Evaluate(DataPoint * measurement)
{
	// Observables
	m23       = measurement->GetObservable( m23Name )->GetValue();
	cosTheta1 = measurement->GetObservable( cosTheta1Name )->GetValue();
	cosTheta2 = measurement->GetObservable( cosTheta2Name )->GetValue();
	phi       = measurement->GetObservable( phiName )->GetValue();
	//m13       = measurement->GetObservable( m13Name )->GetValue();
	//cosZ = measurement->GetObservable( cosZName )->GetValue();
	//cosPsiZ = measurement->GetObservable( cosPsiZName )->GetValue();
	//phiZ       = measurement->GetObservable( phiZName )->GetValue();
	//alpha      = measurement->GetObservable( alphaName )->GetValue();
	pionID    = measurement->GetObservable( pionIDName )->GetValue();
    double m23_mapped = (m23 - 0.64)/(1.59 - 0.64)*2. + (-1); // should really do this in a generic way
    //double m23_mapped = (m23 - 0.64)/(1.68 - 0.64)*2. + (-1); // should really do this in a generic way

#ifdef __RAPIDFIT_USE_GSL

	double angularAcc(0.);
    double Q_l(0.);
    double P_i(0.);
    double Y_jk(0.);
	if ( useAngularAcceptance )
	{
        for ( int l = 0; l < l_max+1; l++ )
        {
        for ( int i = 0; i < i_max+1; i++ )
        {
            for ( int k = 0; k < 3; k++ )
            {
                for ( int j = 0; j < 3; j+=2 ) // limiting the loop here to only look at terms we need
                {
                    if (j < k) continue; // must have l >= k
                    Q_l  = gsl_sf_legendre_Pl     (l,    m23_mapped);
                    P_i  = gsl_sf_legendre_Pl     (i,    cosTheta2);
                    // only consider case where k >= 0
                    // these are the real valued spherical harmonics
                    if ( k == 0 ) Y_jk =           gsl_sf_legendre_sphPlm (j, k, cosTheta1);
                    else          Y_jk = sqrt(2) * gsl_sf_legendre_sphPlm (j, k, cosTheta1) * cos(k*phi);
                    angularAcc += c[l][i][k][j]*(Q_l * P_i * Y_jk);
                }
            }
        }
        }
    }
    else {
        angularAcc = 1.;//0.0195;
    }
    //if (angularAcc <= 0.) cout << "angular acc " << angularAcc << " " << m23 << " " << m23_mapped << " " << cosTheta1 << " " << phi << " " << cosTheta2 << endl;

    // Calculate the Z angles
    //DPHelpers::calculateFinalStateMomenta(massB, m23, massPsi,
	//cosTheta1,  cosTheta2, phi, pionID, 0.1056583715, 0.1056583715, 0.13957018, 0.493677,
	//pMuPlus, pMuMinus, pPi, pK);
    DPHelpers::calculateFinalStateMomentaBelle(massB, m23, massPsi,
                            cosTheta1, cosTheta2, phi,
                            0.1056583715, 0.139570, 0.493677,
                            pMuPlus, pMuMinus, pPi, pK);

        double belle_m23(0.);
        double belle_cosKPi(0.);
        double belle_cosPsi(0.);
        double belle_phiKPiPsi(0.);
        double belle_m13(0.);
        double belle_cosZ(0.);
        double belle_cosPsi_Z(0.);
        double belle_phiPsiZ(0.);
        double belle_phiZPsiPsi(0.);
        const TLorentzVector myMuPlus(pMuPlus);
        const TLorentzVector myMuMinus(pMuMinus);
        const TLorentzVector myPi(pPi);
        const TLorentzVector myK(pK);
        DPHelpers::Belle(myMuPlus, myMuMinus, myPi, myK
            , belle_m23
            , belle_cosKPi
            , belle_cosPsi
            , belle_phiKPiPsi
            , belle_m13
            , belle_cosZ
            , belle_cosPsi_Z
            , belle_phiPsiZ
            , belle_phiZPsiPsi
            );

    //cout << "K* ntuple " << m23 << " " << cosTheta2 << " " << cosTheta1 << " " << phi << " " << pionID << endl;
    //cout << "K* calcul " << belle_m23 << " " << belle_cosKPi << " " << belle_cosPsi << " " << belle_phiKPiPsi << endl;
	//cout << "Z  ntuple " << m13 << " " << cosZ << " " << cosPsiZ << " " << phiZ << " " << alpha << endl;
	//cout << "Z  calcul " << belle_m13 << " " << belle_cosZ << " " << belle_cosPsi_Z << " " << belle_phiPsiZ << " " << belle_phiZPsiPsi << endl;

	double result = 0.;
	TComplex tmp(0,0);

	// This deals with the separate Kpi components
	unsigned int lower = (unsigned)(componentIndex - 1);
	unsigned int upper = (unsigned)componentIndex;
	unsigned int lowerZ = 0;
	unsigned int upperZ = 0;

	// And this switchs things to deal with the Z components.
    if ( (unsigned)componentIndex > KpiComponents.size() )
	{
		lower = 0;
		upper = 0;
		lowerZ = (unsigned)(componentIndex - KpiComponents.size() - 1);
                upperZ = (unsigned)(componentIndex - KpiComponents.size());

	}

    // just plot the background
    if ( componentIndex == 13 )
    {
		lower  = 0;
		upper  = 0;
		lowerZ = 0;
		upperZ = 0;
    }

    if ( componentIndex == 100 ){
		lower  = 0;
		upper  = (unsigned)KpiComponents.size();
		lowerZ = 0;
		upperZ = 0;
    }

    // Finally, for the total of all components
	if ( componentIndex == 0 ) {
		lower  = 0;
		upper  = (unsigned)KpiComponents.size();
		lowerZ = 0;
		upperZ = (unsigned)ZComponents.size();
	}

	// Now sum over final state helicities (this is not general code, but
	// knows about internals of components
	for (int twoLambda = -2; twoLambda <= 2; twoLambda += 4) // Sum over +-1
	{
		tmp = TComplex(0,0);
		for (int twoLambdaPsi = -2; twoLambdaPsi <= 2; twoLambdaPsi += 2) // Sum over -1,0,+1
		{
			    for (unsigned int i = lower; i < upper; ++i) // sum over all components
			    {
				    tmp += KpiComponents[i]->amplitude(m23, cosTheta1, cosTheta2, phi,
						twoLambda, twoLambdaPsi);
			    }
			    // Now comes sum over Z+ components and lambdaPsiPrime
			    for (unsigned int i = lowerZ; i < upperZ; ++i)
			    {
                    //tmp += ZComponents[i]->amplitudeProperVars(belle_m13, belle_cosZ, belle_cosPsi_Z, belle_phiPsiZ, pionID, twoLambda, twoLambdaPsi);
                    tmp += TComplex::Exp(0.5*TComplex::I()*twoLambda*belle_phiZPsiPsi)*(ZComponents[i]->amplitudeProperVars(belle_m13, belle_cosZ, belle_cosPsi_Z, belle_phiPsiZ, pionID, twoLambda, twoLambdaPsi));
                    //std::cout << TComplex::Exp(0.5*TComplex::I()*twoLambda*belle_phiZPsiPsi) << " " << TComplex(cos(belle_phiZPsiPsi), 0.5*twoLambda*sin(belle_phiZPsiPsi)) << std::endl;
			    }
		}
		result += tmp.Rho2();
	}


	//momenta are defined on eq 39.20a/b of the 2010 PDG
	const double m1 = 0.493677;    // kaon mass
	const double m2 = 0.139570; // pion mass
	const double MB0= massB; // B0 mass

	double t1 = m23*m23-(m1+m2)*(m1+m2);
	double t2 = m23*m23-(m1-m2)*(m1-m2);

	double t31 = MB0*MB0 - (m23 + massPsi)*(m23 + massPsi);
	double t32 = MB0*MB0 - (m23 - massPsi)*(m23 - massPsi);

	double p1_st = sqrt(t1*t2)/m23/2.;
	double p3    = sqrt(t31*t32)/MB0/2.;

	//std::cout << result << " " << angularAcc << " " << p1_st << " " << p3 << " " << result*p1_st*p3 << std::endl;

    double returnable_value = result * angularAcc * p1_st * p3;

    double background(0.);

    if ( (componentIndex == 0 || componentIndex == 13) && fraction > 0. )
    {
    for ( int l = 0; l < l_max_b+1; l++ )
    {
        for ( int i = 0; i < i_max_b+1; i++ )
        {
            for ( int k = 0; k < k_max_b+1; k++)
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
                    background += b[l][i][k][j]*(Q_l * P_i * Y_jk);
                }
            }
        }
    }
    }
    //cout << background << " " << fraction << " " << returnable_value << endl;
    //returnable_value = background;
    returnable_value = returnable_value + fraction*background;

	if( std::isnan(returnable_value) || returnable_value < 0. ) return 0.;
	else return returnable_value;

    #endif

    return 0.;
}

vector<string> DPTotalAmplitudePDF_withAcc_withBkg::PDFComponents()
{
        vector<string> components_list;
        components_list.push_back( "892" );
        components_list.push_back( "1410" );
        components_list.push_back( "1680" );
        components_list.push_back( "1430" );
        components_list.push_back( "1430_2" );
        components_list.push_back( "1780_3" );
        components_list.push_back( "2045_4" );
        components_list.push_back( "2380_5" );
        components_list.push_back( "800" );
        components_list.push_back( "LASS" );
        components_list.push_back( "NR" );
        components_list.push_back( "Z4430" );
        components_list.push_back( "background" );
        components_list.push_back( "0" );
        components_list.push_back( "1-Z" );
        return components_list;
}

//Calculate the function value
double DPTotalAmplitudePDF_withAcc_withBkg::EvaluateComponent(DataPoint * measurement, ComponentRef* Component)
{
        componentIndex = Component->getComponentNumber();
        if( componentIndex == -1 )
        {
                string ComponentName = Component->getComponentName();
                if( ComponentName.compare( "892" ) == 0 )
                {
                        Component->setComponentNumber( 1 );
                        componentIndex = 1;
                }
                else if( ComponentName.compare( "1410" ) == 0 )
                {
                        Component->setComponentNumber( 2 );
                        componentIndex = 2;
                }
                else if( ComponentName.compare( "1680" ) == 0 )
                {
                        Component->setComponentNumber( 3 );
                        componentIndex = 3;
                }
                else if( ComponentName.compare( "1430" ) == 0 )
                {
                        Component->setComponentNumber( 4 );
                        componentIndex = 4;
                }
                else if( ComponentName.compare( "1430_2" ) == 0 )
                {
                        Component->setComponentNumber( 5 );
                        componentIndex = 5;
                }
                else if( ComponentName.compare( "1780_3" ) == 0 )
                {
                        Component->setComponentNumber( 6 );
                        componentIndex = 6;
                }
                else if( ComponentName.compare( "2045_4" ) == 0 )
                {
                        Component->setComponentNumber( 7 );
                        componentIndex = 7;
                }
                else if( ComponentName.compare( "2380_5" ) == 0 )
                {
                        Component->setComponentNumber( 8 );
                        componentIndex = 8;
                }
                else if( ComponentName.compare( "800" ) == 0 )
                {
                        Component->setComponentNumber( 9 );
                        componentIndex = 9;
                }
                else if( ComponentName.compare( "LASS" ) == 0 )
                {
                        Component->setComponentNumber( 10 );
                        componentIndex = 10;
                }
                else if( ComponentName.compare( "NR" ) == 0 )
                {
                        Component->setComponentNumber( 11 );
                        componentIndex = 11;
                }
                else if( ComponentName.compare( "Z4430" ) == 0 )
                {
                        Component->setComponentNumber( 12 );
                        componentIndex = 12;
                }
                else if( ComponentName.compare( "background" ) == 0 )
                {
                        Component->setComponentNumber( 13 );
                        componentIndex = 13;
                }
                else if( ComponentName.compare( "1-Z" ) == 0 )
                {
                        Component->setComponentNumber( 100 );
                        componentIndex = 100;
                }
                else
                {
                        Component->setComponentNumber( 0 );
                        componentIndex = 0;
                }
        }

        double return_value = this->Evaluate( measurement );

        return return_value;
}

double DPTotalAmplitudePDF_withAcc_withBkg::Normalisation(PhaseSpaceBoundary * boundary)
{
        (void) boundary;
	return -1.;
}

