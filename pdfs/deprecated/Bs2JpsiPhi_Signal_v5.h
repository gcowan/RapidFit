// $Id: Bs2JpsiPhi_Signal_v5.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2JpsiPhi_Signal_v5
 *
 *  Bs2JpsiPhi_SignalAlt series with mistag as observable
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2011-11-05
 */

#ifndef Bs2JpsiPhi_Signal_v5_H
#define Bs2JpsiPhi_Signal_v5_H

#include <exception>

#include "BasePDF.h"
#include "PDFConfigurator.h"
#include "SlicedAcceptance.h"
#include "AngularAcceptance.h"
#include "Mathematics.h"
#include <iostream>
#include <cstdlib>
#include <float.h>
#include <vector>

#include "TH1.h"
#include "TCanvas.h"

#define DOUBLE_TOLERANCE 1E-6

using namespace::std;

//============================================================================================
class Bs2JpsiPhi_Signal_v5 : public BasePDF
{
	public:
		Bs2JpsiPhi_Signal_v5( PDFConfigurator* );
		Bs2JpsiPhi_Signal_v5( const Bs2JpsiPhi_Signal_v5& );
		~Bs2JpsiPhi_Signal_v5();

		// Mandatory RapidFit Methods
		virtual double EvaluateForNumericIntegral(DataPoint*);
		virtual double Evaluate(DataPoint*);
		virtual double EvaluateTimeOnly(DataPoint*);
		virtual bool SetPhysicsParameters(ParameterSet*);
		virtual vector<string> GetDoNotIntegrateList();

		vector<string> PDFComponents();

		double EvaluateComponent( DataPoint* input, ComponentRef* );

	private:
		//Bs2JpsiPhi_Signal_v5& operator=( const Bs2JpsiPhi_Signal_v5& );
		void MakePrototypes();
		double normalisationCacheUntagged;
		void prepareCDS();

		void prepareTimeFac();
		void SetupAngularTerms();

		PseudoObservable _expLObs;
		PseudoObservable _expHObs;
		PseudoObservable _expSinObs;
		PseudoObservable _expCosObs;

		PseudoObservable _intexpLObs;
		PseudoObservable _intexpHObs;
		PseudoObservable _intexpSinObs;
		PseudoObservable _intexpCosObs;

		vector<PseudoObservable> _intexpLObs_vec;
		vector<PseudoObservable> _intexpHObs_vec;
		vector<PseudoObservable> _intexpSinObs_vec;
		vector<PseudoObservable> _intexpCosObs_vec;

		int timeBinNum;

		DataPoint* _datapoint;

		bool DebugFlag_v5;

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

		int componentIndex;

		//PELC For debugging purposes
		//TH1D * histOfPdfValues ;
		//TCanvas * c0 ;
		//mutable int histCounter ;
		//~PELC

		// Parameters
		ObservableRef gammaName;		// gamma
		ObservableRef deltaGammaName;	// delta gamma
		ObservableRef deltaMName;		// delta mass
		ObservableRef Azero_sqName;	// amplitude
		ObservableRef Apara_sqName;	// amplitude
		bool _fitDirectlyForApara;
		ObservableRef Aperp_sqName;	// amplitude
		ObservableRef As_sqName;		// amplitude
		ObservableRef CspName;		// amplitude
		ObservableRef delta_zeroName;	// strong phase, set to 0
		ObservableRef delta_paraName;	// strong phase
		ObservableRef delta_perpName;	// strong phase
		ObservableRef delta_sName;		// strong phase for S-wave
		ObservableRef cosdparName;		//PELC-COSDPAR Special for fitting cosdpar separately


		ObservableRef Phi_sName;		// what we want to measure!
		ObservableRef cosphisName;		// fitting cosphis and sinphis independently
		ObservableRef sinphisName;		// fitting cosphis and sinphis independently
		ObservableRef lambdaName;		// magnitude of lambda

		ObservableRef mistagName;		// mistag fraction  - may be used as observable also
		ObservableRef mistagP1Name;		// mistag calib
		ObservableRef mistagP0Name;		// mistag calib
		ObservableRef mistagSetPointName;// mistag calib
		ObservableRef mistagDeltaP1Name;		// mistag calib
		ObservableRef mistagDeltaP0Name;		// mistag calib
		ObservableRef mistagDeltaSetPointName;// mistag calib

		ObservableRef eventResolutionName;			// Scale to multiply all Gaussians with
		ObservableRef resScaleName;			// Scale to multiply all Gaussians with
		ObservableRef res1Name;				// time resolution narrow
		ObservableRef res2Name;				// time resolution wide
		ObservableRef res3Name;				// time resolution tail
		ObservableRef res2FractionName;		// fraction of wide
		ObservableRef res3FractionName;		// fraction of tail
		ObservableRef timeOffsetName;		// time offset

		ObservableRef angAccI1Name;
		ObservableRef angAccI2Name;
		ObservableRef angAccI3Name;
		ObservableRef angAccI4Name;
		ObservableRef angAccI5Name;
		ObservableRef angAccI6Name;
		ObservableRef angAccI7Name;
		ObservableRef angAccI8Name;
		ObservableRef angAccI9Name;
		ObservableRef angAccI10Name;

		// Observables
		ObservableRef timeName;		// proper time
		ObservableRef tagName;			// B tag
		ObservableRef cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		ObservableRef cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
		ObservableRef phiName;			// azimuthal angle of the mu+ in Jpsi frame
		ObservableRef cthetakName;
		ObservableRef cthetalName;
		ObservableRef phihName;

		// PseudoObservables which allow for some calculations to be cached per event during Runtime
		vector<ObservableRef> angularTermDependencies;

		PseudoObservable A0A0_Obs;
		PseudoObservable APAP_Obs;
		PseudoObservable ATAT_Obs;
		PseudoObservable ASAS_Obs;
		PseudoObservable ImAPAT_Obs;
		PseudoObservable ReA0AP_Obs;
		PseudoObservable ImA0AT_Obs;
		PseudoObservable ReASAP_Obs;
		PseudoObservable ImASAT_Obs;
		PseudoObservable ReASA0_Obs;

		double A0A0_value;
		double APAP_value;
		double ATAT_value;
		double ASAS_value;
		double ImAPAT_value;
		double ReA0AP_value;
		double ImA0AT_value;
		double ReASAP_value;
		double ImASAT_value;
		double ReASA0_value;

		// Measured Event Observables
		double t;
		int tag;
		// Transversity angles
		double ctheta_tr;
		double phi_tr;
		double ctheta_1;
		// Helicity angles
		double ctheta_k;
		double ctheta_l;
		double phi_h;
		bool _useHelicityBasis;


		// Physics Fit Parameters
		double _gamma;
		double dgam;

		double Aperp_sq;
		double Apara_sq;
		double Azero_sq;
		double As_sq;
		double Csp;
		double CachedA1;
		double CachedA2;
		double CachedA3;
		double CachedA4;
		double CachedA5;
		double CachedA6;
		double CachedA7;
		double CachedA8;
		double CachedA9;
		double CachedA10;
		void CacheAmplitudesAndAngles();

		double delta_para;
		double delta_perp;
		double delta_zero;
		double delta_s;
		double delta1;
		double delta2;
		double cosdpar; //PELC-COSDPAR Special for fitting cosdpar separately

		double delta_ms;
		double phi_s;
		double _cosphis;
		double _sinphis;
		double lambda;
		double _CC, _DD, _SS;

		double _mistag;
		double _mistagP1;
		double _mistagP0;
		double _mistagSetPoint;
		double _mistagDeltaP1;
		double _mistagDeltaP0;
		double _mistagDeltaSetPoint;

		double resolution;
		double eventResolution;
		double resolutionScale;
		double resolution1;
		double resolution2;
		double resolution3;
		double resolution2Fraction;
		double resolution3Fraction;
		double timeOffset;
		bool _useEventResolution;
		inline bool useEventResolution() const {return _useEventResolution; }
		inline bool useTimeAcceptance() const { return _useTimeAcceptance; }


		double angAccI1;
		double angAccI2;
		double angAccI3;
		double angAccI4;
		double angAccI5;
		double angAccI6;
		double angAccI7;
		double angAccI8;
		double angAccI9;
		double angAccI10;
		AngularAcceptance * angAcc;
		bool _angAccIgnoreNumerator;

		// Othere things calculated later on the fly
		double tlo, thi;

		// stored time primitives
		mutable double expL_stored;
		mutable double expH_stored;
		mutable double expSin_stored;
		mutable double expCos_stored;
		mutable double intExpL_stored;
		mutable double intExpH_stored;
		mutable double intExpSin_stored;
		mutable double intExpCos_stored;
		void preCalculateTimeFactors();
		void preCalculateTimeIntegrals();

		bool timeIntegralCacheValid;
		vector< vector<double> > storeExpL;
		vector< vector<double> > storeExpH;
		vector< vector<double> > storeExpSin;
		vector< vector<double> > storeExpCos;
		void CacheTimeIntegrals();
		void deCacheTimeIntegrals( unsigned int ires, unsigned int islice );

		//Time acceptance
		SlicedAcceptance* timeAcc;

		//Configurationparameters
		bool _useTimeAcceptance;
		bool _numericIntegralForce;
		bool _numericIntegralTimeOnly;
		bool _useCosAndSin;
		bool _useCosDpar;
		bool _usePunziSigmat;
		bool _usePunziMistag;
		bool allowNegativeAsSq;
		bool _usePlotComponents;
		bool _usePlotAllComponents;
		bool performingComponentProjection;
		double _offsetToGammaForBetaFactor;

		double sin_delta_perp_s;
		double cos_delta_perp_s;
		double sin_delta_zero_s;
		double cos_delta_zero_s;
		double sin_delta_para_s;
		double cos_delta_para_s;

		double sin_delta1;
		double cos_delta1;
		double sin_delta2;
		double cos_delta2;
		double sin_delta_2_1;
		double cos_delta_2_1;

		//....................................
		//Internal helper functions

		double stored_AT;
		inline double AT() const { return stored_AT; }
		/*	if( Aperp_sq <= 0. ) return 0. ;
			else return sqrt(Aperp_sq) ;
		}*/
		double stored_AP;
		inline double AP() const { return stored_AP; }
		/*	if( Apara_sq <= 0. ) return 0. ;
			else return sqrt(Apara_sq) ;
		}*/
		double stored_A0;
		inline double A0() const { return stored_A0; }
		/*	if( Azero_sq <= 0. ) return 0. ;
			else return sqrt(Azero_sq) ;
		}*/
		double stored_AS;
		inline double AS() const { return stored_AS; }
		/*	if( As_sq <= 0. ) return 0. ;
			else return sqrt(As_sq) ;
		}*/
		double stored_ASint;
		inline double ASint() const { return stored_ASint; }
		/*	if( As_sq <= 0. ) return 0. ;
			else return sqrt(As_sq) * Csp ;
		}*/

		//....................
		// Safety for gammas
		double stored_gammal;
		inline double gamma_l() const { return stored_gammal; }
		/*	const double gl = gamma() + ( dgam *0.5 ) ;
			if( gl < 0. ) {
				PDF_THREAD_LOCK
				cerr << " In Bs2JpsiPhi_Signal_v5 : gamma_l() < 0 so setting it to 0.0000001 " << endl ;
				PDF_THREAD_UNLOCK
				return 0.0000001 ;
			}
			else
				return gl ;
		}*/

		double stored_gammah;
		inline double gamma_h() const { return stored_gammah; }
		/*	const double gh = gamma() - ( dgam *0.5 ) ;
			if( gh < 0. ) {
				PDF_THREAD_LOCK
				cerr << " In Bs2JpsiPhi_Signal_v5 : gamma_h() < 0 so setting it to 0.0000001 " << endl ;
				PDF_THREAD_UNLOCK
				return 0.0000001 ;
			}
			else
				return gh ;
		}*/

		inline double gamma() const { return _gamma ; }

		//...............
		//tagging


		inline double q() const { return tag ;}

		inline double mistag() const {
			double returnValue = -1000.;

			if( (fabs(q()) < 0.5) || (fabs(q()) > 1.) ) {
				returnValue = 0.5;
			}
			else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
				//Normal case
				returnValue =  _mistagP0 + _mistagP1*(_mistag - _mistagSetPoint ) ;
				if( returnValue < 0 )  returnValue = 0;
				if( returnValue > 0.5) returnValue = 0.5;
			}
			else if( _mistag < 0.0 ) {
				PDF_THREAD_LOCK
				cout << "Bs2JpsiPhi_Signal_v5::mistag() : _mistag < 0 so set to 0 " << endl;
				PDF_THREAD_UNLOCK
				returnValue = 0;
			}
			else if( _mistag > 0.5 ) {
				PDF_THREAD_LOCK
				cout << "Bs2JpsiPhi_Signal_v5::mistag() : _mistag > 0.5 so set to 0.5 "  << endl;
				PDF_THREAD_UNLOCK
				returnValue = 0.5;
			}
			else {
				PDF_THREAD_LOCK
				cout << "Bs2JpsiPhi_Signal_v5::mistag() : WARNING ******If you got here you dont know what you are doing  "  << endl;
				PDF_THREAD_UNLOCK
				exit(1);
			}
			return returnValue;
		}

		inline double mistagB() const {
			double returnValue = -1000.;

			if( (fabs(q()) < 0.5) || (fabs(q()) > 1.) ) {
				returnValue = 0.5;
			}
			else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
				//Normal case
				returnValue =  _mistagP0+(_mistagDeltaP0*0.5) + (_mistagP1+(_mistagDeltaP1*0.5))*(_mistag - (_mistagSetPoint+(_mistagDeltaSetPoint*0.5)) );
				//if( true ) returnValue =  _mistagP0 + (_mistagP1)*(_mistag - (_mistagSetPoint) ) ;  // to mock up independent P1/P0 for each tag
				if( returnValue < 0 )  returnValue = 0;
				if( returnValue > 0.5) returnValue = 0.5;
			}
			else if( _mistag < 0.0 ) {
				PDF_THREAD_LOCK
				cout << "Bs2JpsiPhi_Signal_v5::mistagB() : _mistag < 0 so deltaMistag set to 0 also " << endl;
				PDF_THREAD_UNLOCK
				returnValue = 0;
			}
			else if( _mistag > 0.5 ) {
				PDF_THREAD_LOCK
				cout << "Bs2JpsiPhi_Signal_v5::mistagB() : _mistag > 0.5 so so deltaMistag set to 0.5 also "  << endl;
				PDF_THREAD_UNLOCK
				returnValue = 0.5;
			}
			else {
				PDF_THREAD_LOCK
				cout << "Bs2JpsiPhi_Signal_v5::mistagB() : WARNING ******If you got here you dont know what you are doing  "  << endl;
				PDF_THREAD_UNLOCK
				exit(1);
			}
			return returnValue;
		}

		inline double mistagBbar() const {
			double returnValue = -1000.;

			if( fabs(q()) < 0.5 ) {
				returnValue = 0.5 ;
			}
			else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
				//Normal case
				returnValue =  _mistagP0-(_mistagDeltaP0*0.5) + (_mistagP1-(_mistagDeltaP1*0.5))*(_mistag - (_mistagSetPoint-(_mistagDeltaSetPoint*0.5)) );
				//if( true ) returnValue =   _mistagDeltaP0 + (_mistagDeltaP1)*(_mistag - (_mistagDeltaSetPoint) ) ;// to mock up independent P1/P0 for each tag
				if( returnValue < 0 )  returnValue = 0;
				if( returnValue > 0.5) returnValue = 0.5;
			}
			else if( _mistag < 0.0 ) {
				PDF_THREAD_LOCK
				cout << "Bs2JpsiPhi_Signal_v5::mistagBbar() : _mistag < 0 so deltaMistag set to 0 also " << endl;
				PDF_THREAD_UNLOCK
				returnValue = 0;
			}
			else if( _mistag > 0.5 ) {
				PDF_THREAD_LOCK
				cout << "Bs2JpsiPhi_Signal_v5::mistagBbar() : _mistag > 0.5 so so deltaMistag set to 0.5 also "  << endl;
				PDF_THREAD_UNLOCK
				returnValue = 0.5;
			}
			else {
				PDF_THREAD_LOCK
				cout << "Bs2JpsiPhi_Signal_v5::mistagBbar() : WARNING ******If you got here you dont know what you are doing  "  << endl;
				PDF_THREAD_UNLOCK
				exit(1);
			}
			return returnValue;
		}

		inline double D1() const {  return 1.0 - q()*(mistagB()-mistagBbar()); }
		inline double D2() const {  return q()*( 1.0 - mistagB() -mistagBbar() ); }
		//inline double D1() const {  return 1.0 ; }
		//inline double D2() const {  return  ( q()*( 1.0 - mistagB() -mistagBbar() ) ) / ( 1.0 - q()*(mistagB()-mistagBbar())  )  ; }


		//.....................
		// C, D, S
		inline double cosphis() const { return _DD ; } //  _cosphis ; }
		inline double sinphis() const { return _SS ; } //  _sinphis ; }
		inline double CC() const { return _CC ; } //  _sinphis ; }


		//......................................................
		// Time primitives

		inline double expL() const { return expL_stored; }
		inline double intExpL( ) const { return intExpL_stored; }

		inline double expH() const { return expH_stored; }
		inline double intExpH( ) const { return intExpH_stored; }

		inline double expSin() const  { return expSin_stored; }
		inline double intExpSin( ) const { return intExpSin_stored;  }

		inline double expCos() const { return expCos_stored; }
		inline double intExpCos( ) const { return intExpCos_stored; }



		//---------------------------------------------------------
		//............. Differential cross sections and normalisations
		double diffXsec();
		double diffXsecTimeOnly();
		double diffXsecNorm1();
		double diffXsecCompositeNorm1( int resolutionIndex );

		bool normalisationCacheValid;
		double normalisationCacheValue[3];
		//double normalisationCacheValueRes2[3] ;

		void DebugPrint( string , double ) const;
		void DebugPrintXsec( string , double ) const;
		void DebugPrintNorm( string , double ) const;


		//------------------------------------------------------------------------------
		// These are the time factors and their analytic integrals for the one angle PDF

		//..................................
		inline double timeFactorEven()  const
		{
			//if( t < 0.0 ) return 0.0 ;
			const double result =
				D1() * (
						( 1.0 + cosphis() ) * expL( )
						+ ( 1.0 - cosphis() ) * expH( )
				       ) +
				D2() * (
						( 2.0 * sinphis() ) * expSin( )
						+ ( 2.0 * CC()      ) * expCos( )
				       );

			//DEBUG
			if( DebugFlag_v5 )
			{
				if( result < 0 )
				{
					PDF_THREAD_LOCK
					cout << " Bs2JpsiPhi_Signal_v5::timeFactorEven() : result < 0 " << endl ;
					cout << " ->term1 " << ( 1.0 + cosphis() ) * expL( ) << endl ;
					cout << " ->term2 " << ( 1.0 - cosphis() ) * expH( ) << endl ;
					cout << " ->term3 " << q() * ( 2.0 * sinphis()   ) * expSin( ) * (1.0 - 2.0*mistag()) << endl ;
					cout << "   -->sin(phis) "  << sinphis() << endl ;
					cout << "   -->expSin    "  << expSin() << endl ;
					cout << "   -->tagFrac   "  << mistag() << endl ;
					cout << "   -->delta_ms  "  << delta_ms << endl ;
					PDF_THREAD_UNLOCK
				}
			}
			return result;
		}

		inline double timeFactorEvenInt()  const
		{
			return
				D1() * (
						( 1.0 + cosphis() )  * intExpL()
						+ ( 1.0 - cosphis() )  * intExpH()
				       ) +
				D2() * (
						( 2.0 * sinphis() ) * intExpSin( )
						+  ( 2.0 * CC()      ) * intExpCos( )
				       ) ;
		}


		//..................................
		inline double timeFactorOdd(  )   const
		{
			//if( t < 0.0 ) return 0.0 ;
			return
				D1() * (
						( 1.0 - cosphis() ) * expL( )
						+ ( 1.0 + cosphis() ) * expH( )
				       ) +
				D2() * (
						-  ( 2.0 * sinphis() ) * expSin( )
						+  ( 2.0 * CC()      ) * expCos( )
				       ) ;
		}

		inline double timeFactorOddInt(  )  const
		{
			return
				D1() * (
						( 1.0 - cosphis() ) * intExpL()
						+ ( 1.0 + cosphis() ) * intExpH()
				       ) +
				D2() * (
						-  ( 2.0 * sinphis() ) * intExpSin( )
						+  ( 2.0 * CC()      ) * intExpCos( )
				       ) ;
		}


		//----------------------------------------------------------
		// These are the time factors and their analytic integrals for the three angle PDF

		//...........................
		inline double timeFactorA0A0( )    const { return timeFactorEven( ) ; }
		inline double timeFactorA0A0Int( ) const { return timeFactorEvenInt( ) ; }

		//...........................
		inline double timeFactorAPAP( )    const { return timeFactorEven( ) ; }
		inline double timeFactorAPAPInt( ) const { return timeFactorEvenInt( ) ; }

		//...........................
		inline double timeFactorATAT( )    const { return timeFactorOdd( ) ; }
		inline double timeFactorATATInt( ) const { return timeFactorOddInt( ) ; }

		//...........................
		inline double timeFactorImAPAT( ) const
		{
			return
				D1() * (
						( expL( ) - expH( ) ) * cos_delta1 * sinphis()
						+ ( expL( ) + expH( ) ) * sin_delta1 * CC()
				       ) +
				D2() * (
						2.0  * ( sin_delta1*expCos( ) - cos_delta1*cosphis()*expSin( ) )
				       ) ;
		}

		inline double timeFactorImAPATInt( ) const
		{
			return
				D1() * (
						( intExpL() - intExpH() ) * cos_delta1 * sinphis()
						+ ( intExpL() + intExpH() ) * sin_delta1 * CC()
				       ) +
				D2() * (
						2.0  * ( sin_delta1*intExpCos() - cos_delta1*cosphis()*intExpSin() )
				       ) ;
		}


		//...........................
		inline double timeFactorReA0AP( )  const
		{
			if( _useCosDpar ) return cosdpar * this->timeFactorEven(  ) ;//PELC-COSDPAR Special for fitting cosdpar separately
			else return cos_delta_2_1 * this->timeFactorEven(  ) ;
		}

		inline double timeFactorReA0APInt( ) const
		{
			if( _useCosDpar ) return cosdpar * this->timeFactorEvenInt( ) ;//PELC-COSDPAR Special for fitting cosdpar separately
			else return cos_delta_2_1 * this->timeFactorEvenInt( ) ;
		}


		//...........................
		inline double timeFactorImA0AT(  ) const
		{
			return
				D1() * (
						( expL( ) - expH( ) ) * cos_delta2 * sinphis()
						+ ( expL( ) + expH( ) ) * sin_delta2 * CC()
				       ) +
				D2() * (
						2.0  * ( sin_delta2*expCos( ) - cos_delta2*cosphis()*expSin( ) )
				       ) ;
		}

		inline double timeFactorImA0ATInt( ) const
		{

			return
				D1() * (
						( intExpL() - intExpH()  ) * cos_delta2 * sinphis()
						+ ( intExpL() + intExpH()  ) * sin_delta2 * CC()
				       ) +
				D2() * (
						2.0  * ( sin_delta2*intExpCos() - cos_delta2*cosphis()*intExpSin()  )
				       ) ;
		}

		//.... S wave additions.......

		//...........................
		inline double timeFactorASAS( )    const { return timeFactorOdd( ) ; }
		inline double timeFactorASASInt( ) const { return timeFactorOddInt( ) ; }


		//...........................
		inline double timeFactorReASAP( ) const
		{
			return
				D1() * (
						( expL( ) - expH( ) ) * sin_delta_para_s * sinphis()
						+ ( expL( ) + expH( ) ) * cos_delta_para_s * CC()
				       ) +
				D2() * (
						2.0  * ( cos_delta_para_s*expCos( ) - sin_delta_para_s*cosphis()*expSin( ) )
				       ) ;
		}

		inline double timeFactorReASAPInt( ) const
		{
			return
				D1() * (
						( intExpL() - intExpH() ) * sin_delta_para_s * sinphis()
						+ ( intExpL() + intExpH() ) * cos_delta_para_s * CC()
				       ) +
				D2() * (
						2.0  * ( cos_delta_para_s*intExpCos() - sin_delta_para_s*cosphis()*intExpSin() )
				       ) ;
		}


		//...........................
		inline double timeFactorImASAT( )  const
		{
			return sin_delta_perp_s * this->timeFactorOdd(  ) ;
		}

		inline double timeFactorImASATInt( ) const
		{
			return sin_delta_perp_s * this->timeFactorOddInt( ) ;
		}

		//...........................
		inline double timeFactorReASA0( ) const
		{
			return
				D1() * (
						( expL( ) - expH( ) ) * sin_delta_zero_s * sinphis()
						+ ( expL( ) + expH( ) ) * cos_delta_zero_s * CC()
				       ) +
				D2() * (
						2.0  * ( cos_delta_zero_s*expCos( ) - sin_delta_zero_s*cosphis()*expSin( ) )
				       ) ;
		}

		inline double timeFactorReASA0Int( ) const
		{
			return
				D1() * (
						( intExpL() - intExpH() ) * sin_delta_zero_s * sinphis()
						+ ( intExpL() + intExpH() ) * cos_delta_zero_s * CC()
				       ) +
				D2() * (
						2.0  * ( cos_delta_zero_s*intExpCos() - sin_delta_zero_s*cosphis()*intExpSin() )
				       ) ;
		}

};

#endif

