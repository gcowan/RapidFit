// $Id: Bs2Jpsifzero_Signal_v6.h,v 1.1 2013/07/17  Dianne Ferguson Exp $
/** @class Bs2Jpsifzero_Signal_v6
 *
 *  Bs2Jpsifzero 
 *
 *  @author Dianne Ferguson dferguso@cern.ch
 *  @date 2013-07-17
 */

#ifndef Bs2Jpsifzero_Signal_v6_H
#define Bs2Jpsifzero_Signal_v6_H

#include <exception>

#include "BasePDF.h"
#include "PDFConfigurator.h"
#include "SlicedAcceptance.h"
#include "ResolutionModel.h"
#include "Mathematics.h"
#include <iostream>
#include <cstdlib>
#include <float.h>
#include <vector>

#include "TH1.h"
#include "TCanvas.h"

#define DOUBLE_TOLERANCE 1E-6
#define DEBUGFLAG true

using namespace::std;

//============================================================================================
class Bs2Jpsifzero_Signal_v6 : public BasePDF
{
	public:
		Bs2Jpsifzero_Signal_v6( PDFConfigurator* );
		Bs2Jpsifzero_Signal_v6( const Bs2Jpsifzero_Signal_v6& );
		~Bs2Jpsifzero_Signal_v6();

		// Mandatory RapidFit Methods
		virtual double EvaluateForNumericIntegral(DataPoint*) ;
		virtual double Evaluate(DataPoint*);
		virtual double EvaluateTimeOnly(DataPoint*) ;
		virtual bool SetPhysicsParameters(ParameterSet*);
		virtual vector<string> GetDoNotIntegrateList();

	private:
		void MakePrototypes();
		double normalisationCacheUntagged ;
		void prepareCDS();

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

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);


		//PELC For debugging purposes
		//TH1D * histOfPdfValues ;
		//TCanvas * c0 ;
		//mutable int histCounter ;
		//~PELC

		// Parameters
		ObservableRef gammaName;		// gamma
		ObservableRef deltaGammaName;	// delta gamma
		ObservableRef deltaMName;		// delta mass
		ObservableRef Aperp_sqName;	// amplitude
		ObservableRef CspName;		// amplitude

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

		// Observables
		ObservableRef timeName;		// proper time
		ObservableRef tagName;			// B tag

		PseudoObservable ATAT_Obs;

		double ATAT_value;

		// Measured Event Observables
		double t ;
		int tag ;

		// Physics Fit Parameters
		double _gamma ;
		double dgam ;

		double Aperp_sq ;
		double Csp ;

		double delta_ms ;
		double phi_s ;
		double _cosphis ;
		double _sinphis ;
		double lambda ;
		double _CC, _DD, _SS;

		double _mistag ;
		double _mistagP1 ;
		double _mistagP0 ;
		double _mistagSetPoint ;
		double _mistagDeltaP1 ;
		double _mistagDeltaP0 ;
		double _mistagDeltaSetPoint ;

		double timeOffset ;
		bool _useEventResolution ;
		inline bool useEventResolution() const {return _useEventResolution ; }
		inline bool useTimeAcceptance() const { return _useTimeAcceptance ; }

	        ResolutionModel * resolutionModel ;

		// Othere things calculated later on the fly
		double tlo, thi ;

		// stored time primitives
		mutable double expL_stored ;
		mutable double expH_stored ;
		mutable double expSin_stored ;
		mutable double expCos_stored ;
		mutable double intExpL_stored ;
		mutable double intExpH_stored ;
		mutable double intExpSin_stored ;
		mutable double intExpCos_stored ;
		void preCalculateTimeFactors();
		void preCalculateTimeIntegrals();

		bool timeIntegralCacheValid ;

		//Time acceptance
		SlicedAcceptance * timeAcc ;

		//Configurationparameters
		bool _useTimeAcceptance ;
		bool _numericIntegralForce ;
		bool _numericIntegralTimeOnly ;
		bool _useCosAndSin ;
		bool _useCosDpar ;
		bool _usePunziSigmat ;
		bool _usePunziMistag ;


		//....................................
		//Internal helper functions

		double stored_AT;
		inline double AT() const { return stored_AT; }

		//....................
		// Safety for gammas
		double stored_gammal;
		inline double gamma_l() const { return stored_gammal; }

		double stored_gammah;
		inline double gamma_h() const { return stored_gammah; }

		inline double gamma() const { return _gamma ; }

		inline double q() const { return tag ;}

		inline double mistag() const {
			double returnValue = -1000.;

			if( fabs((q()-0.0)) < DOUBLE_TOLERANCE ) {
				returnValue = 0.5 ;
			}
			else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
				//Normal case
				returnValue =  _mistagP0 + _mistagP1*(_mistag - _mistagSetPoint ) ;
				if( returnValue < 0 )  returnValue = 0 ;
				if( returnValue > 0.5) returnValue = 0.5 ;
			}
			else if( _mistag < 0.0 ) {
				PDF_THREAD_LOCK
				cout << "Bs2Jpsifzero_Signal_v5::mistag() : _mistag < 0 so set to 0 " << endl ;
				PDF_THREAD_UNLOCK
				returnValue = 0 ;
			}
			else if( _mistag > 0.5 ) {
				PDF_THREAD_LOCK
				cout << "Bs2Jpsifzero_Signal_v5::mistag() : _mistag > 0.5 so set to 0.5 "  << endl ;
				PDF_THREAD_UNLOCK
				returnValue = 0.5 ;
			}
			else {
				PDF_THREAD_LOCK
				cout << "Bs2Jpsifzero_Signal_v5::mistag() : WARNING ******If you got here you dont know what you are doing  "  << endl ;
				PDF_THREAD_UNLOCK
				exit(1);
			}
			return returnValue ;
		}

		inline double mistagB() const {
			double returnValue = -1000.;

			if( fabs(q()) < 0.5 ) {
				returnValue = 0.5 ;
			}
			else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
				//Normal case
				returnValue =  _mistagP0+(_mistagDeltaP0/2.0) + (_mistagP1+(_mistagDeltaP1/2.0))*(_mistag - (_mistagSetPoint+(_mistagDeltaSetPoint/2.0)) ) ;
				if( returnValue < 0 )  returnValue = 0 ;
				if( returnValue > 0.5) returnValue = 0.5 ;
			}
			else if( _mistag < 0.0 ) {
				PDF_THREAD_LOCK
				cout << "Bs2Jpsifzero_Signal_v6::mistagB() : _mistag < 0 so deltaMistag set to 0 also " << endl ;
				PDF_THREAD_UNLOCK
				returnValue = 0 ;
			}
			else if( _mistag > 0.5 ) {
				PDF_THREAD_LOCK
				cout << "Bs2Jpsifzero_Signal_v6::mistagB() : _mistag > 0.5 so so deltaMistag set to 0.5 also "  << endl ;
				PDF_THREAD_UNLOCK
				returnValue = 0.5 ;
			}
			else {
				PDF_THREAD_LOCK
				cout << "Bs2Jpsifzero_Signal_v6::mistagB() : WARNING ******If you got here you dont know what you are doing  "  << endl ;
				PDF_THREAD_UNLOCK
				exit(1);
			}
			return returnValue ;
		}

		inline double mistagBbar() const {
			double returnValue = -1000.;

			if( fabs(q()) < 0.5 ) {
				returnValue = 0.5 ;
			}
			else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
				//Normal case
				returnValue =  _mistagP0-(_mistagDeltaP0/2.0) + (_mistagP1-(_mistagDeltaP1/2.0))*(_mistag - (_mistagSetPoint-(_mistagDeltaSetPoint/2.0)) ) ;
				if( returnValue < 0 )  returnValue = 0 ;
				if( returnValue > 0.5) returnValue = 0.5 ;
			}
			else if( _mistag < 0.0 ) {
				PDF_THREAD_LOCK
				cout << "Bs2Jpsifzero_Signal_v6::mistagBbar() : _mistag < 0 so deltaMistag set to 0 also " << endl ;
				PDF_THREAD_UNLOCK
				returnValue = 0 ;
			}
			else if( _mistag > 0.5 ) {
				PDF_THREAD_LOCK
				cout << "Bs2Jpsifzero_Signal_v6::mistagBbar() : _mistag > 0.5 so so deltaMistag set to 0.5 also "  << endl ;
				PDF_THREAD_UNLOCK
				returnValue = 0.5 ;
			}
			else {
				PDF_THREAD_LOCK
				cout << "Bs2Jpsifzero_Signal_v6::mistagBbar() : WARNING ******If you got here you dont know what you are doing  "  << endl ;
				PDF_THREAD_UNLOCK
				exit(1);
			}
			return returnValue ;
		}

		inline double D1() const {  return 1.0 - q()*(mistagB()-mistagBbar()) ; }
		inline double D2() const {  return q()*( 1.0 - mistagB() -mistagBbar() ) ; }

		//.....................
		// C, D, S
		inline double cosphis() const { return _DD ; } //  _cosphis ; }
		inline double sinphis() const { return _SS ; } //  _sinphis ; }
		inline double CC() const { return _CC ; } //  _sinphis ; }


		//......................................................
		// Time primitives

		inline double expL() const { return expL_stored ; }
		inline double intExpL( ) const { return intExpL_stored ; }

		inline double expH() const { return expH_stored ; }
		inline double intExpH( ) const { return intExpH_stored ; }

		inline double expSin() const  { return expSin_stored ; }
		inline double intExpSin( ) const { return intExpSin_stored ;  }

		inline double expCos() const { return expCos_stored ; }
		inline double intExpCos( ) const { return intExpCos_stored ; }



		//---------------------------------------------------------
		//............. Differential cross sections and normalisations
		double diffXsec();
		double diffXsecTimeOnly();
		double diffXsecNorm1();
		double diffXsecCompositeNorm1( )  ;

		bool normalisationCacheValid ;
		double normalisationCacheValue[3] ;

		void DebugPrint( string , double ) const ;
		void DebugPrintXsec( string , double ) const ;
		void DebugPrintNorm( string , double ) const ;


		//------------------------------------------------------------------------------
		// These are the time factors and their analytic integrals for the one angle PDF

		//..................................
		inline double timeFactorEven()  const
		{
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
			if( DEBUGFLAG && (result < 0) ) {
				PDF_THREAD_LOCK
				cout << " Bs2Jpsifzero_Signal_v5::timeFactorEven() : result < 0 " << endl ;
				cout << " ->term1 " << ( 1.0 + cosphis() ) * expL( ) << endl ;
				cout << " ->term2 " << ( 1.0 - cosphis() ) * expH( ) << endl ;
				cout << " ->term3 " << q() * ( 2.0 * sinphis()   ) * expSin( ) * (1.0 - 2.0*mistag()) << endl ;
				cout << "   -->sin(phis) "  << sinphis() << endl ;
				cout << "   -->expSin    "  << expSin() << endl ;
				cout << "   -->tagFrac   "  << mistag() << endl ;
				cout << "   -->delta_ms  "  << delta_ms << endl ;
				PDF_THREAD_UNLOCK
			}
			return result ;
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
		inline double timeFactorATAT( )    const { return timeFactorOdd( ) ; }
		inline double timeFactorATATInt( ) const { return timeFactorOddInt( ) ; }
};

#endif

