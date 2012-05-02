// $Id: Bs2JpsiPhi_Signal_v5_Robs.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2JpsiPhi_Signal_v5_Robs
 *
 *  Bs2JpsiPhi_SignalAlt series with mistag as observable
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2011-11-05
 */

#ifndef Bs2JpsiPhi_Signal_v5_Robs_H
#define Bs2JpsiPhi_Signal_v5_Robs_H

#include <exception>

#include "BasePDF.h"
#include "PDFConfigurator.h"
#include "SlicedAcceptance.h"
#include "AngularAcceptance.h"
#include "Mathematics.h"
#include "RooComplex.h"
#include <iostream>
#include <cstdlib>
#include <float.h>
#include <vector>

#include "TH1.h"
#include "TCanvas.h"

#define DOUBLE_TOLERANCE 1E-6
#define DEBUGFLAG true


//============================================================================================
class Bs2JpsiPhi_Signal_v5_Robs : public BasePDF
{
	public:
		Bs2JpsiPhi_Signal_v5_Robs( PDFConfigurator* );
		//Bs2JpsiPhi_Signal_v5_Robs( const Bs2JpsiPhi_Signal_v5_Robs& );
		~Bs2JpsiPhi_Signal_v5_Robs();

		//Mandatory RapidFit Methods
		virtual double EvaluateForNumericIntegral(DataPoint*);
		virtual double Evaluate(DataPoint*);
		virtual double EvaluateTimeOnly(DataPoint*);
		virtual bool SetPhysicsParameters(ParameterSet*);
		virtual vector<string> GetDoNotIntegrateList();

		vector<string> PDFComponents();

		double EvaluateComponent( DataPoint* input, ComponentRef* );

	private:
		Bs2JpsiPhi_Signal_v5_Robs& operator=( const Bs2JpsiPhi_Signal_v5_Robs& );
		void MakePrototypes();
		double normalisationCacheUntagged;

		void SetupAngularTerms();

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

		int componentIndex;

		//PELC For debugging purposes
		//TH1D * histOfPdfValues;
		//TCanvas * c0;
		//mutable int histCounter;
		//~PELC

		// Parameters
		ObservableRef gammaName;		// gamma
		ObservableRef deltaGammaName;	// delta gamma
		ObservableRef deltaMName;		// delta mass
		ObservableRef Azero_sqName;	// amplitude
		ObservableRef Apara_sqName;	// amplitude
		ObservableRef Aperp_sqName;	// amplitude
		ObservableRef As_sqName;		// amplitude
		ObservableRef delta_zeroName;	// strong phase, set to 0
		ObservableRef delta_paraName;	// strong phase
		ObservableRef delta_perpName;	// strong phase
		ObservableRef delta_sName;		// strong phase for S-wave
		ObservableRef cosdparName;		//PELC-COSDPAR Special for fitting cosdpar separately


		ObservableRef Phi_sName;		// what we want to measure!
		ObservableRef cosphisName;		// fitting cosphis and sinphis independently
		ObservableRef sinphisName;		// fitting cosphis and sinphis independently

		ObservableRef mistagName;		// mistag fraction  - may be used as observable also
		ObservableRef mistagP1Name;		// mistag calib
		ObservableRef mistagP0Name;		// mistag calib
		ObservableRef mistagSetPointName;// mistag calib

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

		double _mistag;
		double _mistagP1;
		double _mistagP0;
		double _mistagSetPoint;

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
		bool useEventResolution() const {return _useEventResolution; }

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
		void preCalculateTimeFactors() const;
		void preCalculateTimeIntegrals() const;

		bool timeIntegralCacheValid;
		vector< vector<double> > storeExpL;
		vector< vector<double> > storeExpH;
		vector< vector<double> > storeExpSin;
		vector< vector<double> > storeExpCos;
		void CacheTimeIntegrals();
		void deCacheTimeIntegrals( unsigned int ires, unsigned int islice );

		//Time acceptance
		SlicedAcceptance * timeAcc;

		//Configurationparameters
		bool _useTimeAcceptance;
		bool _numericIntegralForce;
		bool _numericIntegralTimeOnly;
		bool _useCosAndSin;
		bool _useCosDpar;
		bool allowNegativeAsSq;

		//....................................
		//Internal helper functions

		inline double AT() const {
			if( Aperp_sq <= 0. ) return 0.;
			else return sqrt(Aperp_sq);
		}
		inline double AP() const {
			if( Apara_sq <= 0. ) return 0.;
			else return sqrt(Apara_sq);
		}
		inline double A0() const {
			if( Azero_sq <= 0. ) return 0.;
			else return sqrt(Azero_sq);
		}
		inline double AS() const {
			if( As_sq <= 0. ) return 0.;
			else return sqrt(As_sq);
		}

		inline double gamma_l() const {
			const double gl = gamma() + ( dgam *0.5 );
			if( gl < 0. ) {
				cerr << " In Bs2JpsiPhi_SignalAlt_BaseClass_v4 : gamma_l() < 0 so setting it to 0.0000001 " << endl;
				return 0.0000001;
			}
			else
				return gl;
		}

		inline double gamma_h() const {
			const double gh = gamma() - ( dgam *0.5 );
			if( gh < 0. ) {
				cerr << " In Bs2JpsiPhi_SignalAlt_BaseClass_v4 : gamma_h() < 0 so setting it to 0.0000001 " << endl;
				return 0.0000001;
			}
			else
				return gh;
		}

		inline double gamma() const { return _gamma; }

		inline double q() const { return tag;}

		inline double mistag() const {
			double returnValue = -1000.;

			if( fabs((q()-0.0)) < DOUBLE_TOLERANCE ) {
				returnValue = 0.5;
			}
			else if( (_mistag>=0.0) && (_mistag <= 0.5) ) {
				//Normal case
				returnValue =  _mistagP0 + _mistagP1*(_mistag - _mistagSetPoint );
				if( returnValue < 0 )  returnValue = 0;
				if( returnValue > 0.5) returnValue = 0.5;
			}
			else if( _mistag < 0.0 ) {
				cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistag() : _mistag < 0 so set to 0 " << endl;
				returnValue = 0;
			}
			else if( _mistag > 0.5 ) {
				cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistag() : _mistag > 0.5 so set to 0.5 "  << endl;
				returnValue = 0.5;
			}
			else {
				cout << "Bs2JpsiPhi_SignalAlt_BaseClass_v4::mistag() : WARNING ******If you got here you dont know what you are doing  "  << endl;
				exit(1);
			}
			return returnValue;
		}

		inline double cosphis() const { return _cosphis; }
		inline double sinphis() const { return _sinphis; }

		inline bool useTimeAcceptance() const { return _useTimeAcceptance; }

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
		double diffXsec(  )  const;
		double diffXsecTimeOnly(  ) const;
		double diffXsecNorm1(  ) const;
		double diffXsecCompositeNorm1( int resolutionIndex ) ;

		bool normalisationCacheValid;
		double normalisationCacheValue[3];
		//double normalisationCacheValueRes2[3];

		void DebugPrint( string , double ) const;
		void DebugPrintXsec( string , double ) const;
		void DebugPrintNorm( string , double ) const;


		//------------------------------------------------------------------------------
		// These are the time factors and their analytic integrals for the one angle PDF

		//..................................
		inline double timeFactorEven(  )  const
		{
			//if( t < 0.0 ) return 0.0;
			const double result =
				( 1.0 + cosphis() ) * expL( )
				+ ( 1.0 - cosphis() ) * expH( )
				+ q() * ( 2.0 * sinphis()   ) * expSin( ) * (1.0 - 2.0*mistag());

			//DEBUG
			if( DEBUGFLAG && (result < 0) ) {
				cout << " Bs2JpsiPhi_SignalAlt_BaseClass_v4::timeFactorEven() : result < 0 " << endl;
				cout << " ->term1 " << ( 1.0 + cosphis() ) * expL( ) << endl;
				cout << " ->term2 " << ( 1.0 - cosphis() ) * expH( ) << endl;
				cout << " ->term3 " << q() * ( 2.0 * sinphis()   ) * expSin( ) * (1.0 - 2.0*mistag()) << endl;
				cout << "   -->sin(phis) "  << sinphis() << endl;
				cout << "   -->expSin    "  << expSin() << endl;
				cout << "   -->tagFrac   "  << mistag() << endl;
				cout << "   -->delta_ms  "  << delta_ms << endl;
			}
			return result;
		}

		inline double timeFactorEvenInt(  )  const
		{
			return
				( 1.0 + cosphis() )  * intExpL()
				+ ( 1.0 - cosphis() )  * intExpH()
				+ q() * ( 2.0 * sinphis()   ) * intExpSin( ) * (1.0 - 2.0*mistag());
		}


		//..................................
		inline double timeFactorOdd(  )   const
		{
			//if( t < 0.0 ) return 0.0;
			return
				( 1.0 - cosphis() ) * expL( )
				+ ( 1.0 + cosphis() ) * expH( )
				- q() * ( 2.0 * sinphis()   ) * expSin( ) * (1.0 - 2.0*mistag());
		}

		inline double timeFactorOddInt(  )  const
		{
			return
				( 1.0 - cosphis() ) * intExpL()
				+ ( 1.0 + cosphis() ) * intExpH()
				- q() * ( 2.0 * sinphis()   ) * intExpSin( ) * (1.0 - 2.0*mistag());
		}


		//----------------------------------------------------------
		// These are the time factors and their analytic integrals for the three angle PDF

		//...........................
		inline double timeFactorA0A0( )    const { return timeFactorEven( ); }
		inline double timeFactorA0A0Int( ) const { return timeFactorEvenInt( ); }

		//...........................
		inline double timeFactorAPAP( )    const { return timeFactorEven( ); }
		inline double timeFactorAPAPInt( ) const { return timeFactorEvenInt( ); }

		//...........................
		inline double timeFactorATAT( )    const { return timeFactorOdd( ); }
		inline double timeFactorATATInt( ) const { return timeFactorOddInt( ); }

		//...........................
		inline double timeFactorImAPAT( ) const
		{
			return
				q() * 2.0  * ( sin(delta1)*expCos( ) - cos(delta1)*cosphis()*expSin( ) ) * (1.0 - 2.0*mistag())
				- ( expH( ) - expL( ) ) * cos(delta1) * sinphis() ;
		}

		inline double timeFactorImAPATInt( ) const
		{
			//double _tlo = tlo;
			//if(_tlo < 0.) _tlo = 0.;

			return
				q() * 2.0  * ( sin(delta1)*intExpCos() - cos(delta1)*cosphis()*intExpSin() ) * (1.0 - 2.0*mistag())
				- ( intExpH() - intExpL() ) * cos(delta1) * sinphis();
		}


		//...........................
		inline double timeFactorReA0AP( )  const
		{
			if( _useCosDpar ) return cosdpar * this->timeFactorEven(  );//PELC-COSDPAR Special for fitting cosdpar separately
			else return cos(delta2-delta1) * this->timeFactorEven(  );
		}

		inline double timeFactorReA0APInt( ) const
		{
			if( _useCosDpar ) return cosdpar * this->timeFactorEvenInt( );//PELC-COSDPAR Special for fitting cosdpar separately
			else return cos(delta2-delta1) * this->timeFactorEvenInt( );
		}


		//...........................
		inline double timeFactorImA0AT(  ) const
		{
			return
				q() * 2.0  * ( sin(delta2)*expCos( ) - cos(delta2)*cosphis()*expSin( ) ) * (1.0 - 2.0*mistag())
				- ( expH( ) - expL( ) ) * cos(delta2) * sinphis();
		}

		inline double timeFactorImA0ATInt( ) const
		{
			//double _tlo = tlo;
			//if(_tlo < 0.) _tlo = 0.;

			return
				q() * 2.0  * ( sin(delta2)*intExpCos() - cos(delta2)*cosphis()*intExpSin()  ) * (1.0 - 2.0*mistag())
				- ( intExpH() - intExpL()  ) * cos(delta2) * sinphis();
		}

		//.... S wave additions.......

		//...........................
		inline double timeFactorASAS( )    const { return timeFactorOdd( ); }
		inline double timeFactorASASInt( ) const { return timeFactorOddInt( ); }


		//...........................
		inline double timeFactorReASAP( ) const
		{
			double delta = delta_para - delta_s;
			return
				q() * 2.0  * ( cos(delta)*expCos( ) - sin(delta)*cosphis()*expSin( ) ) * (1.0 - 2.0*mistag())
				- ( expH( ) - expL( ) ) * sin(delta) * sinphis() ;
		}

		inline double timeFactorReASAPInt( ) const
		{
			//double _tlo = tlo;
			//if(_tlo < 0.) _tlo = 0.;

			double delta = delta_para - delta_s;

			return
				q() * 2.0  * ( cos(delta)*intExpCos() - sin(delta)*cosphis()*intExpSin() ) * (1.0 - 2.0*mistag())
				- ( intExpH() - intExpL() ) * sin(delta) * sinphis();
		}


		//...........................
		inline double timeFactorImASAT( )  const
		{
			return sin(delta_perp-delta_s) * this->timeFactorOdd(  );
		}

		inline double timeFactorImASATInt( ) const
		{
			return sin(delta_perp-delta_s) * this->timeFactorOddInt( );
		}


		//...........................
		inline double timeFactorReASA0( ) const
		{
			double delta = delta_zero - delta_s;
			return
				q() * 2.0  * ( cos(delta)*expCos( ) - sin(delta)*cosphis()*expSin( ) ) * (1.0 - 2.0*mistag())
				- ( expH( ) - expL( ) ) * sin(delta) * sinphis() ;
		}

		inline double timeFactorReASA0Int( ) const
		{
			//double _tlo = tlo;
			//if(_tlo < 0.) _tlo = 0.;

			double delta = delta_zero - delta_s;

			return
				q() * 2.0  * ( cos(delta)*intExpCos() - sin(delta)*cosphis()*intExpSin() ) * (1.0 - 2.0*mistag())
				- ( intExpH() - intExpL() ) * sin(delta) * sinphis();
		}


};

#endif

