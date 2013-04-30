/** @class Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.h
 *
 *  RapidFit PDF for Bd2JpsiKstar
 *
 *  @author Ailsa Sparkes
 *  @date 2011-01-26
 */

#ifndef Bd2JpsiKstar_sWave_Fscopy_H
#define Bd2JpsiKstar_sWave_Fscopy_H

//#ifndef __CINT__
#include "BasePDF.h"
#include "ComponentRef.h"
#include "SlicedAcceptance.h"
#include "Bs2JpsiPhi_Angluar_Terms.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TAxis.h"
#include "TH3D.h"
//#endif
//#ifdef __CINT__
//#endif

class Bd2JpsiKstar_sWave_Fscopy : public BasePDF
{
	public:
		Bd2JpsiKstar_sWave_Fscopy( PDFConfigurator* );
		Bd2JpsiKstar_sWave_Fscopy( const Bd2JpsiKstar_sWave_Fscopy& );
		~Bd2JpsiKstar_sWave_Fscopy();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);
		virtual double EvaluateComponent(DataPoint*, ComponentRef*);
		virtual bool SetPhysicsParameters(ParameterSet*);
		double NormAnglesOnlyForAcceptanceWeights(DataPoint*, PhaseSpaceBoundary*);
		//Return a list of parameters not to be integrated
        virtual vector<string> GetDoNotIntegrateList();
		double angularFactor();
		string GetComponentName( ComponentRef* Component );
		vector<string> PDFComponents();

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		//Cached values
		double cachedAzeroAzeroIntB, cachedAparaAparaIntB, cachedAperpAperpIntB;
		double cachedAparaAperpIntB, cachedAzeroAparaIntB, cachedAzeroAperpIntB;
		double cachedAsAsIntB, cachedAparaAsIntB, cachedAperpAsIntB, cachedAzeroAsIntB ;
		double AzeroAzeroB, AparaAparaB, AperpAperpB, AsAsB;
		double ImAparaAperpB, ReAzeroAparaB, ImAzeroAperpB;
        double ReAparaAsB, ImAperpAsB, ReAzeroAsB;

		double cachedSinDeltaPerpPara, cachedCosDeltaPara, cachedSinDeltaPerp, cachedCosDeltaParaS, cachedSinDeltaPerpS, cachedCosDeltaS;


		double cachedAzero, cachedApara, cachedAperp, cachedAs;
        int componentIndex;

		// These contain the strings that correspond
		// to the physics parameter names that will be
		// used in the minimiser.
		ObservableRef gammaName;		// gamma
		ObservableRef As_sqName;
		ObservableRef delta_zeroName;	// strong phase, set to 0
		ObservableRef Rzero_sqName;	// amplitude
		ObservableRef Rpara_sqName;	// amplitude
		ObservableRef Rperp_sqName;	// amplitude
		ObservableRef delta_paraName;	// strong phase
		ObservableRef delta_perpName;
		ObservableRef delta_sName;
		ObservableRef CspName;
		ObservableRef angAccI1Name;		// Pre-calculated angular integrals including acceptance
		ObservableRef angAccI2Name;		//
		ObservableRef angAccI3Name;		//
		ObservableRef angAccI4Name;		//
		ObservableRef angAccI5Name;		//
		ObservableRef angAccI6Name;
		ObservableRef angAccI7Name;
		ObservableRef angAccI8Name;
		ObservableRef angAccI9Name;
		ObservableRef angAccI10Name;
		ObservableRef timeRes1Name;
		ObservableRef timeRes2Name;
		ObservableRef timeRes1FractionName;

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF.

		bool normalisationCacheValid, evaluationCacheValid;

		ObservableRef timeName;		// proper time
		ObservableRef cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		ObservableRef phiName;		// azimuthal angle of the mu+ in Jpsi frame
		ObservableRef cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
        ObservableRef helcosthetaKName;
        ObservableRef helcosthetaLName;
        ObservableRef helphiName;
		ObservableRef KstarFlavourName;

		ObservableRef timeconstraintName;
		// Member variables that will contain the parameter values
        double gamma;
        double Rzero_sq;
        double Rpara_sq;
        double Rperp_sq;
        double Azero_sq;
        double Apara_sq;
        double Aperp_sq;
		double As_sq;
        double AzeroApara;
        double AzeroAperp;
        double AparaAperp;
		double AparaAs;
		double AperpAs;
		double AzeroAs;
        double delta_zero;
        double delta_para;
        double delta_perp;
		double delta_s;
		double Csp;
		double CspAs;
        double omega;
        double timeRes;	 // This is the member variable used in the "builder" functions
        double timeRes1; // These are the physics parameters varied in the fit and passed from the XML
        double timeRes2;
        double timeRes1Frac;
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
		double Ap_sq, Ap;

        bool _useNumericalNormalisation ;
        bool _useHelicityBasis ;
        inline bool useHelicityBasis() const { return _useHelicityBasis ; }
		bool _useTimeAcceptance ;
		inline bool useTimeAcceptance() const { return _useTimeAcceptance ; }
		SlicedAcceptance * timeAcc ;
        bool _plotAllComponents ;

		 double tlo, thi ;
		// Member variables for the observables
		double time;
		double cosTheta;
		double phi;
		double cosPsi;
        double helcosthetaL ;
        double helcosthetaK ;
        double helphi ;
		double KstarFlavour;
		double q() const ;
        // Acceptance parameterisation
        static const int i_max = 4;
        static const int j_max = 4;
        static const int k_max = 4;
        double c[i_max+1][k_max+1][j_max+1];


    /* PELC Work in progress
        // PseudoObservables which allow for some calculations to be cached per event during Runtime
        vector<ObservableRef> angularTermDependencies;
        PseudoObservable f1_Obs;
        PseudoObservable f2_Obs;
        PseudoObservable f3_Obs;
        PseudoObservable f1_Obs;
        PseudoObservable f1_Obs;
        PseudoObservable f1_Obs;
        PseudoObservable f1_Obs;
        PseudoObservable f1_Obs;
        PseudoObservable f1_Obs;
        PseudoObservable f1_Obs;
    */
		//double tlow, thigh; // Integration limits

		TH3D *histo;
        TAxis *xaxis, *yaxis, *zaxis;
        int nxbins, nybins, nzbins;
        double xmin, xmax, ymin, ymax, zmin, zmax, deltax, deltay, deltaz;
        double total_num_entries;


		double buildPDFnumerator();
		double buildPDFdenominator();
		double buildCompositePDFdenominator();
		double buildPDFdenominatorAngles();
		void getTimeDependentAmplitudes( double&, double&, double&, double&, double&, double&, double&, double&, double&, double&);
		void getTimeAmplitudeIntegrals(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&);
		bool useFlatAngularDistribution;

        void getAngularFunctionsHelicity( double &,double &,double &,double &,double &,double &, double &, double &, double &, double &, double, double , double  );
        void getAngularFunctionsTransversity( double &,double &,double &,double &,double &,double &, double &, double &, double &, double &, double, double , double  );

		DataPoint* _datapoint;
};

#endif
