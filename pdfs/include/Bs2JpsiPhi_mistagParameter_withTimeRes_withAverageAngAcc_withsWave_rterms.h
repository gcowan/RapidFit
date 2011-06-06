// $Id: Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms.h,v 1.0 2011/01/28 10:35:49 cofitzpa Exp $
/** @class Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave.h
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Conor Fitzpatrick conor.fitzpatrick@cern.ch
 *  @date 2011-01-28
 */

#ifndef Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms_H
#define Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif

class Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms : public BasePDF
{
	public:
		Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms();
		~Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_rterms();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);
		virtual bool SetPhysicsParameters(ParameterSet*);
		//Return a list of parameters not to be integrated
                virtual vector<string> GetDoNotIntegrateList();

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();

		//Cached values
		double cachedAzeroAzeroIntB, cachedAparaAparaIntB, cachedAperpAperpIntB, cachedAparaAperpIntB, cachedAzeroAparaIntB , cachedAzeroAperpIntB , cachedAsAsIntB, cachedAsAparaIntB, cachedAsAperpIntB, cachedAsAzeroIntB;
		double cachedAzeroAzeroIntBbar, cachedAparaAparaIntBbar, cachedAperpAperpIntBbar, cachedAparaAperpIntBbar, cachedAzeroAparaIntBbar, cachedAzeroAperpIntBbar, cachedAsAsIntBbar, cachedAsAparaIntBbar, cachedAsAperpIntBbar, cachedAsAzeroIntBbar;

		double cachedSinDeltaPerpPara, cachedCosDeltaPerpPara, cachedSinDeltaPerp, cachedCosDeltaPerp, cachedCosDeltaPara, cachedSinDeltaPerpS, cachedSinDeltaParaS, cachedCosDeltaParaS, cachedSinDeltaZeroS, cachedCosDeltaZeroS, cachedSinPhis, cachedCosPhis;

		double cachedExpCosh, cachedExpSinh, cachedExpCos, cachedExpSin;


		bool normalisationCacheValid, evaluationCacheValid;

		// These contain the ObservableRefs that correspond
		// to the physics parameter names that will be
		// used in the minimiser.
		ObservableRef gammaName;		// gamma
		ObservableRef deltaGammaName;	// delta gamma
		ObservableRef deltaMName;		// delta mass
		ObservableRef Phi_sName;		// what we want to measure!
		ObservableRef R_alphaName;	// amplitude
		ObservableRef R_betaName;	// amplitude
		ObservableRef R_gammaName;	// amplitude
		ObservableRef delta_zeroName;	// strong phase, set to 0
		ObservableRef delta_paraName;	// strong phase
		ObservableRef delta_perpName;	// strong phase
		ObservableRef delta_sName;	// strong phase
		ObservableRef angAccI1Name;		// Pre-calculated angular integrals including acceptance
		ObservableRef angAccI2Name;		// 
		ObservableRef angAccI3Name;		// 
		ObservableRef angAccI4Name;		// 
		ObservableRef angAccI5Name;		// 
		ObservableRef angAccI6Name;		// 
		ObservableRef angAccI7Name;		// 
		ObservableRef angAccI8Name;		// 
		ObservableRef angAccI9Name;		// 
		ObservableRef angAccI10Name;		// 
		ObservableRef mistagName;
		ObservableRef timeRes1Name;
		ObservableRef timeRes2Name;
		ObservableRef timeRes1FractionName;

		// These contain the ObservableRefs that correspond
		// to the observable names that are used in the
		// PDF.
		ObservableRef timeName;		// proper time
		ObservableRef cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		ObservableRef phiName;		// azimuthal angle of the mu+ in Jpsi frame
		ObservableRef cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
		// in phi rest frame
		ObservableRef tagName;		// B tag

		ObservableRef timeconstraintName;

		// Member variables that will contain the parameter values
                double gamma;
                double deltaGamma;
                double deltaMs;
                double Phi_s;
		double R_alpha;
		double R_beta;
		double R_gamma;
                double Azero_sq;
                double Apara_sq;
                double Aperp_sq;
		double As_sq;
                double AzeroApara;
                double AzeroAperp;
                double AparaAperp;
		double AsApara;
		double AsAperp;
		double AsAzero;
                double delta_zero;
                double delta_para;
                double delta_perp;
		double delta_s;
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

		// Member variables for the observables
		double time;
		double cosTheta;
		double phi;
		double cosPsi;
		int q; // flavour tag

		double tlow, thigh; // Integration limits
	
		double buildPDFnumerator();
		double buildPDFdenominator();
		void getTimeDependentAmplitudes( double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, int);
		void getTimeAmplitudeIntegrals(double&, double&, double&, double&, double&, double&, double&, double&, double&, double&, int);
		void evaluateCache();
};

#endif
