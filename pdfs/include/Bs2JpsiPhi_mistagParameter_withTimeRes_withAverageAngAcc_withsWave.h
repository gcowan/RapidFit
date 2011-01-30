// $Id: Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave.h,v 1.0 2011/01/28 10:35:49 cofitzpa Exp $
/** @class Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave.h
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Conor Fitzpatrick conor.fitzpatrick@cern.ch
 *  @date 2011-01-28
 */

#ifndef Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_H
#define Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave_H

#include "BasePDF.h"

class Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave : public BasePDF
{
	public:
		Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave();
		~Bs2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc_withsWave();

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
		bool normalisationCacheValid, evaluationCacheValid;

		double cachedExpCosh, cachedExpSinh, cachedExpCos, cachedExpSin;



		// These contain the strings that correspond
		// to the physics parameter names that will be
		// used in the minimiser.
		string gammaName;		// gamma
		string deltaGammaName;	// delta gamma
		string deltaMName;		// delta mass
		string Phi_sName;		// what we want to measure!
		string Azero_sqName;	// amplitude
		string Aperp_sqName;	// amplitude
		string As_sqName;	// amplitude
		string delta_zeroName;	// strong phase, set to 0
		string delta_paraName;	// strong phase
		string delta_perpName;	// strong phase
		string delta_sName;	// strong phase
		string angAccI1Name;		// Pre-calculated angular integrals including acceptance
		string angAccI2Name;		// 
		string angAccI3Name;		// 
		string angAccI4Name;		// 
		string angAccI5Name;		// 
		string angAccI6Name;		// 
		string angAccI7Name;		// 
		string angAccI8Name;		// 
		string angAccI9Name;		// 
		string angAccI10Name;		// 
		string mistagName;
		string timeRes1Name;
		string timeRes2Name;
		string timeRes1FractionName;

		// Member variables that will contain the parameter values
                double gamma;
                double deltaGamma;
                double deltaMs;
                double Phi_s;
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


		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		string timeName;		// proper time
		string cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		string phiName;		// azimuthal angle of the mu+ in Jpsi frame
		string cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
					// in phi rest frame
		string tagName;		// B tag

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
