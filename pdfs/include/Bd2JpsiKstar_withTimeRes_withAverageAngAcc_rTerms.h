// $Id: Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.h
 *
 *  RaKstarFlavourFit PDF for Bd2JpsiKstar
 *
 *  @author Ailsa Sparkes
 *  @date 2010-02-10
 */

#ifndef Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms_H
#define Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms_H

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif

class Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms : public BasePDF
{
	public:
		Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms();
		~Bd2JpsiKstar_withTimeRes_withAverageAngAcc_rTerms();

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
		double cachedAzeroAzeroIntB, cachedAparaAparaIntB, cachedAperpAperpIntB;
		double cachedAparaAperpIntB, cachedAzeroAparaIntB, cachedAzeroAperpIntB;

		double cachedSinDeltaPerpPara, cachedCosDeltaPara, cachedSinDeltaPerp;
		double cachedCosSqDeltaM;

		double cachedAzero, cachedApara, cachedAperp;
		bool normalisationCacheValid, evaluationCacheValid;
//		double q;
		// These contain the strings that correspond
		// to the physics parameter names that will be
		// used in the minimiser.
		string gammaName;		// gamma
		string deltaMName;		// delta mass
		//string Azero_sqName;	// amplitude
		//string Apara_sqName;	// amplitude
		//string Aperp_sqName;	// amplitude
		string R_alphaName;     // amplitude
                string R_betaName;      // amplitude
		string delta_zeroName;	// strong phase, set to 0
		string delta_paraName;	// strong phase
		string delta_perpName;	// strong phase
		string angAccI1Name;		// Pre-calculated angular integrals including acceptance
		string angAccI2Name;		//
		string angAccI3Name;		//
		string angAccI4Name;		//
		string angAccI5Name;		//
		string angAccI6Name;		//
		string timeRes1Name;
		string timeRes2Name;
		string timeRes1FractionName;

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF.
		string timeName;		// proper time
		string cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		string phiName;		// azimuthal angle of the mu+ in Jpsi frame
		string cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
		string KstarFlavourName;
		//kaonIDName;				// in phi rest frame


		// Member variables that will contain the parameter values
                double gamma;
                double deltaMs;
		 double R_alpha;
                double R_beta;
		double Azero_sq;
                double Apara_sq;
                double Aperp_sq;
                double AzeroApara;
                double AzeroAperp;
                double AparaAperp;
                double delta_zero;
                double delta_para;
                double delta_perp;
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


		// Member variables for the observables
		double time;
		double cosTheta;
		double phi;
		double cosPsi;
		double KstarFlavour;
		//int kaonID;
		double q() const ;

		double tlow, thigh; // Integration limits

		double buildPDFnumerator();
		double buildPDFdenominator();
		void getTimeDependentAmplitudes( double&, double&, double&, double&, double&, double&);
		void getTimeAmplitudeIntegrals(double&, double&, double&, double&, double&, double&);
		double AT() const ;
                double AP() const ;
                double A0() const ;
		double DT() const ;
                double DP() const ;
                double D0() const ;
		double Gamma() const ;

};


#endif
