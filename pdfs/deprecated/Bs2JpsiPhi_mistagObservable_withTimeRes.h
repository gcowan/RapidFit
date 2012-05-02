// $Id: Bs2JpsiPhi_mistagObservable_withTimeRes.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhi_mistagObservable_withTimeRes Bs2JpsiPhi_mistagObservable_withTimeRes.h
 *
 *  RapidFit PDF for Bs2JpsiPhi with time resolution
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-07-30
 */

#ifndef Bs2JpsiPhi_mistagObservable_withTimeRes_H
#define Bs2JpsiPhi_mistagObservable_withTimeRes_H

#include "BasePDF.h"
#include "RooComplex.h"

class Bs2JpsiPhi_mistagObservable_withTimeRes : public BasePDF
{
	public:
		Bs2JpsiPhi_mistagObservable_withTimeRes();
		~Bs2JpsiPhi_mistagObservable_withTimeRes();

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
		double cachedv1, cachedv2;
		double cachedAzero, cachedApara, cachedAperp, cachedsinDeltaPerpPara, cachedcosDeltaPerpPara, cachedsinDeltaPerp;
		double cachedcosDeltaPerp, cachedcosDeltaPara, cachedsinPhis, cachedcosPhis;
		double cachedExpCosh, cachedExpSinh, cachedExpCos, cachedExpSin;	
		bool normalisationCacheValid, evaluationCacheValid;

		// These contain the strings that correspond
		// to the physics parameter names that will be used in the minimiser.
		string gammaName;	// gamma
		string deltaGammaName;	// delta gamma
		string deltaMName;	// delta mass
		string Phi_sName;	// what we want to measure!
		string Azero_sqName;	// amplitude
		string Apara_sqName;	// amplitude
		string Aperp_sqName;	// amplitude
		string delta_zeroName;	// strong phase, set to 0
		string delta_paraName;	// strong phase
		string delta_perpName;	// strong phase
		string mean_time_res1Name; // propertime resolution parameters. Double Gaussian used.
        	string mean_time_res2Name;
        	string sigma_time_res1Name;
        	string sigma_time_res2Name;
        	string frac_time_res1Name;

		// These contain the strings that correspond
		// to the observable names that are used in the PDF. 
		string timeName;	// proper time
		string timeErrName;	// proper time
		string cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		string phiName;		// azimuthal angle of the mu+ in Jpsi frame
		string cosPsiName;	// helicity angle between K+ and -ve Jpsi direction in phi rest frame
		string tagName;		// B tag
		string mistagName;	// B mistag
	
		void getPhysicsParameters( double&, double&, double&, double&, double&
					 , double&, double&, double&, double&, double&
					 , double&, double&, double&, double&, double&);

		// Function to get the 6 angular and 6 time dependent terms
		void getAngularFunctions( double&, double&, double&, double&, double&, double&, DataPoint*);
		void getTimeDependentAmplitudes( double&, double&, double&, double&, double&, double&, DataPoint*, int);

		// The non-zero time integrals components
		void getTimeAmplitudeIntegrals(double&, double&, double&, double&, double&, double&, PhaseSpaceBoundary*, int, double);
		inline double getAzeroAzeroInt(double, double, double, double, double, int);
		inline double getAparaAparaInt(double, double, double, double, double, int);
		inline double getAperpAperpInt(double, double, double, double, double, int);
		inline double getAparaAperpInt(double, double, double, double, double, int, double, double);
		inline double getAzeroAparaInt(double, double, double, double, double, int, double);
		inline double getAzeroAperpInt(double, double, double, double, double, int, double);

		inline void getExpTrigIntegrals(double&, double&, double&, double&, double, 
						double, double, double, double, double, 
						double, double, double, double);
		inline void getErfPart(double, double, double, double, double, double, double, double&, double&, double&, double&);

		// Need these for the analytic convolution of the propertime resolution	
		inline double ExpCosh(double, double, double, double, double, double, double, double);
		inline double ExpSinh(double, double, double, double, double, double, double, double);
		inline double ExpCos(double, double, double, double, double, double, double, double);
		inline double ExpSin(double, double, double, double, double, double, double, double);
		
		// Helper functions for including the propertime resolution
		inline double erfc(double, double, double, double, double, int);
		inline double cerfRe(double, double, double, double, double);
		inline double cerfIm(double, double, double, double, double);

		// Some more helpers
		inline double evalCerfRe(double, double, double);
		inline double evalCerfIm(double, double, double);
		inline RooComplex evalCerf(double, double, double);
		inline RooComplex evalCerfApprox(double, double, double);
};

#endif
