// $Id: Bs2JpsiPhi_mistagObservable.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhi_mistagObservable Bs2JpsiPhi_mistagObservable.h
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-07-30
 */

#ifndef Bs2JpsiPhi_mistagObservable_H
#define Bs2JpsiPhi_mistagObservable_H

#include "BasePDF.h"

class Bs2JpsiPhi_mistagObservable : public BasePDF
{
	public:
		Bs2JpsiPhi_mistagObservable();
		~Bs2JpsiPhi_mistagObservable();

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
		bool normalisationCacheValid, evaluationCacheValid;

		// These contain the strings that correspond
		// to the physics parameter names that will be
		// used in the minimiser.
		string gammaName;		// gamma
		string deltaGammaName;	// delta gamma
		string deltaMName;		// delta mass
		string Phi_sName;		// what we want to measure!
		string Azero_sqName;	// amplitude
		string Apara_sqName;	// amplitude
		string Aperp_sqName;	// amplitude
		string delta_zeroName;	// strong phase, set to 0
		string delta_paraName;	// strong phase
		string delta_perpName;	// strong phase

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		string timeName;		// proper time
		string cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		string phiName;		// azimuthal angle of the mu+ in Jpsi frame
		string cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
					// in phi rest frame
		string tagName;		// B tag
		string mistagName;		// B mistag
	
		void getPhysicsParameters( double&, double&, double&, double&, double&, double&, double&, double&, double&, double&);
		void getAngularFunctions( double&, double&, double&, double&, double&, double&, DataPoint*);
		
		void getTimeDependentAmplitudes( double&, double&, double&, double&, double&, double&, DataPoint*, int);
		void getTimeAmplitudeIntegrals(double&, double&, double&, PhaseSpaceBoundary*, int);

		inline double getAzeroAzeroInt(double, double, double, double, double, double, double, double, int);
		inline double getAparaAparaInt(double, double, double, double, double, double, double, double, int);
		inline double getAperpAperpInt(double, double, double, double, double, double, double, double, int);
	
        	// Work out what interference terms these correspond to.        
		inline double A4def(const double, const double, double, double, double, double, double, double, double, double, int);
		inline double A5def(const double, const double, double, double, double, double, double, double, double, double, int);
};

#endif
