// $Id: Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc Bd2JpsiPhi_mistagParameter_withTimeRes_withAverageAngAcc.h
 *
 *  RaKstarFlavourFit PDF for Bd2JpsiKstar
 *
 *  @author Ailsa Sparkes
 *  @date 2009-07-30
 */

#ifndef Bd2JpsiKstar_withTimeRes_withAverageAngAcc_H
#define Bd2JpsiKstar_withTimeRes_withAverageAngAcc_H

#include "BasePDF.h"

class Bd2JpsiKstar_withTimeRes_withAverageAngAcc : public BasePDF
{
	public:
		Bd2JpsiKstar_withTimeRes_withAverageAngAcc( PDFConfigurator* );
		~Bd2JpsiKstar_withTimeRes_withAverageAngAcc();

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
		//double q;
		// These contain the strings that correspond
		// to the physics parameter names that will be
		// used in the minimiser.
		ObservableRef gammaName;		// gamma
		ObservableRef deltaMName;		// delta mass
		ObservableRef Azero_sqName;	// amplitude
		ObservableRef Apara_sqName;	// amplitude
		ObservableRef Aperp_sqName;	// amplitude
		ObservableRef delta_zeroName;	// strong phase, set to 0
		ObservableRef delta_paraName;	// strong phase
		ObservableRef delta_perpName;	// strong phase
		ObservableRef angAccI1Name;		// Pre-calculated angular integrals including acceptance
		ObservableRef angAccI2Name;		//
		ObservableRef angAccI3Name;		//
		ObservableRef angAccI4Name;		//
		ObservableRef angAccI5Name;		//
		ObservableRef angAccI6Name;		//
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
		ObservableRef KstarFlavourName;
		//kaonIDName;				// in phi rest frame

		ObservableRef timeconstraintName;

		// Member variables that will contain the parameter values
                double gamma;
                double deltaMs;
                double Azero_sq;
                double Apara_sq;
                double Aperp_sq;
                double delta_zero;
                double delta_para;
                double delta_perp;
		double AzeroApara;
                double AzeroAperp;
                double AparaAperp;
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


};


#endif
