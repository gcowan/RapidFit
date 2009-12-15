// $Id: Bs2JpsiPhi.h,v 1.1 2009/11/10 10:35:49 gcowan Exp $
/** @class Bs2JpsiPhi Bs2JpsiPhi.h
 *
 *  RapidFit PDF for Bs2JpsiPhi
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-07-30
 */

#ifndef RAPDF_BS2JPSIPHI_H
#define RAPDF_BS2JPSIPHI_H

#include "BasePDF.h"

class Bs2JpsiPhi : public BasePDF
{
	public:
		Bs2JpsiPhi();
		~Bs2JpsiPhi();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);
		//virtual double Normalisation(PhaseSpaceBoundary*); // BEN_SAYS: it should look like this

	private:
		void MakePrototypes();

		// These contain the strings that correspond
		// to the physics parameter names that will be
		// used in the minimiser.
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

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF. 
		string timeName;	// proper time
		string tagName;		// B tag
		string cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		string phiName;		// azimuthal angle of the mu+ in Jpsi frame
		string cosPsiName;	// helicity angle between K+ and -ve Jpsi direction
					// in phi rest frame

		void getPhysicsParameters( double&, double&, double&, double&, double&, double&, double&, double&, double&, double&);
		void getAngularFunctions( double&, double&, double&, double&, double&, double&, DataPoint*);
		
		void getTimeDependentAmplitudes( double&, double&, double&, double&, double&, double&, DataPoint*);
		void getTimeAmplitudeIntegrals(double&, double&, double&, double&, double&, double&,  DataPoint*, PhaseSpaceBoundary*);

		double getExpInt( double, double, double ) const;
		double getExpCosInt( double, double, double, double) const;
		double getExpSinInt( double, double, double, double) const;


		double getAzeroAzeroTimeInt(double, double, double, double, double, double, double, int) const;
		double getAparaAparaTimeInt(double, double, double, double, double, double, double, int) const;
		double getAperpAperpTimeInt(double, double, double, double, double, double, double, int) const;
                
		double getAparaAperpTimeInt(double, double, double, double, double, double, double, double, double, double, int) const;
                double getAzeroAparaTimeInt(double, double, double, double, double, double, double, double, double, int) const;
                double getAzeroAperpTimeInt(double, double, double, double, double, double, double, double, double, int) const;

		double getEvenTimeComponentInt(double, double, double, double, double, double, int) const;
		double getOddTimeComponentInt (double, double, double, double, double, double, int) const;
};

#endif
