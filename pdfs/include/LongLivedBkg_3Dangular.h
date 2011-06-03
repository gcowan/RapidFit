// $Id: LongLivedBkg_3Dangular.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class LongLivedBkg_3Dangular LongLivedBkg_3Dangular_withAngDist.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution + non-trivial angular distribution
 *
 *  @author Ailsa Sparkes
 *  @date 2011-05-30
 */

#ifndef LongLivedBkg_3Dangular_H
#define LongLivedBkg_3Dangular_H

#include "TROOT.h"
#include "TFile.h"
#include "TH3D.h"

#ifndef __CINT__
#include "BasePDF.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#endif

class LongLivedBkg_3Dangular : public BasePDF
{
	public:
		LongLivedBkg_3Dangular();
		~LongLivedBkg_3Dangular();

		//Calculate the PDF value

		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();
		double angularFactor(double, double, double);

		// Physics parameters
		string tauLL1Name;		// decay constant 1
		string tauLL2Name;		// decay constant 2
		string f_LL1Name;		// fraction
		string sigmaLL1Name;		// time res sigma 1
		string sigmaLL2Name;		// time res sigma 2
                string timeResLL1FracName;
//		string f_JpsiName;
//		string f_NoJpsiName;

		double tauLL1;
		double tauLL2;
		double f_LL1;
		double sigmaLL; // This is the member variable used in the "builder" functions
		double sigmaLL1; // These are the physics parameters varied in the fit and passed from the XML;
		double sigmaLL2;
                double timeResLL1Frac;
//		double f_Jpsi;
//		double f_NoJpsi;
		double tlow, thigh; // integration limits

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF.
		string timeName;	// proper time
		string cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		string phiName;		// azimuthal angle of the mu+ in Jpsi frame
		string cosPsiName;		// helicity angle between K+ and -ve Jpsi direction
		double time;
		double cosTheta;
		double phi;
		double cosPsi;


		TH3D *histo;
       		TFile* f;
		TAxis *xaxis, *yaxis, *zaxis;
		int nxbins, nybins, nzbins;
		double xmin, xmax, ymin, ymax, zmin, zmax, deltax, deltay, deltaz;
		double total_num_entries;
};

#endif
