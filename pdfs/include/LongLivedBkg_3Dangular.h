// $Id: LongLivedBkg_3Dangular.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class LongLivedBkg_3Dangular LongLivedBkg_3Dangular_withAngDist.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution + non-trivial angular distribution realised by a histogram
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
		//LongLivedBkg_3Dangular();
		LongLivedBkg_3Dangular(PDFConfigurator);
		~LongLivedBkg_3Dangular();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:
		//	Can't be copied
		LongLivedBkg_3Dangular ( const LongLivedBkg_3Dangular& );
		LongLivedBkg_3Dangular& operator = ( const LongLivedBkg_3Dangular& );		

		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();
		double angularFactor( );

		// Physics parameters
		pair<string,int> tauLL1Name;		// decay constant 1
		pair<string,int> tauLL2Name;		// decay constant 2
		pair<string,int> f_LL1Name;		// fraction of decay const 1
	
	    //Detector parameters
		pair<string,int> sigmaLL1Name;	// time res sigma 1
		pair<string,int> sigmaLL2Name;	// time res sigma 2
		pair<string,int> timeResLL1FracName; //fraction of timeres 1

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF.
		pair<string,int> timeName;		// proper time
		pair<string,int> cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
		pair<string,int> phiName;			// azimuthal angle of the mu+ in Jpsi frame
		pair<string,int> cosPsiName;		// helicity angle between K+ and -ve Jpsi direction

		double tauLL1;
		double tauLL2;
		double f_LL1;
		double sigmaLL; 
		double sigmaLL1; 
		double sigmaLL2;
		double timeResLL1Frac;

		double tlow, thigh; // time integration limits

		double time;
		double cosTheta;
		double phi;
		double cosPsi;

        //Additions to deal with 3-D angular distribution via a histogram
		bool useFlatAngularDistribution;
		TH3D *histo;
		TAxis *xaxis, *yaxis, *zaxis;
		int nxbins, nybins, nzbins;
		double xmin, xmax, ymin, ymax, zmin, zmax, deltax, deltay, deltaz;
		double total_num_entries;
};

#endif
