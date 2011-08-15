// $Id: LongLivedBkg.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class LongLivedBkg LongLivedBkg_withAngDist.h
 *
 *  PDF for Bs2JpsiPhi long lived background with time resolution + non-trivial angular distribution realised by a histogram
 *
 *  @author Ailsa Sparkes
 *  @date 2011-05-30
 */

#ifndef LongLivedBkg_H
#define LongLivedBkg_H

#include "TROOT.h"
#include "TFile.h"
#include "TH3D.h"

#ifndef __CINT__
#include "BasePDF.h"
#include "SlicedAcceptance.h"
#endif
#ifdef __CINT__
#include "framework/include/BasePDF.h"
#include "framework/include/SlicedAcceptance.h"
#endif

class LongLivedBkg : public BasePDF
{
	public:
		//LongLivedBkg();
		LongLivedBkg(PDFConfigurator);
		~LongLivedBkg();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);
		virtual double Norm(DataPoint*, PhaseSpaceBoundary*);

	private:
		//	Can't be copied
		LongLivedBkg ( const LongLivedBkg& );
		LongLivedBkg& operator = ( const LongLivedBkg& );		

		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);
		double buildPDFnumerator();
		double buildPDFdenominator();
		double angularFactor( );

		// Physics parameters
		ObservableRef f_LL1Name;		// fraction of decay const 1
                ObservableRef tauLL1Name;            // decay constant 1
		ObservableRef tauLL2Name;            // decay constant 2

		//Detector parameters
		ObservableRef timeResLL1FracName; //fraction of timeres 1
		ObservableRef sigmaLL1Name;	// time res sigma 1
		ObservableRef sigmaLL2Name;	// time res sigma 2

		// These contain the strings that correspond
		// to the observable names that are used in the
		// PDF.
		ObservableRef timeName;		// proper time
//		ObservableRef cosThetaName;	// cos of angle of mu+ wrt z-axis in Jpsi frame
//		ObservableRef phiName;			// azimuthal angle of the mu+ in Jpsi frame
//		ObservableRef cosPsiName;		// helicity angle between K+ and -ve Jpsi direction

		ObservableRef timeconstName;

		double tauLL1;
		double tauLL2;
		double f_LL1;
		double sigmaLL; 
		double sigmaLL1; 
		double sigmaLL2;
		double timeResLL1Frac;

		double tlow, thigh; // time integration limits

		double time;
	//	double cosTheta;
	//	double phi;
	//	double cosPsi;

		//Additions to deal with 3-D angular distribution via a histogram
		TH3D *histo;
		TAxis *xaxis, *yaxis, *zaxis;
		int nxbins, nybins, nzbins;
		double xmin, xmax, ymin, ymax, zmin, zmax, deltax, deltay, deltaz;
		double total_num_entries;
		bool useFlatAngularDistribution;
                bool _useTimeAcceptance;
                SlicedAcceptance * timeAcc;
};

#endif
