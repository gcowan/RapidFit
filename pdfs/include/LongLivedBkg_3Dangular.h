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

#include "IResolutionModel.h"
#include "BasePDF.h"

#include <vector>
#include <string>

using namespace::std;

class LongLivedBkg_3Dangular : public BasePDF
{
	public:
		//LongLivedBkg_3Dangular();
		LongLivedBkg_3Dangular( PDFConfigurator* );
		//LongLivedBkg_3Dangular( const LongLivedBkg_3Dangular& );
		~LongLivedBkg_3Dangular();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		//LongLivedBkg_3Dangular& operator=( const LongLivedBkg_3Dangular& );
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);

		double angularFactor( DataPoint* );

		bool _useTimeAcceptance;

		// Physics parameters
		ObservableRef f_LL1Name;		// fraction of decay const 1
		ObservableRef tauLL1Name;            // decay constant 1
		ObservableRef tauLL2Name;            // decay constant 2

		// For transversity angles
		ObservableRef cosThetaName;	
		ObservableRef phiName;			
		ObservableRef cosPsiName;	
		// For helicity angles
		ObservableRef cthetakName;	
		ObservableRef cthetalName;		
		ObservableRef phihName;			

		ObservableRef timeName;
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
		//double cosTheta;
		//double phi;
		//double cosPsi;
		//double ctheta_k;
		//double ctheta_l;
		//double phi_h;
		double cos1;
		double cos2;
		double phi;
		bool _useHelicityBasis ;
		double useHelicityBasis() { return _useHelicityBasis; }


		TH1* histo;
		TAxis* xaxis;
		TAxis* yaxis;
		TAxis* zaxis;
		//Additions to deal with 3-D angular distribution via a histogram
		int nxbins, nybins, nzbins;
		double xmin, xmax, ymin, ymax, zmin, zmax, deltax, deltay, deltaz;
		double total_num_entries;
		bool useFlatAngularDistribution;

		DataPoint* _datapoint;

		IResolutionModel* resolutionModel;
};

#endif

