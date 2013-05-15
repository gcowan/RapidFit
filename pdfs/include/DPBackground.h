// $Id: DPBackground.h,v 1.2 2009/11/13 15:31:51 gcowan Exp $
/** @class DPBackground DPBackground.h
 *
 *  PDF for Bd2JpsiKpi dalitz background
 *
 *  @author Greig A Cowan greig.alan.cowan@cern.ch
 *  @date 2009-10-04
 */

#ifndef DPBackground_H
#define DPBackground_H

#include "BasePDF.h"

class DPBackground : public BasePDF
{
	public:
		DPBackground( PDFConfigurator* );
		DPBackground( const DPBackground & );
		~DPBackground();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

    protected:
        double Normalisation( PhaseSpaceBoundary * );

	private:
		void MakePrototypes();
		bool SetPhysicsParameters(ParameterSet*);

		// Experimental observables
		ObservableRef m23Name;
		ObservableRef cosTheta1Name;
		ObservableRef cosTheta2Name;
		ObservableRef phiName;

		double m23;
		double cosTheta1; // This is the member variable used in the "builder" functions
		double cosTheta2; // These are the physics parameters varied in the fit and passed from the XML;
		double phi;

        // Background parameterisation
        static const int l_max = 6;
        static const int i_max = 4;
        static const int k_max = 2;
        static const int j_max = 2;
        double c[l_max+1][i_max+1][k_max+1][j_max+1];
};

#endif
