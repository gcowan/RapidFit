// $Id: Bs2JpsiPhi_SignalAlt_MO_v2.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2JpsiPhi_SignalAlt_MO_v2 
 *
 *  Bs2JpsiPhi_SignalAlt series with mistag as observable
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-13
 */

#ifndef Bs2JpsiPhi_SignalAlt_MO_v2_H
#define Bs2JpsiPhi_SignalAlt_MO_v2_H

#include "BasePDF.h"

#include "Bs2JpsiPhi_SignalAlt_BaseClass_v2.h"
#include "RooComplex.h"

#include <exception>

class Bs2JpsiPhi_SignalAlt_MO_v2 : public BasePDF,  public Bs2JpsiPhi_SignalAlt_BaseClass_v2
{
	public:
		Bs2JpsiPhi_SignalAlt_MO_v2( PDFConfigurator* );
		~Bs2JpsiPhi_SignalAlt_MO_v2();

		//Mandatory method to evaluate the PDF value:
		virtual double EvaluateForNumericIntegral(DataPoint*) ;
		virtual double Evaluate(DataPoint*);
		virtual double EvaluateTimeOnly(DataPoint*) ;

		//Other operating methods
		virtual bool SetPhysicsParameters(ParameterSet*);
		virtual vector<string> GetDoNotIntegrateList();

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(DataPoint*, PhaseSpaceBoundary*);

	private:

		void MakePrototypes();
		
};

#endif
