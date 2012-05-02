// $Id: Bs2Jpsifzero_SignalAlt_MO_dev.h,v 1.1 2009/12/06  Pete Clarke Exp $
/** @class Bs2Jpsifzero_SignalAlt_MO_dev 
 *
 *  Bs2Jpsifzero_SignalAlt series with mistag as observable
 *
 *  @author Pete Clarke peter.clarke@ed.ac.uk
 *  @date 2011-02-13
 */

#ifndef Bs2Jpsifzero_SignalAlt_MO_dev_H
#define Bs2Jpsifzero_SignalAlt_MO_dev_H

#include "BasePDF.h"

#include "Bs2Jpsifzero_SignalAlt_BaseClass_dev.h"

#include <exception>

class Bs2Jpsifzero_SignalAlt_MO_dev : public BasePDF,  public Bs2Jpsifzero_SignalAlt_BaseClass_dev
{
	public:
		Bs2Jpsifzero_SignalAlt_MO_dev( PDFConfigurator* );
		~Bs2Jpsifzero_SignalAlt_MO_dev();

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
