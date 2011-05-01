/**
        @class LinearPDF

        An example PDF implementing IPDF directly, without BasePDF

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef LINEAR_PDF_H
#define LINEAR_PDF_H

//	RapidFit Headers
#include "BasePDF.h"

class LinearPDF : public BasePDF
{
	public:
		LinearPDF();
		LinearPDF(string, string);
		~LinearPDF();

		//Calculate the PDF value
		virtual double Evaluate(DataPoint*);

	protected:
		//Calculate the PDF normalisation
		virtual double Normalisation(PhaseSpaceBoundary*);

	private:
		//	Uncopyable!
		LinearPDF ( const LinearPDF& );
		LinearPDF& operator = ( const LinearPDF& );
		void MakePrototypes();

		string gradientName;
		string interceptName;
		string observableName;
};

#endif
