/**
        @class RooJPsiPhiWrapper

        A wrapper for a RooFit format PDF describing Bs decay to J/Psi Phi

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef ROO_JPSIPHI_WRAPPER_H
#define ROO_JPSIPHI_WRAPPER_H

#include "RooPDFWrapper.h"

class RooJPsiPhiWrapper : public RooPDFWrapper
{
	public:
		RooJPsiPhiWrapper();
		RooJPsiPhiWrapper( vector<string>, vector<string> );
		~RooJPsiPhiWrapper();
};

#endif
