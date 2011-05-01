/**
        @class RooGaussianWrapper

        A wrapper for the RooFit gaussian PDF

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#ifndef ROO_GAUSSIAN_WRAPPER_H
#define ROO_GAUSSIAN_WRAPPER_H

#include "RooPDFWrapper.h"

class RooGaussianWrapper : public RooPDFWrapper
{
	public:
		RooGaussianWrapper();
		RooGaussianWrapper( vector<string>, vector<string> );
		~RooGaussianWrapper();
};

#endif
