
///	ROOT Headers
#include "TMatrixDSym.h"
///	RapidFit Headers
#include "IMinimiser.h"

using namespace::std;

class CorrectedCovariance
{
	public:
		static TMatrixDSym* GetCorrectedCovarianceMatrix( IMinimiser* thisMinimiser );
};

