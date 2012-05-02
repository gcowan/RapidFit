
#include "TFile.h"

#include <vector>
#include <string>

#include "Template_Functions.h"

using namespace::std;

int Component_Projections( TFile* input_file, vector<pair<string,string> > Directories )
{
	print( Directories );

	(void) input_file;

	return 0;
}

