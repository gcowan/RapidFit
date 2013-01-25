
#pragma once

#ifndef _CORRMATRIX_HDR
#define _CORRMATRIX_HDR

#include "TTree.h"

#include <vector>

using namespace::std;

class CorrMatrix
{

	public:

	static void Analyse( const vector<TTree*> corr_trees, const vector<string> other_params );

	static void Help();
};

#endif

