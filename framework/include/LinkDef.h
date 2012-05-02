//	Voodoo required for compiling RapidFit as a CINT library
//	It should now again work :D

#pragma once
#ifdef __CINT__

//	Can be defined globally or just for CINT, no harm either way
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

//	This is the Class that I want to be able to access from CINT
#pragma link C++ class RapidRun;


//	Are these required
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ global gROOT;
#pragma link C++ global gEnv;

#endif

