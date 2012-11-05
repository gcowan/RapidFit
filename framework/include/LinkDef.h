//	Voodoo required for compiling RapidFit as a CINT library
//	It should now again work :D

//-------------------------------------------- LinkDef.h
#pragma once
#ifndef _LINKERDEF_RAPIDFIT
#define _LINKERDEF_RAPIDFIT
#ifdef __CINT__

//	Can be defined globally or just for CINT, no harm either way
#pragma link off all class;
#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;



//	This is the Class that I want to be able to access from CINT
#pragma link C++ class RapidRun+;


//	Are these required
//#pragma link C++ nestedclasses;
//#pragma link C++ nestedtypedefs;
//#pragma link C++ global gROOT;
//#pragma link C++ global gEnv;

#endif

#endif

