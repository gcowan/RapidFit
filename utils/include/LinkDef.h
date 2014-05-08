//	Voodoo required for compiling RapidFit as a CINT library
//	It should now again work :D

//-------------------------------------------- LinkDef.h
#pragma once
#ifndef _LINKERDEF_RAPIDFIT_UTIL
#define _LINKERDEF_RAPIDFIT_UTIL
#ifdef __CINT__

//	Can be defined globally or just for CINT, no harm either way
#pragma link off all class;
#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;

//	This is the Class that I want to be able to access from CINT
#pragma link C++ class ROOT_File_Processing+;
#pragma link C++ class TTree_Processing+;
#pragma link C++ class Histogram_Processing+;
#pragma link C++ class RapidFit_Output_File+;
#pragma link C++ class StringOperations+;
#pragma link C++ class Template_Functions+;
#pragma link C++ class XMLUtilFunctions+;
#pragma link C++ class Mathematics+;
#pragma link C++ class RapidLL+;
#pragma link C++ class Rapid2DLL+;
#pragma link C++ class ToyStudyAnalysis+;

#endif

#endif

