# $Id: Makefile,v 1.30 2009/11/11 17:18:14 gcowan Exp $
SHELL = /bin/sh
UNAME = $(shell uname)

# Root variables
TEMPCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags | awk -F "-I" '{print $$1" -isystem"$$2}' )
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

 #GSLPATH      =/usr/local#/afs/cern.ch/sw/lcg/external/GSL/1.8/$(CMTCONFIG)
 #GSLINC       =$(GSLPATH)/include
 #GSLLIB       =$(GSLPATH)/lib

################
##linux
CXX          = g++
RM           = rm -f
AR           = ar cru

##Flags
#CXXFLAGS     += -D_DEBUG
CXXFLAGS     = -O3 -msse -msse2 -m3dnow -g -ansi -fPIC -funroll-all-loops -D__ROOFIT_NOBANNER -Wconversion -Wextra -Wsign-compare -Wfloat-equal -Wmissing-noreturn -Wall -Wno-non-virtual-dtor
#CXXFLAGS    = -Os -ansi -Wconversion -Wextra -Wsign-compare -Wfloat-equal -Wmissing-noreturn -Wall -Wno-non-virtual-dtor -fPIC

SRCEXT   = cpp
HDREXT   = h
SRCDIR   = framework/src
SRCPDFDIR= pdfs/src
INCDIR   = framework/include
INCPDFDIR= pdfs/include
OBJDIR   = framework/build
OBJPDFDIR= pdfs/build
EXEDIR   = bin
LIBDIR   = lib
UTILSSRC = utils/src
SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)' -not -name 'Roo*.cpp' | grep -v 'RapidRun' )
PDFSRCS := $(shell find $(SRCPDFDIR) -name '*.$(SRCEXT)' -not -name 'Roo*.cpp' )
HEADERS := $(shell find $(PWD)/$(INCDIR) -name '*.$(HDREXT)' | grep -v 'LinkDef' )
PDFHEAD := $(shell find $(PWD)/$(INCPDFDIR) -name '*.$(HDREXT)' )
SRCDIRS := $(shell find . -name '*.$(SRCEXT)' -exec dirname {} \; | uniq)
OBJS    := $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))
PDFOBJS := $(patsubst $(SRCPDFDIR)/%.$(SRCEXT),$(OBJPDFDIR)/%.o,$(PDFSRCS))

#ALL_HEADERS = ""
ALL_HEADERS += $(HEADERS)
ALL_HEADERS += $(PDFHEAD)

GARBAGE  = $(OBJDIR)/*.o $(OBJPDFDIR)/*.o $(EXEDIR)/fitting pdfDict.h pdfDict.cpp *.so *.rootmap

#################
##Dependencies
# Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) -I$(INCPDFDIR) $(ROOTCFLAGS) #-I$(GSLINC)
LINKFLAGS    =
endif

# OS X
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) -I$(INCPDFDIR) $(ROOTCFLAGS) #-I$(GSLINC)
LINKFLAGS    =
endif

##Libraries
#LIBS       += $(ROOTLIBS) -lHtml -lThread -lMinuit -lRooFit -lRooFitCore -lMathCore -lMinuit2 -lFoam #-lboost_thread-xgcc40-mt #-lMathMore
LIBS       += $(ROOTLIBS) -lHtml -lThread -lMinuit -lMathCore -lMinuit2 -lRooFit -lRooFitCore -lFoam #-lProof #-lboost_thread #-lMathMore

#HEADERS    =  $(INCDIR)/RooBs2PhiPhiFullPdf.h $(INCDIR)/RooPdf_Bs2JPsiPhi.h $(INCDIR)/LinkDef.h

all : $(EXEDIR)/fitting utils

$(EXEDIR)/fitting : $(OBJS) $(PDFOBJS) $(OBJDIR)/RootFileDataSet.o
	$(CXX) -o $@ $(OBJS) $(PDFOBJS) $(OBJDIR)/RootFileDataSet.o $(LINKFLAGS) $(LIBS)

# Not using RooFit PDFs anymore so we can ignore this
#$(EXEDIR)/fitting : pdfDict.o $(OBJS)
#	$(CXX) -o $@ pdfDict.o $(OBJS) $(LINKFLAGS) $(LIBS)

$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPDFDIR)/%.o : $(SRCPDFDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

#pdfDict.cpp : $(HEADERS)
#	rootcint -f $*.cpp -c $^

doc : $(OBJS) $(PDFOBJS)
	cd doc; doxygen RapidFit_doxygen.cfg; tar cvfz RapidFit_html.tgz html/; scp RapidFit_html.tgz ph-ppe:~/WWW/RapidFit/doc/; cd ..

clean   :
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE)

$(EXEDIR)/rapidfit_toyresults: $(OBJDIR)/rapidfit_toyresults.o
	$(CXX) -o $@ $< $(CXXFLAGS) $(ROOTLIBS)

$(OBJDIR)/rapidfit_toyresults.o: $(UTILSSRC)/rapidfit_toyresults.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(EXEDIR)/rapidfit_fcscanresults: $(OBJDIR)/rapidfit_fcscanresults.o
	$(CXX) -o $@ $< $(ROOTLIBS)

$(OBJDIR)/rapidfit_fcscanresults.o: $(UTILSSRC)/rapidfit_fcscanresults.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(EXEDIR)/rapidfit_fcscanresults_2: $(OBJDIR)/rapidfit_fcscanresults_2.o
	$(CXX) -o $@ $< $(ROOTLIBS)

$(OBJDIR)/rapidfit_fcscanresults_2.o: $(UTILSSRC)/rapidfit_fcscanresults_2.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#$(EXEDIR)/rapidfit_llscanresults: $(OBJDIR)/rapidfit_llscanresults.o
#	$(CXX) -o $@ $< $(ROOTLIBS)

#$(OBJDIR)/rapidfit_llscanresults.o: $(UTILSSRC)/rapidfit_llscanresults.cc
#	$(CXX) $(CXXFLAGS) -o $@ -c $<

#utils: $(EXEDIR)/rapidfit_toyresults $(EXEDIR)/rapidfit_llscanresults  $(EXEDIR)/rapidfit_fcscanresults
utils: $(EXEDIR)/rapidfit_toyresults $(EXEDIR)/rapidfit_fcscanresults $(EXEDIR)/rapidfit_fcscanresults_2


#	This command will generate a C++ file which interfaces the rest of humanity with root...
#	It requires the explicit paths of all files, or that you remain in the same working directory at all times during the build process
$(OBJDIR)/rapidfit_dict.cpp:
	@echo "Building Root Dictionary:"
	@echo "rootcint -f $(OBJDIR)/rapidfit_dict.cpp -c $(ALL_HEADERS) framework/include/LinkDef.h"
	@rootcint -f $(OBJDIR)/rapidfit_dict.cpp -c $(ALL_HEADERS) framework/include/LinkDef.h

#	Compile the class that root has generated for us which is the linker interface to root	(i.e. dictionaries & such)
$(OBJDIR)/rapidfit_dict.o: $(OBJDIR)/rapidfit_dict.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/RapidRun.o: $(SRCDIR)/RapidRun.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) -o $@ -c $<

#	Compile RapidFit as a library making use of the existing binaries for other classes
$(LIBDIR)/libRapidRun.so: $(OBJDIR)/RapidRun.o $(OBJDIR)/rapidfit_dict.o $(OBJS) $(PDFOBJS)
	$(CXX) $(LIBS) -shared -fPIC $(CXXFLAGS) $(OBJS) $(PDFOBJS) $(OBJDIR)/rapidfit_dict.o $(OBJDIR)/RapidRun.o -o $@

