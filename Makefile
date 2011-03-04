# $Id: Makefile,v 1.30 2009/11/11 17:18:14 gcowan Exp $
SHELL = /bin/sh
UNAME = $(shell uname)

# Root variables
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
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
CXXFLAGS     = -O3 -g -fPIC -funroll-loops -D__ROOFIT_NOBANNER  -Wall -Wno-non-virtual-dtor

SRCEXT   = cpp
SRCDIR   = framework/src
SRCPDFDIR= pdfs/src
INCDIR   = framework/include
INCPDFDIR= pdfs/include
OBJDIR   = framework/build
OBJPDFDIR= pdfs/build
EXEDIR   = bin
UTILSSRC = utils/src
SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)' -not -name 'Roo*.cpp')
PDFSRCS := $(shell find $(SRCPDFDIR) -name '*.$(SRCEXT)' -not -name 'Roo*.cpp')
SRCDIRS := $(shell find . -name '*.$(SRCEXT)' -exec dirname {} \; | uniq)
OBJS    := $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))
PDFOBJS := $(patsubst $(SRCPDFDIR)/%.$(SRCEXT),$(OBJPDFDIR)/%.o,$(PDFSRCS))

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

HEADERS    =  $(INCDIR)/RooBs2PhiPhiFullPdf.h $(INCDIR)/RooPdf_Bs2JPsiPhi.h $(INCDIR)/LinkDef.h

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
	$(CXX) -o $@ $< $(ROOTLIBS)

$(OBJDIR)/rapidfit_toyresults.o: $(UTILSSRC)/rapidfit_toyresults.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(EXEDIR)/rapidfit_llscanresults: $(OBJDIR)/rapidfit_llscanresults.o
	$(CXX) -o $@ $< $(ROOTLIBS)

$(OBJDIR)/rapidfit_llscanresults.o: $(UTILSSRC)/rapidfit_llscanresults.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

utils: $(EXEDIR)/rapidfit_toyresults $(EXEDIR)/rapidfit_llscanresults
	
