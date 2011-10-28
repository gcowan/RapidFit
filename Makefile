# $Id: Makefile,v 1.30 2009/11/11 17:18:14 gcowan Exp $
SHELL = /bin/bash
UNAME = $(shell uname)

#		ROOT
TEMPCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)

#		Include Root files as system headers as they're NOT standards complient and we do not want to waste time fixing them!
#		ROOT has some broken backwards compatability for OSX so won't claim to be a set of system headers
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags | awk -F "-I" '{print $$1" -isystem"$$2}' )
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)
EXTRA_ROOTLIBS=-lHtml -lThread -lMinuit -lMinuit2 -lRooFit -lRooStats -lRooFitCore -lFoam -lMathMore

#		Command Line Tools
CXX          = g++
RM           = rm -f

SVN_REV = $(shell svnversion -n .)

#		Compiler Flags
CXXFLAGS     = -DSVN_REV=$(SVN_REV) -fPIC -Wabi -Weffc++ -O3 -msse -msse2 -msse3 -m3dnow -g -ansi -fmerge-all-constants -funroll-all-loops -D__ROOFIT_NOBANNER -Wconversion -Wextra -Wsign-compare -Wfloat-equal -Wmissing-noreturn -Wall -Wno-non-virtual-dtor -Wno-reorder -pthread

#		Some Useful global variables, makes this file MUCH easier to maintain
SRCEXT   = cpp
HDREXT   = h
SRCDIR   = framework/src
SRCPDFDIR= pdfs/src
UTILSSRC  = utils/src
INCDIR   = framework/include
INCPDFDIR= pdfs/include
INCUITLS = utils/include
OBJDIR   = framework/build
OBJPDFDIR= pdfs/build
EXEDIR   = bin
LIBDIR   = lib


#	Source Files to be Built	ignoring all files in 'unused' and the RapidRun source for ROOT linking
SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)' | grep -v 'unused' )
#	PDF source files which will be required for building libpdf.so
PDFSRCS := $(shell find $(SRCPDFDIR) -name '*.$(SRCEXT)' | grep -v 'unused' )


#	Absolute Paths of headers	ignoring the LinkDef written for ROOT and ignoring unused code
HEADERS := $(shell find $(PWD)/$(INCDIR) -name '*.$(HDREXT)' | grep -v 'unused' | grep -v 'LinkDef' )
PDFHEAD := $(shell find $(PWD)/$(INCPDFDIR) -name '*.$(HDREXT)' )


#	Binary Objects to make in the build
OBJS    := $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))
#	Binary objects to be linked into libpdf.so
PDFOBJS := $(patsubst $(SRCPDFDIR)/%.$(SRCEXT),$(OBJPDFDIR)/%.o,$(PDFSRCS))



#	All Headers in their absolute paths as required by CINT to construct a dictionary of RapidFit
ALL_HEADERS += $(HEADERS)
ALL_HEADERS += $(PDFHEAD)

UTIL_HEADERS = $(shell find $(PWD)/$(INCUITLDIR) -name '*.$(HDREXT)' )

#	BUILD OUTPUT
OUTPUT  = $(OBJDIR)/*.o $(OBJPDFDIR)/*.o $(EXEDIR)/fitting $(LIBDIR)/*.so $(OBJDIR)/rapidfit_dict.* *.so *.rootmap $(EXEDIR)/rapidfit_toyresults $(EXEDIR)/rapidfit_fcscanresults $(EXEDIR)/rapidfit_fcscanresults_2 $(EXEDIR)/betas_sweightfitter $(EXEDIR)/merge_plot $(EXEDIR)/RapidLL $(EXEDIR)/RapidPlot $(EXEDIR)/print



#################
##Dependencies

LINKER=ld
LINKFLAGS=-lpthread

CXXFLAGS     += -I$(INCDIR) -I$(INCPDFDIR) $(ROOTCFLAGS)
CXXFLAGS_LIB += -I$(INCDIR) -I$(INCPDFDIR) $(ROOTCFLAGS)

# Linux
ifeq "$(UNAME)" "Linux"
	CXXFLAGS     += -fPIE
	LINKFLAGS    += -Wl,--export-dynamic -pie
endif

# OS X
ifeq "$(UNAME)" "Darwin"
	
endif





#	Default build command when someone asks for 'make'
all : $(EXEDIR)/fitting utils lib 


$(OBJPDFDIR)/%.o : $(SRCPDFDIR)/%.$(SRCEXT) $(INCPDFDIR)/%.$(HDREXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

#	Some fairly cool Makefile code for linking the build together :D
$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(INCDIR)/%.$(HDREXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

#	Main Build of RapidFit Binary
$(EXEDIR)/fitting : $(OBJS) $(PDFOBJS) $(OBJDIR)/rapidfit_dict.o
	$(CXX) $(LINKFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(EXTRA_ROOTLIBS)







#	Does anyone use this any more?
doc : $(OBJS) $(PDFOBJS)
	cd doc; doxygen RapidFit_doxygen.cfg; tar cvfz RapidFit_html.tgz html/; scp RapidFit_html.tgz ph-ppe:~/WWW/RapidFit/doc/; cd ..



#	Cleanup
clean   :	distclean
cleanall:	distclean
distclean:
	$(RM) $(OUTPUT)


clang: override CXX=clang++
clang: all





#	RapidFit Utils:

#	Tool for plotting Toy Study Results
$(EXEDIR)/rapidfit_toyresults: $(OBJDIR)/rapidfit_toyresults.o
	$(CXX) -o $@ $< $(LINKFLAGS) $(ROOTLIBS)
$(OBJDIR)/rapidfit_toyresults.o: $(UTILSSRC)/rapidfit_toyresults.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<


#	Original Tool for plotting FCScan Results
$(EXEDIR)/rapidfit_fcscanresults: $(OBJDIR)/rapidfit_fcscanresults.o
	$(CXX) -o $@ $< $(LINKFLAGS) $(ROOTLIBS)
$(OBJDIR)/rapidfit_fcscanresults.o: $(UTILSSRC)/rapidfit_fcscanresults.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<


#	Tool for plotting FCScan Results and more
$(EXEDIR)/rapidfit_fcscanresults_2: $(OBJDIR)/rapidfit_fcscanresults_2.o
	$(CXX) -o $@ $< $(LINKFLAGS) $(ROOTLIBS)
$(OBJDIR)/rapidfit_fcscanresults_2.o: $(UTILSSRC)/rapidfit_fcscanresults_2.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#	Tool that sWeights the betas ntuple
$(EXEDIR)/betas_sweightfitter: $(OBJDIR)/betas_sweightfitter.o
	$(CXX) -o $@ $< $(LINKFLAGS) $(ROOTLIBS) $(EXTRA_ROOTLIBS)
$(OBJDIR)/betas_sweightfitter.o: $(UTILSSRC)/betas_sweightfitter.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#	Trying to move people off this tool
#$(EXEDIR)/rapidfit_llscanresults: $(OBJDIR)/rapidfit_llscanresults.o
#	$(CXX) -o $@ $< $(ROOTLIBS)
#$(OBJDIR)/rapidfit_llscanresults.o: $(UTILSSRC)/rapidfit_llscanresults.cc
#	$(CXX) $(CXXFLAGS) -o $@ -c $<

#       A VERY stupid tool designed to overlay the results of 2 LL scans from within Edinburgh
$(EXEDIR)/merge_plot: $(OBJDIR)/merge_plot.o
	$(CXX) -o $@ $< $(LINKFLAGS) $(ROOTLIBS)
$(OBJDIR)/merge_plot.o: $(UTILSSRC)/merge_plot.C
	$(CXX) $(CXXFLAGS) -o $@ -c $<


#	Binaries sharing LOTS of useful code for processing/outputting results

$(OBJDIR)/NTuple_Processing.o: $(UTILSSRC)/NTuple_Processing.C
	$(CXX) $(CXXFLAGS) -I$(INCUITLS) -o $@ -c $<
$(OBJDIR)/TString_Processing.o: $(UTILSSRC)/TString_Processing.C
	$(CXX) $(CXXFLAGS) -I$(INCUITLS) -o $@ -c $<
$(OBJDIR)/Histo_Processing.o: $(UTILSSRC)/Histo_Processing.C
	$(CXX) $(CXXFLAGS) -I$(INCUITLS) -o $@ -c $<

#       New tool for plotting 1D LL and overlaying multiple copies of the same
$(EXEDIR)/RapidLL: $(OBJDIR)/RapidLL.o $(OBJDIR)/EdStyle.o $(OBJDIR)/NTuple_Processing.o $(OBJDIR)/TString_Processing.o $(OBJDIR)/StringProcessing.o
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)
$(OBJDIR)/RapidLL.o: $(UTILSSRC)/RapidLL.C
	$(CXX) $(CXXFLAGS) -Iutils/include -o $@ -c $<

#       New tool for plotting 1D LL and overlaying multiple copies of the same
$(EXEDIR)/RapidLL2: $(OBJDIR)/RapidLL2.o $(OBJDIR)/EdStyle.o $(OBJDIR)/NTuple_Processing.o $(OBJDIR)/TString_Processing.o $(OBJDIR)/StringProcessing.o
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)
$(OBJDIR)/RapidLL2.o: $(UTILSSRC)/RapidLL2.C
	$(CXX) $(CXXFLAGS) -Iutils/include -o $@ -c $<

#       New tool for plotting 2DLL and FC
$(EXEDIR)/RapidPlot: $(OBJDIR)/RapidPlot.o $(OBJDIR)/EdStyle.o $(OBJDIR)/NTuple_Processing.o $(OBJDIR)/Histo_Processing.o $(OBJDIR)/TString_Processing.o $(OBJDIR)/StringProcessing.o
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)
$(OBJDIR)/RapidPlot.o: $(UTILSSRC)/RapidPlot.C
	$(CXX) $(CXXFLAGS) -Iutils/include -o $@ -c $<

#       New tool for plotting 2DLL and FC
$(EXEDIR)/print: $(OBJDIR)/print.o $(OBJDIR)/EdStyle.o $(OBJDIR)/NTuple_Processing.o $(OBJDIR)/Histo_Processing.o $(OBJDIR)/TString_Processing.o $(OBJDIR)/StringProcessing.o
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)
$(OBJDIR)/print.o: $(UTILSSRC)/print.C
	$(CXX) $(CXXFLAGS) -Iutils/include -o $@ -c $<

#	Tinter tool for analysing old format toy stuides
$(EXEDIR)/tinter: $(OBJDIR)/tinter.o $(OBJDIR)/EdStyle.o $(OBJDIR)/TString_Processing.o $(OBJDIR)/Histo_Processing.o $(OBJDIR)/NTuple_Processing.o $(OBJDIR)/StringProcessing.o
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)
$(OBJDIR)/tinter.o: $(UTILSSRC)/tinter.C
	$(CXX) $(CXXFLAGS) -Iutils/include -o $@ -c $<


utils: $(EXEDIR)/rapidfit_toyresults $(EXEDIR)/rapidfit_fcscanresults $(EXEDIR)/rapidfit_fcscanresults_2 $(EXEDIR)/betas_sweightfitter $(EXEDIR)/print $(EXEDIR)/merge_plot $(EXEDIR)/RapidLL $(EXEDIR)/RapidPlot






#	For building RapidFit as a library to use within CINT which makes life easier on the grid... (supposedly)
#	make lib

lib:    $(LIBDIR)/libRapidRun.so


#	This command will generate a C++ file which interfaces the rest of humanity with root...
#	It requires the explicit paths of all files, or that you remain in the same working directory at all times during the build process
#	We want to place the output dictionary in the Build directory as this is CODE that is NOT to be editted by the $USER!
$(OBJDIR)/rapidfit_dict.cpp: $(ALL_HEADERS) $(PWD)/framework/include/LinkDef.h
	@echo "Building Root Dictionary:"
	#@echo "rootcint -f $(OBJDIR)/rapidfit_dict.cpp -c -I$(PWD)/framework/include $^"
	@rootcint -f $(OBJDIR)/rapidfit_dict.cpp -c -I$(PWD)/framework/include $^

#	Compile the class that root has generated for us which is the linker interface to root	(i.e. dictionaries & such)
$(OBJDIR)/rapidfit_dict.o: $(OBJDIR)/rapidfit_dict.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#	Class which has a dictionary generated for it, think of this as the equivalent to int main() in a CINT-y Universe
$(OBJDIR)/RapidRun.o: $(SRCDIR)/RapidRun.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#	Finally, Compile RapidFit as a library making use of the existing binaries for other classes
$(LIBDIR)/libRapidRun.so: $(OBJDIR)/RapidRun.o $(OBJDIR)/rapidfit_dict.o $(OBJS) $(PDFOBJS) $(OBJDIR)/ClassLookUp.o
	$(CXX) -shared $(LIBS) $(LINKFLAGS) $(ROOTLIBS) $(EXTRA_ROOTLIBS) $(CXXFLAG) $^ -o $@

