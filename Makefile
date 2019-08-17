# $Id: Makefile,v 1.30 2009/11/11 17:18:14 gcowan Exp $
SHELL=/bin/bash
UNAME=$(shell uname -s )
CC=g++



#		ROOT
TEMPCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)

#		Include Root files as system headers as they're NOT standards complient and we do not want to waste time fixing them!
#		ROOT has some broken backwards compatability for OSX so won't claim to be a set of system headers
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags | awk -F "-I" '{print $$1" -isystem"$$2}' )
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

#               On some Systems with Mathmore compiled, sometimes things need to be resolved against it... I don't know why
EXTRA_ROOTLIBS=-lTreePlayer -lHtml -lThread -lMinuit -lMinuit2 -lRooFit -lRooStats -lRooFitCore -lFoam $(shell if [ "$(shell root-config --features | grep mathmore)" == "" ]; then echo "" ; else echo "-lMathMore" ; fi)

#		Command Line Tools
CXX          = $(CC) $(shell if [ "$(shell root-config --arch | grep 32)" = "" ]; then echo ""; else echo "--arch=i386"; fi)
RM           = rm -f

#	We can ignore 99% of Makefile modifications it builds and links or it doesn't normally
#	We can ignore pdf modifications, 95% of all work is exploring the effect of changing a pdf so you don't care about benchmarking the version I claim
SVN_REV ="$(shell svnversion -n ./framework)"
SVN_PDF_REV ="$(shell svnversion -n ./pdfs)"
BUILD_DATE ="$(shell date +%H:%M_%F)"

CXXFLAGS_BASE_MINIMAL = -DSVN_REV=$(SVN_REV) -DSVN_PDF_REV=$(SVN_PDF_REV) -DBUILD_DATE=$(BUILD_DATE) -rdynamic -D_GNU_SOURCE -D__USE_GNU -fPIC -pthread

CXXFLAGS_BASE_WARNINGS = -Wconversion -Wextra -Wsign-compare -Wfloat-equal -Wall -Wno-non-virtual-dtor -Wno-reorder -Wshadow -Wmissing-noreturn -Wcast-align

#		Compiler Flags
CXXFLAGS_BASE_COMMON  = $(CXXFLAGS_BASE_MINIMAL) -D__ROOFIT_NOBANNER  $(CXXFLAGS_BASE_WARNINGS)

CXXFLAGS_BASE_OPT = -O3 -msse2 -msse3 -fmerge-all-constants -funroll-all-loops -fno-common -m3dnow

#CXXFLAGS_BASE = $(CXXFLAGS_BASE_COMMON) -Wmissing-noreturn -Wcast-align -msse -m3dnow

CXXFLAGS_BASE = -std=c++11 $(CXXFLAGS_BASE_COMMON) -O3 -msse2 -msse3 -m3dnow -ftree-vectorize -finline-limit=2000 -fprefetch-loop-arrays -fmerge-all-constants $(CXXFLAGS_BASE_WARNINGS)

CXX_FLAGS_LITE = -DSVN_REV=$(SVN_REV) -DSVN_PDF_REV=$(SVN_PDF_REV) -DBUILD_DATE=$(BUILD_DATE) -rdynamic -D_GNU_SOURCE -D__USE_GNU -fPIC -O3 -msse -msse2 -msse3 -m3dnow -ansi -fmerge-all-constants -funroll-all-loops -fno-common -D__ROOFIT_NOBANNER -Wconversion -Wextra -Wsign-compare -Wfloat-equal -Wmissing-noreturn -Wall -Wno-non-virtual-dtor -Wno-reorder -pthread -Wshadow -Wcast-align

#		Some Useful global variables, makes this file MUCH easier to maintain
SRCEXT    = cpp
HDREXT    = h
UTILSSRCEXT=C
SRCDIR    = framework/src
SRCPDFDIR = pdfs/src
UTILSSRC  = utils/src
INCDIR    = framework/include
INCPDFDIR = pdfs/include
INCUTILS  = utils/include
INCGSL_1  = $(shell if command -v gsl-config >/dev/null 2>&1; then echo "$(shell gsl-config --cflags)"; else echo ""; fi )
LINKGSL_1 = $(shell if command -v gsl-config >/dev/null 2>&1; then echo "$(shell gsl-config --libs)"; else echo ""; fi )
USE_GSL_1 = $(shell if command -v gsl-config >/dev/null 2>&1; then echo "-D__RAPIDFIT_USE_GSL"; else echo ""; fi )
INCGSL    = $(shell if test -d /sw/lib/lcg/external/GSL/1.10/x86_64-slc5-gcc43-opt/include; then echo "-I/sw/lib/lcg/external/GSL/1.10/x86_64-slc5-gcc43-opt/include"; else echo ${INCGSL_1}; fi )
LINKGSL   = $(shell if test -d /sw/lib/lcg/external/GSL/1.10/x86_64-slc5-gcc43-opt/include; then echo "-L/sw/lib/lcg/external/GSL/1.10/x86_64-slc5-gcc43-opt/lib -lgsl -lgslcblas -lm"; else echo ${LINKGSL_1}; fi )
USE_GSL   = $(shell if test -d /sw/lib/lcg/external/GSL/1.10/x86_64-slc5-gcc43-opt/include; then echo "-D__RAPIDFIT_USE_GSL"; else echo ${USE_GSL_1}; fi )
OBJDIR    = framework/build
OBJPDFDIR = pdfs/build
OBJUTILDIR= utils/build
EXEDIR    = bin
LIBDIR    = lib
SRCDALITZEXT = cc
HDRDALITZEXT = hh
SRCDALITZDIR = pdfs/dalitz/src
INCDALITZDIR = pdfs/dalitz/include
OBJDALITZDIR = pdfs/dalitz/build


#	Source Files to be Built	ignoring all files in 'unused' and the RapidRun source for ROOT linking
SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)' | grep -v 'unused' )
#UTILSSRCS := $(shell find $(UTILSSRC) -name '*.$(UTILSSRCEXT)' )
#	PDF source files which will be required for building libpdf.so
PDFSRCS := $(shell find $(SRCPDFDIR) -name '*.$(SRCEXT)' | grep -v 'unused' )
DALITZSRCS := $(shell find $(SRCDALITZDIR) -name '*.$(SRCDALITZEXT)' | grep -v 'unused' )


#	Absolute Paths of headers	ignoring the LinkDef written for ROOT and ignoring unused code
HEADERS := $(shell find $(INCDIR) -name '*.$(HDREXT)' | grep -v 'unused' | grep -v 'LinkDef' )
UTILHEADERS := $(shell find $(INCUTILS) -name '*.$(HDREXT)' | grep -v 'unused' | grep -v 'LinkDef' )
PDFHEAD := $(shell find $(INCPDFDIR) -name '*.$(HDREXT)' )
DALITZHEAD := $(shell find $(INCDALITZDIR) -name '*.$(HDRDALITZEXT)' )


#	Binary Objects to make in the build
OBJS    := $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))
#UTILOBJS := $(patsubst $(UTILSSRC)/%.$(UTILSSRCEXT),$(OBJUTILDIR)/%.o,$(UTILSSRCS))
#	Binary objects to be linked into libpdf.so
PDFOBJS := $(patsubst $(SRCPDFDIR)/%.$(SRCEXT),$(OBJPDFDIR)/%.o,$(PDFSRCS))
DALITZOBJS := $(patsubst $(SRCDALITZDIR)/%.$(SRCDALITZEXT),$(OBJDALITZDIR)/%.o,$(DALITZSRCS))

UTIL_HEADERS = $(shell find $(PWD)/$(INCUITLDIR) -name '*.$(HDREXT)' )

#	BUILD OUTPUT
OUTPUT  = $(OBJDIR)/*.o $(OBJPDFDIR)/*.o $(EXEDIR)/fitting $(LIBDIR)/*.so $(OBJDIR)/rapidfit_dict.* $(EXEDIR)/RapidPlot $(EXEDIR)/print

#################
##Dependencies

LINKER=ld
LINKFLAGS= -lpthread -Wl,-export-dynamic
LIBS=-static-libstdc++

CXXFLAGSUTIL = $(CXXFLAGS_BASE) -I$(INCUTILS) $(ROOTCFLAGS) -Iframework/include
CXXFLAGS     = $(CXXFLAGS_BASE) -I$(INCDIR) -I$(INCPDFDIR) -I$(INCDALITZDIR) -I$(INCGSL) $(ROOTCFLAGS)
CXXFLAGS_LIB = $(CXXFLAGS_BASE) -I$(INCDIR) -I$(INCPDFDIR) -I$(INCDALITZDIR) -I$(INCGSL) $(ROOTCFLAGS)

LIBLINKFLAGS = -pie -m64

# OS X
DARWIN=Darwin
ifeq ($(UNAME),$(Darwin))
	CXXFLAGS+= -fPIE
	CXXFLAGSUTIL+= -fPIE
	LINKFLAGS+= $(shell if [ "$(shell root-config --arch | grep 32)" = "" ]; then echo " -m64"; else echo ""; fi) -Wl,-rpath,$(LD_LIBRARY_PATH)
else
	CXXFLAGS+= -fPIE
	CXXFLAGSUTIL+= -fPIE
	LINKFLAGS+= -pie -m64 -Wl,-rpath,$(LD_LIBRARY_PATH)
endif





#	Default build command when someone asks for 'make'
all : $(EXEDIR)/fitting utils lib

$(OBJDALITZDIR)/%.o : $(SRCDALITZDIR)/%.$(SRCDALITZEXT) $(INCDALITZDIR)/%.$(HDRDALITZEXT)
	$(CXX) $(CXXFLAGS) $(USE_GSL) $(INCGSL) -c $< -o $@

$(OBJPDFDIR)/%.o : $(SRCPDFDIR)/%.$(SRCEXT) $(INCPDFDIR)/%.$(HDREXT)
	$(CXX) $(CXXFLAGS) $(USE_GSL) $(INCGSL) -c $< -o $@

#	Some fairly cool Makefile code for linking the build together :D
$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(INCDIR)/%.$(HDREXT)
	$(CXX) $(CXXFLAGS) $(USE_GSL) $(INCGSL) -c $< -o $@

$(OBJUTILDIR)/%.o : $(UTILSSRC)/%.$(UTILSSRCEXT) $(INCUTILS)/%.$(HDREXT)
	$(CXX) $(CXXFLAGSUTIL) $(USE_GSL) $(INCGSL) -c $< -o $@

#	Main Build of RapidFit Binary
$(EXEDIR)/fitting : $(OBJS) $(PDFOBJS) $(DALITZOBJS) $(OBJDIR)/rapidfit_dict.o
	$(CXX) $(LINKFLAGS) -o $@ $^ $(LIBS) $(USE_GSL) $(ROOTLIBS) $(EXTRA_ROOTLIBS) $(LINKGSL)
	chmod +t $(EXEDIR)/fitting


#	Does anyone use this any more?
doc : $(OBJS) $(PDFOBJS)
	cd doc; doxygen RapidFit_doxygen.cfg; tar cvfz RapidFit_html.tgz html/; scp RapidFit_html.tgz ph-ppe:~/WWW/RapidFit/doc/; cd ..



#	Cleanup
clean   :	distclean
cleanall:	distclean
distclean:
	$(RM) $(EXEDIR)/* $(OBJDIR)/* $(OBJPDFDIR)/* $(OBJDALITZDIR)/* $(OBJUTILDIR)/* $(LIBDIR)/*
#	$(RM) $(OUTPUT)

cleanF  :
	$(RM) $(EXEDIR)/* $(OBJDIR)/* $(LIBDIR)/*
cleanU  :
	$(RM) $(EXEDIR)/* $(OBJUTILDIR)/*
cleanP  :
	$(RM) $(OBJPDFDIR)/*


#	Allow chosing of which compiler to use on systems with mutliple compilers
clang: override CC=clang++
clang: override CXXFLAGS_BASE=-std=c++11 $(CXXFLAGS_BASE_COMMON) -O3 -msse2 -msse3 -m3dnow -ftree-vectorize -fprefetch-loop-arrays -fmerge-all-constants -Wall -Wextra -Wno-reorder
clang: all
clang-utils: override CC=clang++
clang-utils: utils

icc: override CC=icc
icc: override CXXFLAGS_BASE=$(CXXFLAGS_BASE_COMMON)
icc: all

gcc46: override CC=g++-4.6
gcc46: all
gcc47: override CC=g++-4.7
gcc47: all
gcc48: override CC=g++-4.8
gcc48: all
gcc49: override CC=g++-4.9
gcc49: all

valgrind: override CXXFLAGS_BASE+=-D__USE_VALGRIND
valgrind: all

valgrindPDF: override CXXFLAGS_BASE+=-D__USE_VALGRIND_INPDF
valgrind: all

gsl: override CXXFLAGS+= -D__RAPIDFIT_USE_GSL -D__RAPIDFIT_USE_GSL_MATH $(gsl-config --cflags)
gsl: override LINKFLAGS+= -L/sw/lib/lcg/external/GSL/1.10/x86_64-slc5-gcc43-opt/lib -lgsl -lgslcblas -lm $(gsl-config --libs)
gsl: all

#	Have a build option that SCREAMS at the user for potential mistakes!!!
debug: override CXXFLAGS+= -Wall -Wextra -Wabi -Weffc++ -ggdb -Wno-reorder -Wunused-function -Wunused-label -Wunused-value -Wunused-variable -DRAPIDFIT_USETGLTIMER
debug: override EXTRA_ROOTLIBS+= -lRGL
debug: all
debug-utils: override CXXFLAGS+= -Wall -Wextra -Wabi -Weffc++
debug-utils: utils



#	Binaries sharing LOTS of useful code for processing/outputting results
SHARED_UTIL_LIBS=$(OBJDIR)/StringProcessing.o $(OBJUTILDIR)/TTree_Processing.o $(OBJUTILDIR)/Mathematics.o  $(OBJUTILDIR)/ROOT_File_Processing.o $(OBJUTILDIR)/Histo_Processing.o $(OBJDIR)/EdStyle.o $(OBJUTILDIR)/StringOperations.o  $(OBJUTILDIR)/Template_Functions.o $(OBJUTILDIR)/RapidFit_Output_File.o $(OBJUTILDIR)/XMLUtilFunctions.o $(OBJUTILDIR)/utilsDict.o $(OBJUTILDIR)/RapidLL.o $(OBJUTILDIR)/Rapid2DLL.o $(OBJUTILDIR)/Toy_Study.o

#       New mostly automated plotting tool taking the pain out of plotting RapidFit output
$(EXEDIR)/RapidPlot: $(OBJUTILDIR)/RapidPlot.o $(OBJUTILDIR)/DoFCAnalysis.o $(OBJUTILDIR)/OutputPlots.o $(OBJUTILDIR)/CorrMatrix.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)


$(EXEDIR)/RapidToyDiff: $(OBJUTILDIR)/RapidToyDiff.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)

$(EXEDIR)/RapidDiff: $(OBJUTILDIR)/RapidDiff.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)

###   Various tools for the utils directory

#       Tool for printing information about a ROOT file and it's contents
$(EXEDIR)/print: $(OBJUTILDIR)/print.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)

###	Tool for extracting the angular distribution from the sidebands in JpsiPhi
$(EXEDIR)/AngularDist: $(OBJUTILDIR)/AngularDist.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)

###	Tool for making an exact comparison between 2 tuple branches and plotting the output
$(EXEDIR)/tupleDiff: $(OBJUTILDIR)/tupleDiff.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)

$(EXEDIR)/Compare: $(OBJUTILDIR)/Compare.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)

$(EXEDIR)/ApplyWeights: $(OBJUTILDIR)/ApplyWeights.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)

$(EXEDIR)/weighted: $(OBJUTILDIR)/weighted.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)

###	Tool for calculating gamma/deltaGamma from angular momentum analysis
$(EXEDIR)/lifetime_tool: $(OBJUTILDIR)/lifetime_tool.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAFS) $(ROOTLIBS)

$(EXEDIR)/Per-Event: $(OBJUTILDIR)/Per-Event.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)

#       Tool for printing information about a ROOT file and it's contents
$(EXEDIR)/plotDists: $(OBJUTILDIR)/plotDists.o $(SHARED_UTIL_LIBS)
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)


utils:	$(EXEDIR)/print $(EXEDIR)/RapidPlot $(EXEDIR)/RapidDiff $(EXEDIR)/RapidToyDiff

extra:	$(EXEDIR)/Per-Event $(EXEDIR)/lifetime_tool $(EXEDIR)/weighted $(EXEDIR)/ApplyWeights $(EXEDIR)/Compare $(EXEDIR)/tupleDiff $(EXEDIR)/AngularDist $(EXEDIR)/plotDists



#	For building RapidFit as a library to use within CINT which makes life easier on the grid... (supposedly)
#	make lib

lib:    $(LIBDIR)/libRapidRun.so $(LIBDIR)/libUtils.so


#	This command will generate a C++ file which interfaces the rest of humanity with root...
#	It requires the explicit paths of all files, or that you remain in the same working directory at all times during the build process
#	We want to place the output dictionary in the Build directory as this is CODE that is NOT to be editted by the $USER!
$(OBJDIR)/rapidfit_dict.cpp: framework/include/RapidRun.h framework/include/LinkDef.h
	@echo "Building RapidFit Root Dictionary:"
	@echo "rootcling -f $(OBJDIR)/rapidfit_dict.cpp -c -I\"$(PWD)/framework/include\" $^"
	@rootcling -f $(OBJDIR)/rapidfit_dict.cpp -c -I"$(PWD)/framework/include" $^

#	Compile the class that root has generated for us which is the linker interface to root	(i.e. dictionaries & such)
$(OBJDIR)/rapidfit_dict.o: $(OBJDIR)/rapidfit_dict.cpp
	$(CXX) $(CXXFLAGS) -o $@ -I"$(PWD)" -c $<

#	Class which has a dictionary generated for it, think of this as the equivalent to int main() in a CINT-y Universe
$(OBJDIR)/RapidRun.o: $(SRCDIR)/RapidRun.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#	Finally, Compile RapidFit as a library making use of the existing binaries for other classes
$(LIBDIR)/libRapidRun.so: $(OBJDIR)/RapidRun.o $(OBJDIR)/rapidfit_dict.o $(OBJS) $(PDFOBJS) $(DALITZOBJS)
	$(CXX) $(LIBLINKFLAGS) $(LINKFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(EXTRA_ROOTLIBS) $(LINKGSL)


#	This Generates the C++ file which is the interface between ROOT and the utils library

$(OBJUTILDIR)/utilsDict.cpp: $(UTILHEADERS) $(INCUTILS)/LinkDef.h
	@echo "Building Utils Root Dictionary:"
	@echo "rootcling -f $(OBJUTILDIR)/utilsDict.cpp -c -I"$(PWD)" $^"
	@rootcling -f $(OBJUTILDIR)/utilsDict.cpp -c -I"$(PWD)" $^

$(OBJUTILDIR)/utilsDict.o: $(OBJUTILDIR)/utilsDict.cpp
	$(CXX) $(CXXFLAGSUTIL) -o $@ -I"$(PWD)" -c $<

#	This is the Utils library which exposes a LOT of pre-written useful functions to the user
$(LIBDIR)/libUtils.so: $(SHARED_UTIL_LIBS) $(OBJUTILDIR)/print.o
	$(CXX) $(LIBLINKFLAGS) $(LINKFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(EXTRA_ROOTLIBS) $(LINKGSL)

