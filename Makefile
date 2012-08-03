# $Id: Makefile,v 1.30 2009/11/11 17:18:14 gcowan Exp $
SHELL = /bin/bash
UNAME = $(shell uname)
CC = g++

#		ROOT
TEMPCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)

#		Include Root files as system headers as they're NOT standards complient and we do not want to waste time fixing them!
#		ROOT has some broken backwards compatability for OSX so won't claim to be a set of system headers
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags | awk -F "-I" '{print $$1" -isystem"$$2}' )
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

#               On some Systems with Mathmore compiled, sometimes things need to be resolved against it... I don't know why
EXTRA_ROOTLIBS=-lTreePlayer -lHtml -lThread -lMinuit -lMinuit2 -lRooFit -lRooStats -lRooFitCore -lFoam $(shell if [ "$(shell root-config --features | grep mathmore)" = "" ]; then echo "" ; else echo "-lMathMore" ; fi)

#		Command Line Tools
CXX          = $(CC) $(shell if [ "$(shell root-config --arch | grep 32)" = "" ]; then echo ""; else echo "--arch=i386"; fi)
RM           = rm -f

SVN_REV ="$(shell svnversion -n .)"

#		Compiler Flags
CXXFLAGS_BASE  = -DSVN_REV=$(SVN_REV) -rdynamic -D_GNU_SOURCE -D__USE_GNU -fPIC -O3 -msse -msse2 -msse3 -m3dnow -g -ansi -fmerge-all-constants -funroll-all-loops -fno-common -D__ROOFIT_NOBANNER -Wconversion -Wextra -Wsign-compare -Wfloat-equal -Wmissing-noreturn -Wall -Wno-non-virtual-dtor -Wno-reorder -pthread -Wshadow -Wcast-align

CXX_FLAGS_LITE = -DSVN_REV=$(SVN_REV) -rdynamic -D_GNU_SOURCE -D__USE_GNU -fPIC -Os -msse -msse2 -msse3 -m3dnow -g -ansi -fmerge-all-constants -D__ROOFIT_NOBANNER -Wconversion -Wextra -Wsign-compare -Wfloat-equal -Wmissing-noreturn -Wall -Wno-non-virtual-dtor -Wno-reorder -pthread -Wshadow -Wcast-align

#		Some Useful global variables, makes this file MUCH easier to maintain
SRCEXT    = cpp
HDREXT    = h
SRCDIR    = framework/src
SRCPDFDIR = pdfs/src
UTILSSRC  = utils/src
INCDIR    = framework/include
INCPDFDIR = pdfs/include
INCUTILS  = utils/include
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
#	PDF source files which will be required for building libpdf.so
PDFSRCS := $(shell find $(SRCPDFDIR) -name '*.$(SRCEXT)' | grep -v 'unused' )
DALITZSRCS := $(shell find $(SRCDALITZDIR) -name '*.$(SRCDALITZEXT)' | grep -v 'unused' )


#	Absolute Paths of headers	ignoring the LinkDef written for ROOT and ignoring unused code
HEADERS := $(shell find $(INCDIR) -name '*.$(HDREXT)' | grep -v 'unused' | grep -v 'LinkDef' )
PDFHEAD := $(shell find $(INCPDFDIR) -name '*.$(HDREXT)' )
DALITZHEAD := $(shell find $(INCDALITZDIR) -name '*.$(HDRDALITZEXT)' )


#	Binary Objects to make in the build
OBJS    := $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))
#	Binary objects to be linked into libpdf.so
PDFOBJS := $(patsubst $(SRCPDFDIR)/%.$(SRCEXT),$(OBJPDFDIR)/%.o,$(PDFSRCS))
DALITZOBJS := $(patsubst $(SRCDALITZDIR)/%.$(SRCDALITZEXT),$(OBJDALITZDIR)/%.o,$(DALITZSRCS))



#	All Headers in their absolute paths as required by CINT to construct a dictionary of RapidFit
ALL_HEADERS += $(HEADERS)
ALL_HEADERS += $(PDFHEAD)
ALL_HEADERS += $(DALITZHEAD)

UTIL_HEADERS = $(shell find $(PWD)/$(INCUITLDIR) -name '*.$(HDREXT)' )

#	BUILD OUTPUT
OUTPUT  = $(OBJDIR)/*.o $(OBJPDFDIR)/*.o $(EXEDIR)/fitting $(LIBDIR)/*.so $(OBJDIR)/rapidfit_dict.* $(EXEDIR)/rapidfit_toyresults $(EXEDIR)/rapidfit_fcscanresults $(EXEDIR)/rapidfit_fcscanresults_2 $(EXEDIR)/betas_sweightfitter $(EXEDIR)/RapidLL $(EXEDIR)/RapidPlot $(EXEDIR)/print



#################
##Dependencies

LINKER=ld
LINKFLAGS=-lpthread

LIBS=-lstdc++

CXXFLAGSUTIL = $(CXXFLAGS_BASE) -I$(INCUTILS) $(ROOTCFLAGS) -Iframework/include
CXXFLAGS     = $(CXXFLAGS_BASE) -I$(INCDIR) -I$(INCPDFDIR) -I$(INCDALITZDIR) $(ROOTCFLAGS)
CXXFLAGS_LIB = $(CXXFLAGS_BASE) -I$(INCDIR) -I$(INCPDFDIR) -I$(INCDALITZDIR) $(ROOTCFLAGS)

# Linux
ifeq "$(UNAME)" "Linux"
	GCC_V:=$(shell gcc -dumpversion | awk -F '.' '{print $$2}')
	CXX_LTO:=$(shell if [ ${GCC_V} -ge 6 ]; then echo '-flto '; else echo ''; fi)
	CXXFLAGS+=${CXX_LTO}-fPIE
	LINKFLAG+= -flto -pie -m64
endif

# OS X
ifeq "$(UNAME)" "Darwin"
	CXXFLAGS+= -fPIE
	LINKFLAGS+= $(shell if [ "$(shell root-config --arch | grep 32)" = "" ]; then echo " -m64"; else echo ""; fi)
endif





#	Default build command when someone asks for 'make'
all : $(EXEDIR)/fitting utils lib 

$(OBJDALITZDIR)/%.o : $(SRCDALITZDIR)/%.$(SRCDALITZEXT) $(INCDALITZDIR)/%.$(HDRDALITZEXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPDFDIR)/%.o : $(SRCPDFDIR)/%.$(SRCEXT) $(INCPDFDIR)/%.$(HDREXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

#	Some fairly cool Makefile code for linking the build together :D
$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(INCDIR)/%.$(HDREXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

#	Main Build of RapidFit Binary
$(EXEDIR)/fitting : $(OBJS) $(PDFOBJS) $(DALITZOBJS) $(OBJDIR)/rapidfit_dict.o
	$(CXX) $(LINKFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(EXTRA_ROOTLIBS)






#	Does anyone use this any more?
doc : $(OBJS) $(PDFOBJS)
	cd doc; doxygen RapidFit_doxygen.cfg; tar cvfz RapidFit_html.tgz html/; scp RapidFit_html.tgz ph-ppe:~/WWW/RapidFit/doc/; cd ..



#	Cleanup
clean   :	distclean
cleanall:	distclean
distclean:
	$(RM) $(EXEDIR)/* $(OBJDIR)/* $(OBJPDFDIR)/* $(OBJDALITZDIR)/* $(OBJUTILDIR)/* $(LIBDIR)/*
#	$(RM) $(OUTPUT)


#	Allow chosing of which compiler to use on systems with mutliple compilers
clang: override CC=clang++
clang: all
clang-utils: override CC=clang++
clang-utils: utils

gcc46: override CC=g++-4.6
gcc46: all
gcc47: override CC=g++-4.7
gcc47: all
gcc48: override CC=g++-4.8
gcc48: all

gsl: override CXXFLAGS+= -D__USE_GSL_ERR $(gsl-config --cflags) 
gsl: override LINKFLAGS+= -lgsl -lgslcblas -lm $(gsl-config --libs)
gsl: all

#	Have a build option that SCREAMS at the user for potential mistakes!!!
debug: override CXXFLAGS+= -Wall -Wextra -Wabi -Weffc++ -ggdb -Wno-reorder
debug: all
debug-utils: override CXXFLAGS+= -Wall -Wextra -Wabi -Weffc++
debug-utils: utils



#	Binaries sharing LOTS of useful code for processing/outputting results
$(OBJUTILDIR)/TTree_Processing.o: $(UTILSSRC)/TTree_Processing.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/ROOT_File_Processing.o: $(UTILSSRC)/ROOT_File_Processing.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/StringOperations.o: $(UTILSSRC)/StringOperations.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/Histo_Processing.o: $(UTILSSRC)/Histo_Processing.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/RapidFit_Output_File.o: $(UTILSSRC)/RapidFit_Output_File.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/Mathematics.o: $(UTILSSRC)/Mathematics.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<

#       New mostly automated plotting tool taking the pain out of plotting RapidFit output
$(EXEDIR)/RapidPlot: $(OBJUTILDIR)/RapidPlot.o $(OBJUTILDIR)/DoFCAnalysis.o $(OBJUTILDIR)/Mathematics.o $(OBJUTILDIR)/OutputPlots.o $(OBJUTILDIR)/RapidLL.o $(OBJUTILDIR)/Rapid2DLL.o $(OBJUTILDIR)/Toy_Study.o $(OBJDIR)/EdStyle.o $(OBJDIR)/StringProcessing.o $(OBJUTILDIR)/TTree_Processing.o $(OBJUTILDIR)/ROOT_File_Processing.o $(OBJUTILDIR)/Histo_Processing.o $(OBJUTILDIR)/StringOperations.o $(OBJUTILDIR)/Component_Projections.o $(OBJUTILDIR)/RapidFit_Output_File.o $(OBJUTILDIR)/Template_Functions.o
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)
$(OBJUTILDIR)/RapidPlot.o: $(UTILSSRC)/RapidPlot.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/RapidLL.o: $(UTILSSRC)/RapidLL.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/Rapid2DLL.o: $(UTILSSRC)/Rapid2DLL.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/Toy_Study.o: $(UTILSSRC)/Toy_Study.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/Component_Projections.o: $(UTILSSRC)/Component_Projections.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/DoFCAnalysis.o: $(UTILSSRC)/DoFCAnalysis.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/OutputPlots.o: $(UTILSSRC)/OutputPlots.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<
$(OBJUTILDIR)/Template_Functions.o: $(UTILSSRC)/Template_Functions.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<

#       Tool for printing information about a ROOT file and it's contents
$(EXEDIR)/print: $(OBJUTILDIR)/print.o $(OBJDIR)/EdStyle.o $(OBJDIR)/StringProcessing.o $(OBJUTILDIR)/TTree_Processing.o $(OBJUTILDIR)/Mathematics.o  $(OBJUTILDIR)/ROOT_File_Processing.o $(OBJUTILDIR)/Histo_Processing.o $(OBJUTILDIR)/StringOperations.o $(OBJUTILDIR)/Template_Functions.o
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)
$(OBJUTILDIR)/print.o: $(UTILSSRC)/print.C
	$(CXX) $(CXXFLAGS) -I$(INCUTILS) -o $@ -c $<

$(EXEDIR)/weighted: $(OBJUTILDIR)/weighted.o $(OBJUTILDIR)/TTree_Processing.o $(OBJUTILDIR)/Mathematics.o $(OBJUTILDIR)/ROOT_File_Processing.o $(OBJUTILDIR)/Histo_Processing.o $(OBJUTILDIR)/StringOperations.o $(OBJUTILDIR)/Template_Functions.o
	$(CXX) -o $@ $^ $(LINKFLAGS) $(ROOTLIBS)
$(OBJUTILDIR)/weighted.o: $(UTILSSRC)/weighted.C
	$(CXX) $(CXXFLAGSUTIL) -o $@ -c $<

utils: $(EXEDIR)/weighted $(EXEDIR)/print $(EXEDIR)/RapidPlot






#	For building RapidFit as a library to use within CINT which makes life easier on the grid... (supposedly)
#	make lib

lib:    $(LIBDIR)/libRapidRun.so


#	This command will generate a C++ file which interfaces the rest of humanity with root...
#	It requires the explicit paths of all files, or that you remain in the same working directory at all times during the build process
#	We want to place the output dictionary in the Build directory as this is CODE that is NOT to be editted by the $USER!
$(OBJDIR)/rapidfit_dict.cpp: $(ALL_HEADERS) framework/include/LinkDef.h
	@echo "Building Root Dictionary:"
	@echo "rootcint -f $(OBJDIR)/rapidfit_dict.cpp -c -I\"$(PWD)/framework/include\" $^"
	@rootcint -f $(OBJDIR)/rapidfit_dict.cpp -c -I"$(PWD)/framework/include" $^

#	Compile the class that root has generated for us which is the linker interface to root	(i.e. dictionaries & such)
$(OBJDIR)/rapidfit_dict.o: $(OBJDIR)/rapidfit_dict.cpp
	$(CXX) $(CXXFLAGS) -o $@ -I"$(PWD)" -c $<

#	Class which has a dictionary generated for it, think of this as the equivalent to int main() in a CINT-y Universe
$(OBJDIR)/RapidRun.o: $(SRCDIR)/RapidRun.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#	Finally, Compile RapidFit as a library making use of the existing binaries for other classes
$(LIBDIR)/libRapidRun.so: $(OBJDIR)/RapidRun.o $(OBJDIR)/rapidfit_dict.o $(OBJS) $(PDFOBJS) $(DALITZOBJS) $(OBJDIR)/ClassLookUp.o
	$(CXX) $(LINKFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(EXTRA_ROOTLIBS)

