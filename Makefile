# $Id: Makefile,v 1.30 2009/11/11 17:18:14 gcowan Exp $
SHELL = /bin/bash
UNAME = $(shell uname)

#		ROOT
TEMPCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
#		Include Root files as system headers as they're NOT standards complient and we do not want to waste time fixing them!
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags | awk -F "-I" '{print $$1" -isystem"$$2}' )
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)


#		Command Line Tools
CXX          = g++
RM           = rm -f


#		Compiler Flags
CXXFLAGS     = -Os -shared -msse -msse2 -m3dnow -g1 -ansi -fPIC -shared -D__ROOFIT_NOBANNER -Wconversion -Wextra -Wsign-compare -Wfloat-equal -Wmissing-noreturn -Wall -Wno-non-virtual-dtor


#		Some Useful global variables, makes this file MUCH easier to maintain
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

#	Source Files to be Built
SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)' -not -name 'Roo*.cpp' | grep -v 'RapidRun' | grep -v 'ClassLookUp' )
PDFSRCS := $(shell find $(SRCPDFDIR) -name '*.$(SRCEXT)' -not -name 'Roo*.cpp' )

#	Absolute Paths of headers
HEADERS := $(shell find $(PWD)/$(INCDIR) -name '*.$(HDREXT)' | grep -v 'LinkDef' )
PDFHEAD := $(shell find $(PWD)/$(INCPDFDIR) -name '*.$(HDREXT)' )

#	Binary Objects to make in the build
OBJS    := $(patsubst $(SRCDIR)/%.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))
PDFOBJS := $(patsubst $(SRCPDFDIR)/%.$(SRCEXT),$(OBJPDFDIR)/%.o,$(PDFSRCS))

#	All Headers in their absolute paths as required by CINT to construct a dictionary of RapidFit
ALL_HEADERS += $(HEADERS)
ALL_HEADERS += $(PDFHEAD)



#	BUILD OUTPUT
OUTPUT  = $(OBJDIR)/*.o $(OBJPDFDIR)/*.o $(EXEDIR)/fitting $(LIBDIR)/*.so $(OBJDIR)/rapidfit_dict.* *.so *.rootmap



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



#LIBS       += $(ROOTLIBS) -lHtml -lThread -lMinuit -lRooFit -lRooFitCore -lMathCore -lMinuit2 -lFoam #-lboost_thread-xgcc40-mt #-lMathMore
LIBS       += $(ROOTLIBS) -lHtml -lThread -lMinuit -lMathCore -lMinuit2 -lRooFit -lRooFitCore -lFoam #-lProof #-lboost_thread #-lMathMore



#	Default build command when someone asks for 'make'
all : $(EXEDIR)/fitting utils




#	SPECIAL CASE This class SHOULD depend ON ALL changes made in the PDF directory as matter of course
$(OBJDIR)/ClassLookUp.o : $(INCDIR)/ClassLookUp.h $(SRCDIR)/ClassLookUp.cpp $(PDFHEAD)
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/ClassLookUp.cpp -o $(OBJDIR)/ClassLookUp.o



#	Some fairly cool Makefile code for linking the build together :D
$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(INCDIR)/%.$(HDREXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(OBJPDFDIR)/%.o : $(SRCPDFDIR)/%.$(SRCEXT) $(INCPDFDIR)/%.$(HDREXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@


#	Main Build of RapidFit Binary
$(EXEDIR)/fitting : $(OBJS) $(PDFOBJS) $(OBJDIR)/ClassLookUp.o $(OBJDIR)/RootFileDataSet.o
	$(CXX) -o $@ $^ $(LINKFLAGS) $(LIBS)







#	Does anyone use this any more?
doc : $(OBJS) $(PDFOBJS)
	cd doc; doxygen RapidFit_doxygen.cfg; tar cvfz RapidFit_html.tgz html/; scp RapidFit_html.tgz ph-ppe:~/WWW/RapidFit/doc/; cd ..



#	Cleanup
clean   :	distclean
cleanall:	distclean
distclean:
	$(RM) $(OUTPUT)




#	RapidFit Utils:

#	Tool for plotting Toy Study Results
$(EXEDIR)/rapidfit_toyresults: $(OBJDIR)/rapidfit_toyresults.o
	$(CXX) -o $@ $< $(CXXFLAGS) $(ROOTLIBS)
$(OBJDIR)/rapidfit_toyresults.o: $(UTILSSRC)/rapidfit_toyresults.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<


#	Original Tool for plotting FCScan Results
$(EXEDIR)/rapidfit_fcscanresults: $(OBJDIR)/rapidfit_fcscanresults.o
	$(CXX) -o $@ $< $(ROOTLIBS)
$(OBJDIR)/rapidfit_fcscanresults.o: $(UTILSSRC)/rapidfit_fcscanresults.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<


#	Tool for plotting FCScan Results and more
$(EXEDIR)/rapidfit_fcscanresults_2: $(OBJDIR)/rapidfit_fcscanresults_2.o
	$(CXX) -o $@ $< $(ROOTLIBS)
$(OBJDIR)/rapidfit_fcscanresults_2.o: $(UTILSSRC)/rapidfit_fcscanresults_2.cc
	$(CXX) $(CXXFLAGS) -o $@ -c $<


#	Trying to move people off this tool
#$(EXEDIR)/rapidfit_llscanresults: $(OBJDIR)/rapidfit_llscanresults.o
#	$(CXX) -o $@ $< $(ROOTLIBS)
#$(OBJDIR)/rapidfit_llscanresults.o: $(UTILSSRC)/rapidfit_llscanresults.cc
#	$(CXX) $(CXXFLAGS) -o $@ -c $<


#utils: $(EXEDIR)/rapidfit_toyresults $(EXEDIR)/rapidfit_llscanresults  $(EXEDIR)/rapidfit_fcscanresults
utils: $(EXEDIR)/rapidfit_toyresults $(EXEDIR)/rapidfit_fcscanresults $(EXEDIR)/rapidfit_fcscanresults_2






#	For building RapidFit as a library to use within CINT which makes life easier on the grid... (supposedly)
#	make lib

lib:    $(LIBDIR)/libRapidRun.so


#	This command will generate a C++ file which interfaces the rest of humanity with root...
#	It requires the explicit paths of all files, or that you remain in the same working directory at all times during the build process
#	We want to place the output dictionary in the Build directory as this is CODE that is NOT to be editted by the $USER!
$(OBJDIR)/rapidfit_dict.cpp: $(ALL_HEADERS) $(PWD)/framework/include/LinkDef.h
	@echo "Building Root Dictionary:"
	@echo "rootcint -f $(OBJDIR)/rapidfit_dict.cpp -c $^"
	@rootcint -f $(OBJDIR)/rapidfit_dict.cpp -c $^

#	Compile the class that root has generated for us which is the linker interface to root	(i.e. dictionaries & such)
$(OBJDIR)/rapidfit_dict.o: $(OBJDIR)/rapidfit_dict.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#	Class which has a dictionary generated for it, think of this as the equivalent to int main() in a CINT-y Universe
$(OBJDIR)/RapidRun.o: $(SRCDIR)/RapidRun.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

#	Finally, Compile RapidFit as a library making use of the existing binaries for other classes
$(LIBDIR)/libRapidRun.so: $(OBJDIR)/RapidRun.o $(OBJDIR)/rapidfit_dict.o $(OBJS) $(PDFOBJS) $(OBJDIR)/ClassLookUp.o
	$(CXX) $(LIBS) -shared -fPIC $(CXXFLAGS) $^ -o $@

