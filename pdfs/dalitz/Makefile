#--------------------------------------------------------------------------------------------------------------
# get the correct root configuration
#
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

SYSTEM=$(shell uname -s)

HEADERS=$(wildcard include/*.hh) 
SOURCES=$(wildcard src/*.cc)
OBJECTS=$(patsubst src/%.cc, tmp/%.o, $(SOURCES))

#--------------------------------------------------------------------------------------------------------------------------------------------
# lib to be built
#
LIB=lib/libDPAmpCalc.so
TBIN=test/test.exe test/refAxisTest.exe

#--------------------------------------------------------------------------------------------------------------------------------------------
#
# binaries to be built
#
default: lib bin

lib: $(LIB)

bin: lib $(TBIN)

#-----------------------------------------------------------------------------------------
#
# compile flags
#
CXX	    = g++
#CXXFLAGS    = -g -O -fPIC -Wall $(ROOTCFLAGS) -I.  -Iinclude
CXXFLAGS    = -pg -g  -fPIC -Wall $(ROOTCFLAGS) -I.  -Iinclude 

#--------------------------------------------------------------------------------------------------------------------------------------------
#
# set the linker options
#
LD	    = g++
LDFLAGS     = -O -g

LIBS        = -L. $(ROOTLIBS) $(SYSLIBS) -lstdc++ -Llib \
              -lDPAmpCalc
SOFLAGS     = -dynamiclib 
ifeq ($(SYSTEM),Darwin)
	SOFLAGS+=  -undefined dynamic_lookup
endif
#SOFLAGS     = -bundle -undefined dynamic_lookup


#------------------------------------------------------------------------------
# executables
PRINTLINK=@echo "<*** linking $< ***>"

test/%.exe: tmp/%.o
	$(PRINTLINK); $(LD) $(LDFLAGS) $(LIBS) $< -o $@

#-----------------------------------------------------------------------------------------------------------------------------------------
# implicite rule for default object files
#
PRINT=@echo "<*** compiling $< ***>"

tmp/%.o: src/%.cc
	$(PRINT); $(CXX) $(CXXFLAGS) -c $< -o $@

tmp/%.o: test/%.cc
	$(PRINT); $(CXX) $(CXXFLAGS) -c $< -o $@

#-----------------------------------------------------------------------------------------------------------------------------------------
# how to build lib 
#
$(LIB): $(OBJECTS) tmp/dict.o
	$(LD) $(CXXFLAGS) $(SOFLAGS) -o $@ $^ 
#	@echo "<*** build lib ***>"; $(LD) $(CXXFLAGS) $(SOFLAGS) -o $@ $^ 


#-----------------------------------------------------------------------------------------------------------------------------------------
# rootcint stuff
#
CINTINCLUDE=-Iinclude -I$(ROOTSYS)/include

LINKDEFFILE=include/linkdef.h

tmp/dict.o: tmp/dict.cxx $(LINKDEFFILE)
	$(PRINT); $(CXX) $(CXXFLAGS) -c $< -o $@

tmp/dict.cxx: $(HEADERS) $(LINKDEFFILE)
	@echo "<*** rootcint ***>"; rm -f tmp/dict.cxx tmp/dict.h; rootcint tmp/dict.cxx -c $(CINTINCLUDE) $^

#----------------------------------------------------------------------------------------------------------------------------------------

clean:
	rm -f tmp/*.o $(FITLIB) tmp/dict.* test/*.exe

binclean:
	rm -f  bin/*.exe

veryclean: clean binclean

#------------------------------------------------------------------------------------------------------------------------------------------

#doc: macros/doc.cc src/*cc include/*hh
#	root -b -q macros/doc.cc
doc: $(SOURCES) $(HEADERS) doc/doxygen.config
	doxygen doc/doxygen.config
#------------------------------------------------------------------------------------------------------------------------------------------
