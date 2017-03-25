# These should come from environment:
# CXX
# OPTFLAGS

# TMV_DIR
# YAML_DIR
# MKL_DIR (optionally)
# It is assumed that the gbtools/utilities/ directory is at ../utilities 

export CXX
export OPTFLAGS

# ABS_INCLUDES are absolute paths 
ABS_INCLUDES = -I $(TMV_DIR)/include -I $(YAML_DIR)/include

ifdef MKL_DIR
ABS_INCLUDES += -I $(MKL_DIR)/include
endif

SUBDIRS = 

INCLUDES = -I ../utilities

CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)

SRC = $(shell ls *.cpp)

OBJ = Astrometry.o Ephemeris.o PixelMap.o PolyMap.o PixelMapCollection.o SubMap.o Wcs.o \
      TemplateMap.o PiecewiseMap.o AlphaUpdater.o YAMLCollector.o

all:  $(OBJ)

# For building test programs:
UTILITIES := ../utilities
SUBOBJ = $(UTILITIES)/StringStuff.o $(UTILITIES)/Poly2d.o 
LIB_DIRS = -L $(TMV_DIR)/lib -L -L $(YAML_DIR)/lib
TMV_LINK := $(shell cat $(TMV_DIR)/share/tmv/tmv-link)
CXXFLAGS = $(OPTFLAGS) $(ABS_INCLUDES) $(INCLUDES)
LIBS = -lm $(LIB_DIRS) -lyaml-cpp -ltmv_symband $(TMV_LINK)

testSphericalYAML: testSphericalYAML.o Astrometry.o $(SUBOBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@
testCollector: testCollector.o YAMLCollector.o $(SUBOBJ)
	$(CXX) $(CXXFLAGS) $^  $(LIBS) -o $@


###############################################################
## Standard stuff:
###############################################################


subs:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE)); done

depend:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) depend); done
	$(CXX) $(CXXFLAGS) -MM $(SRC) > .$@

clean:
	for dir in $(SUBDIRS); do (cd $$dir && $(MAKE) clean); done
	rm -f *.o *~ *.aux *.log *.dvi core .depend

ifeq (.depend, $(wildcard .depend))
include .depend
endif

export

.PHONY: all install dist depend clean 
