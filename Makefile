CXX := g++
CXXFLAGS := -O3 -g -std=c++11 -frounding-math

GSL_DIR := /Users/arthur/Software/gsl
NLOPT_DIR := /Users/arthur/Software/nlopt
CUBA_DIR := /Users/arthur/Software/cuba

INC_DIRS := -I$(GSL_DIR)/include 	-I$(NLOPT_DIR)/include 	-I$(CUBA_DIR)/include
LIB_DIRS := -L$(GSL_DIR)/lib 		-L$(NLOPT_DIR)/lib		-L$(CUBA_DIR)/lib

LFLAGS := -lm -lgsl -lgslcblas -lcuba -lnlopt

all: fibez

fibez: *.h *.cxx
	$(CXX) $(CXXFLAGS) $(INC_DIRS) $(LIB_DIRS) main.cxx -o fibez $(LFLAGS)

clean:
	rm -rf fibez* *.lammpstrj log.fibez

#You may need to include the NLopt library in the dynamic library path on some architectures, i.e.
#export DYLD_LIBRARY_PATH=$(NLOPT_DIR)/lib:$DYLD_LIBRARY_PATH
