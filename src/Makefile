#enclose multiple programs in quotes
PROGRAMS='$(subst .cxx,, $(shell cd ./; ls *.cxx 2>/dev/null))'
SRCS = $(shell ls *.cxx 2>/dev/null)
OBJ= $(foreach PP, $(PROGRAMS), $(PP).o)

LD=g++

INCLUDEFLAGS = -I./

LD_FLAGS = -lmycode

ROOTLIBS = `root-config --libs` -lFoam -lMinuit -lTreePlayer
ROOTCFLAGS = $(shell root-config --cflags) -g3
ROOTGLIBS = $(shell root-config --glibs)
PWD = $(shell pwd)

.cxx.o:
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(INCLUDEFLAGS) -c $<

all: $(PROGS)
	$(CXX) MMPlots.cxx $(CXXFLAGS) $(ROOTCFLAGS) $(INCLUDEFLAGS) -c $<
	$(LD) -o MMPlots MMPlots.o $(ROOTLIBS)
clean:
	export PROGS=$(PROGRAMS);
	-rm MMPlots *.o *~
