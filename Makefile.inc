################################################################################
#
#	APFELcomb Makefile includes
#
################################################################################

CC= g++
LD= g++

SHELL=/bin/bash

# root
ROOTINCS = $(shell root-config --cflags) 
ROOTLIBS = $(shell root-config --glibs) 

# APFEL
APFELINCS = $(shell apfel-config --cppflags) -I ./src
APFELLIBS = $(shell apfel-config --ldflags) 

#LHAPDF
LHAPDFINCS = -I$(shell lhapdf-config --prefix)/include
LHAPDFDIR  = $(shell lhapdf-config --prefix)/lib
LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF

# applgrid
APPLINCS = -I$(shell applgrid-config --prefix)/include
APPLCLIBS = -L$(shell applgrid-config --prefix)/lib -lAPPLgrid 

# additional libraries to be included 
PRJLDFLAGS = $(LHAPDFLIBS) $(APPLCLIBS) $(ROOTLIBS) $(APFELLIBS)

# scheduling and optimization options
PRJCXXFLAGS = -Wall -O3 $(ALLDIR) $(APPLINCS) $(ROOTINCS) $(APFELINCS) $(LHAPDFINCS)
