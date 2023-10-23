# Path to Pythia8 installation
include Makefile.inc
GFORTRAN=/usr/local/Cellar/gcc/10.2.0_4/lib/gcc/10
PYTHIAXMLDIR=$(PYTHIADIR)/share/Pythia8/xmldoc

CXX=g++
CXXFLAGS=-g -O -std=c++11
INCLUDEDIR=$(PYTHIADIR)/include
LIBDIR=$(PYTHIADIR)/lib
FC=gfortran
FFLAGS=-O0 -g -frecord-marker=8 -fbounds-check

dis: dis.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

collins: collins.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

sivers: sivers.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl
hermes: hermes.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

z: z.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

zsivers: zsivers.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl
pt: pt.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl


ptsivers: ptsivers.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl
def.o: definitions.f90 
	$(FC) -c $< -o $@ 

mc3P0.o: mc3P0.f90 def.o
	$(FC) -c $(FFLAGS) $< def.o

clean:
	rm -f mc3P0.o dis routines.mod m20.o def.o event.o routines_and_functions.mod eventdefinition.mod newvariables.mod parameters.o parameters.mod 

