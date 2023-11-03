# Path to Pythia8 installation
include Makefile.inc
GFORTRAN=/usr/local/Cellar/gcc/10.2.0_4/lib/gcc/10
PYTHIAXMLDIR=$(PYTHIADIR)/share/Pythia8/xmldoc

CXX=g++
CXXFLAGS=-g -O2 -std=c++17
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
hermes_collins: hermes_collins.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

compass_zcollins: compass_zcollins.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

compass_zsivers: compass_zsivers.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl
compass_ptcollins: compass_ptcollins.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

hermes_zcollins: hermes_zcollins.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

hermes_ptcollins: hermes_ptcollins.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

compass_ptsivers: compass_ptsivers.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

hermes_sivers: hermes_sivers.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

hermes_zsivers: hermes_zsivers.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl

hermes_ptsivers: hermes_ptsivers.cc StringSpinner.h Transversity.h VectorMesonDecays.h PrimordialKT.h mc3P0.o def.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDEDIR) -o $@ $< mc3P0.o def.o -L$(GFORTRAN) -lgfortran -L$(LIBDIR) -Wl,-rpath $(LIBDIR) -lpythia8 -ldl
def.o: definitions.f90 
	$(FC) -c $< -o $@ 

mc3P0.o: mc3P0.f90 def.o
	$(FC) -c $(FFLAGS) $< def.o

clean:
	rm -f mc3P0.o dis routines.mod m20.o def.o event.o routines_and_functions.mod eventdefinition.mod newvariables.mod parameters.o parameters.mod 

