The StringSpinner package is composed of:
a) StringSpinner.h:     C++ header file that contains the implementation of the UserHooks class for the introduction of spin effects in the hadronization part of PYTHIA 8.3.
b) Transversity.h:      C++ header file that contains the definitions of the transversity parton distribution functions.
c) mc3P0.f90:           Fortran module that contains the routines necessary for the propagation of the spin effects.
d) dis.cc:              Sample main program for the simulation of the polarized DIS processs.
e) VectorMesonDecays.h: C++ header file with the implementation of the external decay handler class.
f) definitions.f90:     Fortran file containing definitions of 4-vectors and operations on them.
g) PrimordialKT.h:      C++ header file needed for the determination of the string axis when the primordial kT is switched on.
h) Makefile:            Makefile for the compilation of the above files.

Assuming that PYTHIA 8.3 has already been installed, to run a main program with StringSpinner

i.      git clone https://gitlab.com/albikerbizi/stringspinner
ii.     cd stringspinner
iii.    ./configure path/to/pythia/installation/directory 
iv.     make dis
iv.     ./dis

The compilation procedure produces the main program executable "dis", the Fortran module "routines.mod",
and the object files "mc3P0.o" and "def.o". These files can be deleted with "make clean".

Test example:
Running the main program dis with 1M events, one obtains for the average Collins asymmetry for pi+ the value
-0.044 +/- 0.006.

Contact:
Albi Kerbizi: albi.kerbizi@ts.infn.it
Leif Lonnblad: leif.lonnblad@thep.lu.se


