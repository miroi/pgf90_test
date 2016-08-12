==========
pgf90_test
==========

Simple test of the correctness of the PGI/pgf90 compiler.

The TPR 22465 issue is attached to the PGI problem.

The simple compilation works Intel and GNU Fortran compilers:

::

  gfortran  -c wrapper_sorting.F90 
  ifort -c wrapper_sorting.F90


However, it fails with pgf90:




