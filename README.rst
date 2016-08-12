==========
pgf90_test
==========

Simple test of the correctness of the PGI/pgf90 compiler.

The testing software compilation works with Intel (v14) and GNU Fortran (v4.7) compilers:

::

  gfortran  -c src/wrapper_sorting.F90 
  ifort -c src/wrapper_sorting.F90


However, it fails with the pgf90 (here v13.10).
The TPR 22465 issue is attached to this PGI problem.

::

 pgf90 -c src/wrapper_sorting.F90 

 pgf90-Fatal-/mnt/apps/pgi/linux86-64/13.10/bin/pgf901 TERMINATED by signal 11
 Arguments to /mnt/apps/pgi/linux86-64/13.10/bin/pgf901
 /mnt/apps/pgi/linux86-64/13.10/bin/pgf901 wrapper_sorting.F90 -opt 1 -terse 1 -inform warn -nohpf -nostatic -x 19 0x400000 -quad -x 59 4 -x 59 4 -x 15 2 -x 49 0x400004 -x 51 0x20 -x 57 0x4c -x 58 0x10000 -x 124 0x1000 -tp nehalem -x 57 0xfb0000 -x 58 0x78031040 -x 48 4608 -x 49 0x100 -x 120 0x200 -stdinc /mnt/apps/pgi/linux86-64/13.10/include-gcc41:/mnt/apps/pgi/linux86-64/13.10/include:/usr/local/include:/usr/lib/gcc/x86_64-redhat-linux/4.4.7/include:/usr/lib/gcc/x86_64-redhat-linux/4.4.7/include:/usr/include -def unix -def __unix -def __unix__ -def linux -def __linux -def __linux__ -def __NO_MATH_INLINES -def __x86_64 -def __x86_64__ -def __LONG_MAX__=9223372036854775807L -def '__SIZE_TYPE__=unsigned long int' -def '__PTRDIFF_TYPE__=long int' -def __THROW= -def __extension__= -def __amd_64__amd64__ -def __k8 -def __k8__ -def __SSE__ -def __MMX__ -def __SSE2__ -def __SSE3__ -def __SSSE3__ -preprocess -freeform -vect 48 -y 54 1 -x 70 0x40000000 -modexport /tmp/pgf90G4RgGEVNiJ3E.cmod -modindex /tmp/pgf90G4RgGsRvnp4O.cmdx -output /tmp/pgf90G4RgGk_3HKdo.ilm


