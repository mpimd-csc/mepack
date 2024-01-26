Version 1.1.1 (January 26, 2024)
=================================

- Enable Windows building with MingW
- Enable mex files for Windows with MingW
- Fix compilation error in the examples/ folder, see GH #2
- Fix Hessenberg-Triangular Detection in MATLAB/Octave interfaces


Version 1.1.0 (October 18, 2023)
================================

- Changed default block size for level-2 code from 32 to 64
- Add routines to compute the relative residual
- Fix mepack.pc package config file
- Fix OpenMP settings while configuring
- Fix wrong block size setting in examples/triangular/benchmark_trsylv.c
- Standard Equation solvers now allow Hessenberg coefficient matrices on input
- Add version script
- Fix: Error in computing workspace for D/SLA_TRSYLV_LEVEL3_2S
- Fix: Working with Intel ifx / OneAPI
- Fix: D/SLA_TCGSYLV_L2 wrong results in the (T,T) case
- Fix: mepack_memory: wrong results if M or N = 0
- Fix: mepack_memory_frontend: wrong results if M or N = 0
- MATLAB/Octave checks for Hessenberg (Triangular) Matrix automatically.

Version 1.0.3 (January 31, 2023)
================================

- Wrong CMAKE setting for installing Fortran modules

Version 1.0.2 (January 25, 2023)
================================

- MATLAB and GNU Octave has wrong library dependencies inside their .so file.


Version 1.0.1 (January 20, 2023)
================================

- Disable building the static library by default (see BUILD_STATIC in
  install.md)
- Enable -Wmaybe-uninitialized compiler flag if supported
- Enable -Wstringop-trunction compiler flag if supported
- Enable -Wconversion in Fortran compiler flag if supported
- Fix LTO warnings and lazy C-Fortran calls.
- LTO enabled by default
- Fortran modules are installed in the same location as the C headers.

Version 1.0.0 (December 6, 2022)
================================

- Initial Release


