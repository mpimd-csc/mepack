Version 1.0.4 (May 5, 2023)
===========================

- Fix mepack.pc package config file
- Fix OpenMP settings while configuring
- Fix: Error in computing workspace for D/SLA_TRSYLV_LEVEL3_2S
- Fix: Working with Intel ifx / OneAPI

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


