Version 1.0.1 (January 12, 2023)
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


