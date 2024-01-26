# Installation                                                      {#install}

## Requirements

### Linux and Unix-like Operating Systems

The *following* tools and libraries are required for installing MEPACK on Linux
or Unix-like operating systems:

 - A Fortran compiler supporting Fortran 95/2003, e.g. GNU Fortran, Intel ifort,
   IBM XLF, ...
 - A C99 compiler, e.g. GNU GCC, Intel icc, IBM XLC, ...
 - (*optional*) A C++14 compiler for MATLAB/GNU Octave interfaces
 - CMake 3.15 or newer. If Intel's oneAPI icx and ifx are used, CMake>=3.20 is
   required.
 - BLAS and LAPACK, at least version 3.4.2
 - (*optional*) HDF5, at least version 1.8.5 (only for tests and examples)
 - (*optional*) doxygen, at least version 1.8.16 with dot, if the HTML
   documentation is required
 - (*optional*) MATLAB for building the MATLAB interface
 - (*optional*) GNU Octave, at least version 4.4, for building the GNU Octave
   interface
 - (*optional*) lcov and gcov for building the tests with code coverage

The following compilers are known to be unusable due to compiler bugs:

 - AMD AOCC 3.x and 4.0 (Compiler Crash)
 - Intel ICC/IFort Classic 2021 (OpenMP failures)
 - Intel ICX/IFX 2021.x, 2022.x (OpenMP failures)

If the BLAS / LAPACK library does not provide the full interface, i.e. some
routines are not implemented or available for the user, MEPACK builds the
required but missing routines from reference LAPACK and includes them in MEPACK.
This is the if IBM ESSL is used as BLAS library.

### Microsoft Windows

MEPACK can be built on MS Windows using the MSYS2 (https://www.msys2.org/)
environemnt in combination with the MingW64(gcc) compiler. We support both, the
`MSVCRT` and the `UCRT` runtime environment. If you are plan to use `MSVCRT` set
```shell
export MINGWPREFIX=mingw-w64-x86_64
```
In case of the `UCRT` runtime set
```shell
export MINGWPREFIX=mingw-w64-ucrt-x86_64
```

Afterwards, the required packages can be installed using
```shell
pacman -S ${MINGWPREFIX}-gcc ${MINGWPREFIX}-gcc-fortran \
          ${MINGWPREFIX}-cmake ${MINGWPREFIX}-ninja \
          ${MINGWPREFIX}-openblas ${MINGWPREFIX}-openblas64 \
          ${MINGWPREFIX}-hdf5
```

## Installation

The source code is configured using CMake. An out-of-source build is preferred.
Typically, this is done by
```shell
cmake -S . -B build-dir  <OtherCMakeOptions>
make -C build-dir
make -C build-dir install
```

The following options can be added to `CMake` to configure MEPACK properly, for
MATLAB/GNU Octave specific options, see below:

|  Option                          | Description                                  | Default              |
|----------------------------------|----------------------------------------------|----------------------|
| `-DDEBUG=ON/OFF`                 | Enable building the library in debug mode    | `OFF`                |
| `-DBUILD_STATIC=ON/OFF`          | Enable building the static libraries.        | `OFF`                |
| `-DBUILD_DYNAMIC=ON/OFF`         | Enable building the dynamic/shared libraries.| `ON`                 |
| `-DCOVERAGE=ON/OFF`              | Enable the creation of code coverage files   | `OFF`                |
| `-DCMAKE_INSTALL_PREFIX=PREFIX`  | Installation prefix                          | `/usr/local`         |
| `-DBLA_VENDOR=VENDOR_NAME`       | Set the BLAS library to use.                 | `Generic`            |
|                                  | See FindBLAS.cmake for details               |                      |
| `-DHOSTOPT=ON/OFF`               | Enable host-specific optimizations           | `OFF`                |
| `-DLTO=ON/OFF`                   | Enable Link Time Optimization                | `OFF`                |
| `-DVALGRIND_SUPPRESSIONS=PATH`   | Set a suppressions file for valgrind         | empty                |
| `-DINTEGER8=ON/OFF`              | Build with 64-bit integers and 64 bit BLAS   | `OFF`                |
| `-DFORTRAN_BOUND_CHECK=ON/OFF`   | Check the array bounds at runtime (slow)     | `OFF`                |
| `-DFORTRAN_SANITIZE=ON/OFF`      | Check the memory accesses at runtime (slow)  | `OFF`                |
| `-DDOC=ON/OFF`                   | Build the documentation                      | `OFF`                |
| `-DMATLAB=ON/OFF`                | Enable MATLAB / Octave interface building    | `OFF`                |
| `-DEXAMPLES=ON/OFF`              | Build the tests and examples                 | `ON`                 |
| `-DMEPACK_LUA_CONFIG=PATH`       | Path for an alternative configuration file   | empty                |
| `-DRECSY=PATH`                   | Path to a precompiled library containing     | empty                |
|                                  | the RECSY library.                           | empty                |
| `-DGL705=PATH`                   | Path to a precompiled library containing     | empty                |
|                                  | the code from Algorithm 705 as library.      | empty                |
| `-DCMAKE_INSTALL_MODULEDIR=PATH` | Custom destination directory for Fortran     | /usr/include/mepack/ |
|                                  | modules.                                     |                      |

The examples/tests are only built if `BUILD_DYNAMIC=ON`, otherwise they can not
be used.

The test suite is executed by calling
```shell
make -C build-dir test
```
If `valgrind` was found during the configuration procedure, all tets can be
executed with the help of `valgrind` by executing:
```shell
cd build-dir
ctest -D NightlyMemoryCheck
```

The documentation is built in `build-dir/doc` by calling:
```shell
make -C build-dir doc
```

## MATLAB/Octave Interface

The MATLAB/Octave interface is built if `-DMATLAB=ON` is set during
configuration. The build process can be adjusted using the following CMake
options:

|  Option                       | Description                                       |  Default               |
|-------------------------------|---------------------------------------------------|------------------------|
| `-DMEXOCT_MATLAB_DOC=ON/OFF`  | Extract the help texts for MATLAB from Mex-Files, | `OFF`                  |
| `-DMEXOCT_MATLAB=ON/OFF`      | Search for MATLAB                                 | `ON`                   |
| `-DMEXOCT_OCTAVE=ON/OFF`      | Search for Octave                                 | `ON`                   |
| `-DMEXOCT_MATLAB_ROOT=PATH`   | Path of the MATLAB installation                   | empty                  |
| `-DMEXOCT_OCTAVE_CONFIG=PATH` | Full path of the `octave-config` tool             | empty                  |
| `-DMEXOCT_LINK_STATIC=ON/OFF` | Static linking of MEPACK to the MEX/OCT files     | `OFF`                  |
| `-DMEXOCT_SKIP_RPATH=1/0`     | Skip the RPATH setting in MEX/OCT files           | `0`                    |

After compiling, the build directory contains a `matlab/matlab` folder including
the mex files to be used with MATLAB and a `matlab/octave` folder, which
contains the octfiles for GNU Octave. Furthermore, a distribution archive
```
mepack-<MEPACK_VERSION>-matlab-<MATLAB_VERSION>-<ARCH>-dist.{tar.gz,zip}
```
for MATLAB and
```
mepack-<MEPACK_VERSION>-octave-<OCTAVE_VERSION>-<ARCH>-dist.tar.gz
```
for GNU Octave is created in the root of the build directory. These files can be
used to ship the compiled MEX/OCT files to othersystems.

The test suite is also added to the `test`
target of the main makefile. In this way
```shell
make -C build-dir test
```
will also test the MATLAB/Octave interface. In order to test the MATLAB/Octave
interface separately run
```shell
cd build-dir/matlab
make test
```

By default (`MEXOCT_SKIP_RPATH=0`) the MEX and OCT files have set its `RPATH` to
the $ORIGIN and thus will search for the `libmepack.so` file in the current
directory. In this case, the `libmepack.so` is copied to the `matlab` or
`octave` directory. In this case, the `libmepack.so` file will also be added to
the distribution archives. If `MEXOCT_LINK_STATIC` is enabled, no `RPATH` is
required since MEPACK is statically linked to the MEX/OCT files, resulting in
larger files. The `libmepack.so` file is not part of the distribution file in
this case.

**Attention:** The MATLAB part of the interface requires `-DINTEGER8=ON` in most
cases, since MATLAB uses 64-bit integers everywhere. In the case of GNU Octave
the setting `-DINTEGER8=OFF` is required in most cases, since Octave uses 64-bit
integers internally, but uses  BLAS/LAPACK/other Fortran codes with 32-bit
integers.

**Attention:** Enable the static linkage with `MEXOCT_LINK_STATIC` also
requires to build the static library with `BUILD_STATIC=ON`.

**Warning:** The MATLAB interface is implemented using C++ template meta
programming and thus the compiling the code leads to a huge memory consumption.
If MEPACK is build in parallel, e.g. `make -j 4`, this can crash your computer.
For reason the MATLAB interface should be compiled with one job per 8 GB of free
memory. If your are not sure if this is fulfilled, consider compiling the code
with only one job, i.e. `make -j 1`.

### MATLAB Interface on Microsoft Windows

Since MATLAB and its shipped libraries behave slightly different compared to the
Linux/Unix world, building the MEPACK MATLAB interface under Microsoft Windows
requires a special build configuration. Using MEPACK the MSYS2 MingW64
environment using either the `MSVCRT` or the `UCRT` runtime is required.
Therefor, the requirements from above are required to installed beforehand.

The MATLAB interface is built from the MINGW64 command line using:
```shell
cmake -DLTO=OFF -DINTEGER8=ON -DBUILD_STATIC=ON -DBUILD_DYNAMIC=OFF \
      -DMATLAB=ON -DMEXOCT_OCTAVE=OFF -DMEXOCT_LINK_STATIC=ON \
      -DCMAKE_Fortran_FLAGS="-fno-underscoring" -DMEXOCT_MATLAB_DOC=OFF \
      -S . -B build-matlab
cmake --build build-matlab
```
Afterwards, the distribution can be found as zip-file in the `build-dir`
directory.

Building for GNU Octave on Windows is not possible at the moment.

## RECSY and Algorithm 705

The RECSY library and the reference implementation of Algorithm 705 cannot be
shipped due to unclear licenses. If you want to compile the examples with
support for these codes, you have to provide a shared or static library
containing the codes yourself. Perferably, RECSY and Algorithm 705 should be
compiled as static library with the `-fPIC` flag set during compilation. Then
they can be enabled in the during the build via
```
 cmake ....... -DRECSY=/path/to/librecsy.a -DGL705=/path/to/libgl705.a
```

