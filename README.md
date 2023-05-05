Matrix Equations PACKage (MEPACK)
=================================

**Version:** 1.0.4

Copyright 2017-2023 by Martin Köhler, MPI-Magdeburg

The software is licensed under GPLv3 or later. See `LICENSE` for details.

MEPACK is a Fortran software library for the solution of dense Sylvester-like
matrix equations. This includes the following equations:

  * Standard Sylvester Equation (`SYLV`) $`AX\pm XB = Y`$
  * Standard Lyapunov Equation (`LYAP`) $`AX + XA^T = Y`$
  * Standard Stein (Discrete-Time Lyapunov) Equation (`STEIN`) $`AXA^T - X = Y`$
  * Discrete-Time Sylvester Equation (`SYLV2`) $`AXB - X = Y`$
  * Generalized Sylvester Equation (`GSYLV`) $`AXB \pm CXD = Y`$
  * Generalized Lyapunov Equation (`GLYAP`) $`AXB^{T} + BXA^T = Y`$
  * Generalized Stein (Discrete-Time Lyapunov) Equation (`GSTEIN`)
    $`AXA^T - EXE^T = Y`$
  * Generalized Coupled Sylvester Equation (`CSYLV`)
    $`(AR \pm LB=E,CR \pm LD = F)`$ and
    its dual counterpart (`CSYLV_DUAL`) $`(AR + CL = E, \pm RB \pm LD =F)`$.

The implemented solvers use a forward/backward substitution scheme on matrix
equations with (quasi) triangular coefficient matrices. The algorithms are
implemented in different ways, including the following features:

 * Level-3 BLAS enabled solvers with adjustable block/tile sizes.
 * OpenMP 4 accelerated Directed-Acyclic-Graph scheduled solvers with focus on
   current multicore CPUs
 * Optimized BLAS level-2 solvers for problems with right-hand sides smaller
   than 128x128.
 * Iterative refinement solvers to improve the accuracy if the (generalized)
   Schur decompositions are perturbed.

Besides that, the library implements some other solvers for performance
comparison:

  * The recursive blocking approach by Jonsson and Kågström,
  * The Gardiner-Laub approach for the generalized Sylvester equation.

  References for both techniques can be found below.


Installation
------------
See [Installation](doc/install.md).

API Documentation
-----------------
All routines in the `src/` directory are documented in a BLAS/LAPACK style. The
documentation can be extracted using `doxygen` during the build process.
See [Installation](doc/install.md).



Configuration at Runtime
------------------------
Some parameters, like the block sizes of the level-3, OpenMP 4, and the two
stage algorithm, are hardware dependent and need to be tuned for an optimal
runtime behavior. They can either be hard-coded or configured dynamically at
runtime. The dynamic selection is used, if a block size or a parameter is set to
zero at runtime. This is the default situation. The dynamic configuration is
given as a LUA script, either at compile-time in the  `src/default_config.lua`
file or at runtime by passing a configuration file in the the
`MEPACK_LUA_CONFIG` environment variable. Details about this feature can be
found in  [Tuning](doc/tuning.md).

In some cases, MEPACK  prints some debug information to `stderr`. This can be
controlled by setting the `MEPACK_VERBOSE` environment variable. A non-zero
value means, that MEPACK shows some about the runtime environment, like the
loaded configuration file. The computational routines are not effected by these
outputs. Error while reading the configuration file will
always be printed to the screen independent of the value of `MEPACK_VERBOSE`.

Using MEPACK your project
-------------------------

The MEPACK project exports CMake package files and pkg-config files to make
MEPACK usable for other projects. The package files are located in the library
directory in the installation prefix.

In order to use MEPACK with CMake it is included via

```cmake
find_package(MEPACK REQUIRED)
# or in case of 64 bit integers
# find_package(MEPEACK64 REQUIRED)
...
target_link_libraries(
  yourtarget
  PRIVATE
  MEPACK::mepack
)
```

If the static version of MEPACK should be used, replace `MEPACK::mepack` with
`MEPACK::mepack_static`.

To make the installed MEPACK project discoverable, add the MEPACK install prefix
to the `CMAKE_PREFIX_PATH` variable. The usual install location of the
package files is `PREFIX/lib/cmake/mepack` or `PREFIX/lib/cmake/mepack64`.

For non-CMake build systems (like make) you can use the exported pkg-config
file by setting `PKG_CONFIG_PATH` to include the directory containing the
exported pc-file. The usual install location of the pc-file is
`PREFIX/lib/pkgconfig`. In `make`, you can obtain the required compile and link
arguments with

```make
MEPACK_CFLAGS := $(shell pkg-config --cflags mepack)
MEPACK_LIBS := $(shell pkg-config --libs mepack)
```

or

```make
MEPACK_CFLAGS := $(shell pkg-config --cflags mepack64)
MEPACK_LIBS := $(shell pkg-config --libs mepack64)
```

if you require 64 bit integers.

Citation
---------
If you are using MEPACK inside your project, please cite at least one of the
following references:
```
M. Köhler. 2021. Approximate Solution of Non-Symmetric Generalized Eigenvalue
Problems and Linear Matrix Equation on HPC Platforms. Dissertation. Logos
Verlag Berlin, Magdeburg, Germany.
```
or
```
M. Köhler. 2022, Matrix Equations PACKage -- A Fortran library for the solution
of Sylvester-like Matrix equations,
Zendo, DOI: 10.5281/zenodo.7554327
```

Examples
--------

The `examples` directory provide a huge set of demonstration codes how to use
MEPACK. Furthermore, these examples are also constructed to serve as tests for
the library. Thus, most of them are executed when `make test` is run.

For details about the examples see [examples/README.md](examples/README.md).

External Codes included in the Examples
---------------------------------------
The examples contain a set of foreign codes for benchmark purpose.
These codes are subject to the copyrights issued by their original authors.
At the moment this effects the `examples/slicot` directory, which provides
some routines from SLICOT (https://www.slicot.org,
https://github.com/SLICOT/SLICOT-Reference, BSD-3 License):
`MB01RD`, `MB01RW`, `MB02UV`, `MB02UU`, `SB03MD`, `SB03MV`, `SB03MW`, `SB03MX`,
`SB03MY`, `SB04PX`, `SB04PY`, `SG03AD`, `SG03AX`, and `SG03AY`.

Furthermore, the codes of the following publications can be used for comparison
in the examples:
```
Judith D. Gardiner, Matthew R. Wette, Alan J. Laub, James J. Amato, and
Cleve B. Moler. 1992. Algorithm 705; a FORTRAN-77 software package for solving
the Sylvester matrix equation AXB_T_ + CXD_T_ = E. _ACM Trans. Math. Softw._ 18,
2 (June 1992), 232–238.
DOI:https://doi.org/10.1145/146847.146930
```
and
```
Isak Jonsson and Bo Kågström. 2002. Recursive blocked algorithms for solving triangular
systems—Part I: one-sided and coupled Sylvester-type matrix equations. _ACM Trans. Math.
Softw._ 28, 4 (December 2002), 392–415.
DOI:https://doi.org/10.1145/592843.592845

Isak Jonsson and Bo Kågström. 2002. Recursive blocked algorithms for solving triangular
systems—Part II: two-sided and generalized Sylvester and Lyapunov matrix equations. _ACM
Trans. Math. Softw._ 28, 4 (December 2002), 416–435.
DOI:https://doi.org/10.1145/592843.592846
```

Due to unclear licensing, these codes need to be provided by the user.
See [Installation](doc/install.md) for details.



Known Issues
------------
  * If the OpenBLAS 0.3.8 with PTHREADs is used, `OMP_NUM_THREADS=1` needs to be
    set. The pthread scheduler of OpenBLAS conflicts with the one from OpenMP.
  * The MATLAB interface ( not the Octave Interface ) cannot use OpenMP 4
    accelerated solvers due to deadlocks in MATLAB's OpenMP translation
    interface.
  * If doxygen 1.8.17 is installed on the system, like on Ubuntu 20.04, the
    documentation includes MathJax to render the formulas. This is due to a bug
    in 1.8.17 with some versions of Ghostscript, especially Ghostscript 9.50.-

