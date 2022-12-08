Examples
========

This directory contains a set of benchmark and examples codes for MEPACK.  All
examples try to cache intermediate Schur Decompositions to accelerate the
benchmark and test runs. Therefore, the `MEPACK_HDF_PATH` environment variable
can be used to specify an alternative path. If the variable is not set, the
current directory is used. The cached matrices are stored in HDF5 files. The can
be precomputed with the `generate_evp` and `generate_gevp` tools from the
`tools` directory. The HDF5 files are **not** MATLAB mat-file compatible.

The examples consist of the following subdirectories:

 - [`benchmark`](benchmark/README.md):  Overall benchmarks for the solvers with
   triangular coefficient matrices.
 - [`benchmark_lib`](benchmark_lib/README.md): Auxiliary routines for managing
   the benchmarks. This contains no executables and only provides a helper
   library.
 - [`frontend`](frontend/README.md): Tests and examples for the frontend
   routines, i.e. the interfaces, which can be used with general coefficient
   matrices.
 - [`misc`](misc/README.md): TBA
 - [`refine`](refine/README.md): Tests and examples for the iterative refinement
   routines.
 - [`tools`](tools/README.md): Various tools supporting the tests and examples.
 - [`triangular`](triangular/README.md):  Tests and examples for the solvers
   with triangular coefficient matrices.
 - [`data`](data/README.md): Input data for some examples.

Tuning File Example
-------------------

The `mepack_lua_config.lua` file contains a minimal tuning script that can be
used as examples or as base for own optimizations. The default configuration can
be changed via setting the `MEPACK_LUA_CONFIG` environment variable. For details
about the tuning see [Tuning](../doc/tuning.md).

