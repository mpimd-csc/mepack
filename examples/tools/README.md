Tools
=====

This directory contains some tools helping with the examples. Each tool supports
the `-h` or `--help` option on the command line to see the acutal command line
parameters that can be used to specify its work,

- `generate_evp.c`: Generates random general martrices and computes their Schur
  decompositions. The results are stored in `${MEPACK__HDF_PATH}/qr`. The naming
  scheme of the individual files are consistent to the way the matrix equation
  benchmarks will cache their Schur decompositions. In this way, this tool can
  precompute input data for the solvers with triangular coefficient matrices.

- `generate_evp.c`: Generates random general martrix pairs and computes their
  generalized Schur decompositions. The results are stored in
  `${MEPACK__HDF_PATH}/qz`. The naming scheme of the individual files are
  consistent to the way the matrix equation benchmarks will cache their
  generalized Schur decompositions. In this way, this tool can
  precompute input data for the generalized solvers with triangular coefficient
  matrices.

