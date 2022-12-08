Miscellaneous Examples
======================

This directory contains examples and tests for the frontend routines, i.e.
solving equations with general coefficient matrices. The examples are used as
test codes as well.  All examples feature the `--help` command line options to
view the actual parameters that can be used to adjust the example,

The following miscellaneous codes are available:
- `benchmark_eigenvalue_reorder`: Check the influence of sorted eigenvalues on
                                  the performance of the solvers with local
                                  aligned copies for the generalized Sylvester
                                  equation.
- `benchmark_dag`: Benchmark the level-3 and the DAG solver for the generalized
                   Sylvester Equation.
- `benchmark_recsy`: Benchmark the standard RECSY and the parallel RECSY
                     implementation. **This is only available if the examples
                     are compiled with RECSY support. Otherwise, this is not
                     build**.
- `benchmark_row_columnwise`: Benchmark the influence of the row- and column-
                              wise computation of the solution for the
                              generalized Sylvester equation.

Generalized Lyapunov Solvers
----------------------------

The following codes are solvers for the generalized Lyapunov equation that can
read the coefficient matrices from a MATLAB (v7.3) file. The following datasets
have to exist in the MATLAB file:

 - `A`: co
Miscellaneous Examples
======================

This directory contains examples and tests for the frontend routines, i.e.
solving equations with general coefficient matrices. The examples are used as
test codes as well.  All examples feature the `--help` command line options to
view the actual parameters that can be used to adjust the example,

The following miscellaneous codes are available:
- `benchmark_eigenvalue_reorder`: Check the influence of sorted eigenvalues on
                                  the performance of the solvers with local
                                  aligned copies for the generalized Sylvester
                                  equation.
- `benchmark_dag`: Benchmark the level-3 and the DAG solver for the generalized
                   Sylvester Equation.
- `benchmark_recsy`: Benchmark the standard RECSY and the parallel RECSY
                     implementation. **This is only available if the examples
                     are compiled with RECSY support. Otherwise, this is not
                     build**.
- `benchmark_row_columnwise`: Benchmark the influence of the row- and column-
                              wise computation of the solution for the
                              generalized Sylvester equation.

Generalized Lyapunov Solvers
----------------------------

The following codes are solvers for the generalized Lyapunov equation that can
read the coefficient matrices from a MATLAB (v7.3) file. The following datasets
have to exist in the MATLAB file:

 - `A`: coefficient matrix A
 - `B`: coefficient matrix B
 - `ASCHUR`: generalized Schur form of A
 - `BSCHUR`: generalized Schur form of B
 - `Q`: left-hand side transformation matrix for the generalized Schur form
 - `Z`: right-hand side transformation matrix for the generalized Schur form
 - `TIMES`: vector of length two containing the WALL and the CPU time required
            to compute the generalized Schur decomposition.

All codes are executes via:
```shell
 examples inputfile.mat
```

The following examples exist:

 - `solve_slicot_glyap`: Solve the generalized Lyapunov equation with SLICOT.
 - `solve_mepack_glyap`: Solve the generalized Lyapunov equation with MEPACK.
 - `solve_refine_glyap`: Solve the generalized Lyapunov equation with iterative
                         refinement.


effcient matrix A
 - `B`: coeffcient matrix B
 - `ASCHUR`: generalized Schur form of A
 - `BSCHUR`: generalized Schur form of B
 - `Q`: left-hand side transformation matrix for the generalized Schur form
 - `Z`: right-hand side transformation matrix for the generalized Schur form
 - `TIMES`: vector of length two containing the WALL and the CPU time required
            to compute the generalized Schur decomposition.

All codes are executes via:
```shell
 examples inputfile.mat
```

The following examples exist:

 - `solve_slicot_glyap`: Solve the generalized Lyapunov equation with SLICOT.
 - `solve_mepack_glyap`: Solve the generalized Lyapunov equation with MEPACK.
 - `solve_refine_glyap`: Solve the generalized Lyapunov equation with iterative
                         refinement.


