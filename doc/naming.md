# Naming Scheme for the Computational Routines                         {#naming}

The computational subroutines, that are not part of a Fortran 90 module, follow
the Fortran 77 old style calling conventions. The routines have the following
naming scheme:
```
    SUBROUTINE AAA_FFEE_WWW(ARGUMENTS...)

```
Thereby, the three components `AAA`, `EEE`, and `WWW` have to following
meanings:

* `AAA` specifies the data type/precision for the numerical input values. This
  is either `SLA` for **Single Precision Linear Algebra** or `DLA` for **Double
  Precision Linear Algebra**. A few routines begin with different prefixes.
  These routines are only for internal usage and not meant to be called from
  outside. Making them visible is necessary to fix some issues with different
  Fortran compilers.
* `FF` are two characters specifying the coefficient matrices. This can be one
  the following combinations:
  * `TR`, the coefficient matrices are (quasi) upper triangular matrices.
  * `TG`, the coefficient matrices are matrix pairs with (quasi) upper triangular
    matrices.
  * `GE`, the coefficient matrices are general matrices.
  * `GG`, the coefficient matrices are general matrices, that appear in matrix
    pairs.
* `EE` specifies the equations, this can be
  * `SYLV`, for Sylvester equations
  * `SYLV2`, for discrete time Sylvester equations
  * `LYAP`, for Lyapunov equations
  * `STEIN`, for Stein equations
  * `CSYLV`, for the coupled generalized Sylvester equation
  * `CSYLV_DUAL`, for the dual coupled generalized Sylvester equation
* `WWW` (optionally) specifies the details about the solver. The following
  strings can be part of this description:

  * `L3`, for a level-3 solver
  * `L2`, for a level-2 solver
  * `DAG`, for a parallel, graph-scheduled, solver
  * `RECURSIVE`, for a recursive solver
  * `REFINE`, for a solver which uses iterative refinement

  Additionally, the `WWW` contains hints about the optimizations in the
  routines. For details, see the documentation of the individual routines.

