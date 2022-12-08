# Tuning                                                               {#tuning}

MEPACK can be tuned at runtime. Therefore, all computation related parameters
can be adjusted at runtime, see @ref options for details. Another variant is
to provide a configuration file. This configuration file is a LUA script, which
contains a set of functions. Each function gets the problem dimension as input
and returns the optimal value for the desired parameter as return value.

By default, the configuration file `src/defaulf_config.lua` is compiled into the
library. If a user supplied configuration file can be used in two ways:

* The environment variable `MEPACK_LUA_CONFIG` is set to an alternative
  configuration file.
* Set the `-DMEPACK_LUA_CONFIG=PATH` option at compile time to the path of the new
  configuration.

Configuration File Template
---------------------------
The following template includes all parameters and the corresponding function
used by MEPACK at the moment. All function arguments and return values are
integer values. The input `m` represents the number of rows of the solution and
the right-hand side. The input `n` denotes the number of columns of the solution
and the right-hand side. If the algorithm solves a "symmetric equation", like
the Lyapunov or the Stein equation, only `m` is used. The parameter only effect
the solvers with (quasi) triangular coefficient matrices, since the solver for
equations with general coefficient matrices rely on them.

```lua
--
-- Outer block size for the two stage schemes. These values are used by the
-- _2STAGE routines.
--

-- TGCSYLV in double precision
function tgcsylv_double_2stage(m,n)
    return 1024
end

-- TGCSYLV_DUAL in double precision
function tgcsylv_dual_double_2stage(m,n)
    return 1024
end

-- TGCSYLV in single precision
function tgcsylv_single_2stage(m,n)
    return 1024
end

-- TGCSYLV_DUAL in single precision
function tgcsylv_dual_single_2stage(m,n)
    return 1024
end

-- TGSYLV in double precision
function tgsylv_double_2stage(m,n)
    return 1024
end

-- TGSYLV in single precision
function tgsylv_single_2stage(m,n)
    return 1024
end

-- TRSYLV in double precision
function trsylv_double_2stage(m,n)
    return 1024
end

-- TRSYLV in single precision
function trsylv_single_2stage(m,n)
    return 1024
end

-- TRSYLV2 in double precision
function trsylv2_double_2stage(m,n)
    return 1024
end

-- TRSYLV2 in single precision
function trsylv2_single_2stage(m,n)
    return 1024
end

-- TGLYAP for double precision
function tglyap_double_2stage(m)
    return 1024
end

-- TGLYAP for single precision precision
function tglyap_single_2stage(m)
    return 1024
end

-- TRLYAP for double precision
function trlyap_double_2stage(m)
    return 1024
end

-- TRLYAP for single precision precision
function trlyap_single_2stage(m)
    return 1024
end

-- TGSTEIN for double precision
function tgstein_double_2stage(m)
    return 1024
end

-- TGSTEIN for single precision precision
function tgstein_single_2stage(m)
    return 1024
end

--
-- Block-sizes for the Level-3 and DAG schemes. These values are used by
-- routines with the names '_L3' or '_DAG'. The function with the suffix '_mb'
-- returns the block-size with respect to the number of rows of the solution.
-- The function suffixed with '_nb' returns the block-size with respect to the
-- number of columns of the solution. In case of symmetric equations, like
-- the Lypanov or the Stein equation. only the function with the suffix '_mb'
-- is used.

-- MB for TGCSYLV in double precision
function tgcsylv_double_mb(m,n)
    return 32
end

-- NB for TGCSYLV in double precision
function tgcsylv_double_nb(m,n)
    return 32
end

-- MB for TGCSYLV_DUAL in double precision
function tgcsylv_dual_double_mb(m,n)
    return 32
end

-- NB for TGCSYLV_DUAL in double precision
function tgcsylv_dual_double_nb(m,n)
    return 32
end

-- MB for TGCSYLV in single precision
function tgcsylv_single_mb(m,n)
    return 32
end

-- NB for TGCSYLV in single precision
function tgcsylv_single_nb(m,n)
    return 32
end

-- MB for TGCSYLV_DUAL in single precision
function tgcsylv_dual_single_mb(m,n)
    return 32
end

-- NB for TGCSYLV_DUAL in single precision
function tgcsylv_dual_single_nb(m,n)
    return 32
end

-- MB for TGSYLV in double precision
function tgsylv_double_mb(m,n)
    return 32
end


-- NB for TGSYLV in double precision
function tgsylv_double_nb(m,n)
    return 32
end

-- MB for TGSYLV in single precision
function tgsylv_single_mb(m,n)
    return 32
end

-- NB for TGSYLV in single precision
function tgsylv_single_nb(m,n)
    return 32
end

-- MB for TRSYLV in double precision
function trsylv_double_mb(m,n)
    return 32
end


-- NB for TRSYLV in double precision
function trsylv_double_nb(m,n)
    return 32
end

-- MB for TRSYLV in single precision
function trsylv_single_mb(m,n)
    return 32
end

-- NB for TRSYLV in single precision
function trsylv_single_nb(m,n)
    return 32
end

-- MB for TRSYLV2 in double precision
function trsylv2_double_mb(m,n)
    return 32
end

-- NB for TRSYLV2 in double precision
function trsylv2_double_nb(m,n)
    return 32
end

-- MB for TRSYLV2 in single precision
function trsylv2_single_mb(m,n)
    return 32
end

-- NB for TRSYLV2 in single precision
function trsylv2_single_nb(m,n)
    return 32
end

-- MB for TGLYAP in double precision
function tglyap_double_mb(m)
    return 32
end

-- MB for TGLYAP in single precision
function tglyap_single_mb(m)
    return 32
end

-- MB for TRLYAP in double precision
function trlyap_double_mb(m)
    return 32
end

-- MB for TRLYAP in single precision
function trlyap_single_mb(m)
    return 32
end

-- MB for TGSTEIN in double precision
function tgstein_double_mb(m)
    return 32
end

-- MB for TGSTEIN in single precision
function tgstein_single_mb(m)
    return 32
end

-- MB for TRSTEIN in double precision
function trstein_double_mb(m)
    return 32
end

-- MB for TRSTEIN in single precision
function trstein_single_mb(m)
    return 32
end
```

Obtaining Optimal values
------------------------
Optimal block-sizes can be obtained by running benchmarks. Therefore, the
examples from `examples/triangular` can be used. For example, the optimal block
size `MB` for the standard Lyapunov equation can be obtained by executing
```
$ ./examples/triangular/benchmark_trlyap --solver=1 --rows=1000:1000:5000 --mb=32:32:128
# HDF5 Store Path: ./ (set with MEPACK_HDF_PATH)
# Command Line: --solver=1 --rows=1000:1000:5000 --mb=32:32:128
# RUNS:  5
# Number of Matrices: 1
# Rows: 1000 (1000:1000:5000)
# TRANSA: N
# Block Alignment: YES
# MACHINE DEFAULT Config:
# Solver: LEVEL3 - LEVEL2: LOCAL COPY fixed maximum size with alignment
#
#  M   MB  Wall-Time     CPU-Time       Ratio    Forward-Err
 1000  32  1.47723e-01  5.73458e-01  3.88198e+00  1.02451e-12
 1000  64  2.02979e-01  6.52629e-01  3.21526e+00  7.70460e-13
 1000  96  1.96380e-01  7.09709e-01  3.61397e+00  7.76942e-13
 1000 128  2.20266e-01  8.20464e-01  3.72487e+00  7.80125e-13
 2000  32  1.52258e+00  5.40316e+00  3.54869e+00  2.53205e-12
 2000  64  1.37894e+00  5.08727e+00  3.68926e+00  3.81413e-12
 2000  96  1.50218e+00  5.38654e+00  3.58582e+00  6.89304e-12
 2000 128  1.42162e+00  5.25848e+00  3.69895e+00  3.36752e-12
 3000  32  4.94669e+00  1.81289e+01  3.66485e+00  6.80014e-12
 3000  64  5.09770e+00  1.77916e+01  3.49011e+00  5.05035e-12
 3000  96  4.85457e+00  1.75745e+01  3.62020e+00  6.99595e-12
 3000 128  6.37335e+00  2.37100e+01  3.72018e+00  5.60129e-12
 4000  32  1.49584e+01  5.41659e+01  3.62110e+00  3.11764e-12
 4000  64  1.32034e+01  4.91697e+01  3.72402e+00  3.50054e-12
 4000  96  1.55753e+01  5.29380e+01  3.39884e+00  3.79209e-12
 4000 128  1.45271e+01  5.27828e+01  3.63339e+00  2.42395e-12
 5000  32  2.87573e+01  1.03002e+02  3.58176e+00  4.27677e-12
 5000  64  2.20910e+01  8.00796e+01  3.62499e+00  5.20105e-12
 5000  96  1.84365e+01  6.93392e+01  3.76097e+00  4.36991e-12
 5000 128  2.49911e+01  9.18636e+01  3.67585e+00  5.64129e-12
```
Selecting the block-sizes with the minimal runtimes from the output and
interpolating between them gives the following `trlyap_double_mb` function on an
Intel Celeron N3450 with OpenBLAS 0.3.8 in pthread-mode:
```
function trlyap_double_mb(m)
    if ( m < 1500 ) then
        return 32
    elseif ( 1500 <= m and m < 2500 ) then
        return 64
    elseif ( 2500 <= m and m < 3500 ) then
        return 96
    elseif ( 3500 <= m and m < 4500 ) then
        return 64
    else
        return 96
    end
end
```

Predefined Configuration Files
------------------------------

We provide a some preconfigured tuning files inside the `src/config/` directory.
At the moment, we provide them for the following systems:

 * Dual Socket Intel Xeon Silver Edition 4110 (2x 8 cores), Intel
   Parallel Studio 2018 (icc, ifort, mkl)
   (`src/config/intel-xeon-silver-4110-parallel-studio-2018.lua`)
 * Dual Socket Intel Xeon Haswell E5-2640 v3 (2x 8 cores), Intel
   Parallel Studio 2018 (icc, ifort, mkl)
   (`src/config/intel-xeon-e5-2640v3-parallel-studio-2018.lua`)
 * Dual Socket IBM Power8 (2x 10 cores), IBM XLC 16.1, IBM XLF 16.1,
   IBM ESSL 6.3
   (`src/config/ibm-power8-xlf-essl.lua`)



Solver Selection
----------------
MEPACK provides a huge set of Level-3, Level-2 and DAG accelerated solvers. By
default the level-3 solvers use the level-2 solvers with aligned local copies.
The DAG accelerated solvers use the same selection. In the routines for solving
the equations with general coefficient matrices, the level-3 triangular solver
is the default. Everything can be changed with the help of the routines from
[Options - Level 2 Solvers](@ref isolver) and
[Options - Frontend Solvers](@ref frontend_solver). If the GNU Fortran compiler
is used, using the level-2 solvers with aligned local copies is no longer
beneficial. In this case one can select the level-2 with local copies (but
without alignment) to achieve the same or in some cases a better performance.
For such reasons, all the differently optimized solvers are included in MEPACK.
In this way, one can evaluate different solver approaches on newly emerge
hardware platforms to check for the best performance easily.

Recommendations
---------------
Since there are too many different possibilities to solve a given equation with
MEPACK, there are some basic rules, when which solver-type is beneficial. In all
cases 'dimension' is meant in the context of the size of the right-hand side.

 * If the problem is small, e.g. m,n <= 128, use a level-2.
 * If the problem is small to medium sized, e.g. m,n <= 1000, use a level-3
   solver.
 * If the problem is large, e.g. m,n = 1000 ... 5000, and only a few CPU cores
   are available, use a level-3 solver as well.
 * If the problem is large, e.g. m,n = 1000 ... 5000, and many CPU cores are
   available, use a DAG accelerated solver.
 * If the problem is huge, e.g. m,n >= 5000, use the 2-stage solver.

These recommendations are only a personal experience on the systems of the
author had available during development. The situation can change on different
hardware dramatically.
