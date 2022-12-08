Benchmarks
==========

This directory contains comparison benchmarks for the different level-2 solvers.
Each problem will be solved using the standard level-3 solver with different
level-3 solvers, the OpenMP accelerated solver, the level-3/OpenMP two-stage
solver, and RECSY (if provided by the user). Each benchmark features the
`--help` command line options to view the actual parameters that can be used to
adjust the benchmark. All benchmarks only work on matrix equations with
triangular coefficient matrices.

The following benchmarks are available:
- `benchmark_all_tgcsylv`: Benchmark the generalized coupled Sylvester equation
						   (CSYLV)
- `benchmark_all_tgstein`: Benchmark the generalized Stein equation (GSTEIN)
- `benchmark_all_trlyap`: Benchmark the standard Lyapunov equation (LYAP)
- `benchmark_all_trsylv2`: Benchmark the discrete-time Sylvester equations
						   (SYLV2)
- `benchmark_all_tglyap`: Benchmark the generalized Lyapunov equation (GLYAP)
- `benchmark_all_tgsylv`: Benchmark the generalized Sylvester equation (GSYLV)
- `benchmark_all_trstein`: Benchmark the standard Stein equation (STEIN)
- `benchmark_all_trsylv`: Benchmark the standard Sylvester equation (SYLV)


