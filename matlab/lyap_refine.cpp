/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (C) Martin Koehler, 2017-2023
 */


#define MEXOCT_MATLAB_DQ_STRING


#include "mepack.h"
#include "mexoct/interface.hpp"
#include "mepackoptions.hpp"
#include "qtriangular_test.hpp"

using namespace mexoct;

#ifdef MEXOCT_OCTAVE
#define ROWS(a)  (a.rows())
#define COLS(a)  (a.columns())
#define ARG_T      const octave_value &
#else
#define ROWS(a)  (mxGetM(a))
#define COLS(a)  (mxGetN(a))
#define ARG_T    const mxArray *
#endif

template<typename T>
struct lyap_refine_solver {
    constexpr static const char *fn = NULL;
    virtual void solve (char * TRANS, char *GUESS, int M , T * A, int LDA,
        T *X, int LDX, T *Y, int ldy, T *AS, int LDAS, T * Q, int LDQ,
        int *MAXIT, T *TAU, T *CONVLOG, T * WORK, size_t LWORK, int *INFO);
    virtual void blocksize_set(int nb) = 0;
    virtual void blocksize_big_set(int nb) = 0;
    virtual void isolver_set(int is) = 0;
    virtual void fsolver_set(int fs) = 0;
};

template<>
struct lyap_refine_solver<double> {
    constexpr static const char *fn = "DLA_GELYAP_REFINE";

    static void solve (const char * TRANS, const char *GUESS, int M , double * A, int LDA,
        double *X, int LDX, double *Y, int LDY, double *AS, int LDAS, double * Q, int LDQ,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LWORK, int *INFO) {
        mepack_double_gelyap_refine(TRANS, GUESS, M, A, LDA, X, LDX, Y, LDY, AS, LDAS, Q, LDQ, MAXIT, TAU, CONVLOG, WORK, LWORK, INFO);
    }

    static void blocksize_set(int nb) {
        mepack_double_trlyap_blocksize_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_trlyap_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_trlyap_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_double_trlyap_blocksize_2stage_set(nb);
    };

};

template<>
struct lyap_refine_solver<float> {
    constexpr static const char *fn = "SLA_GELYAP_REFINE";

    static void solve (const char * TRANS, const char *GUESS, int M , float * A, int LDA,
        float *X, int LDX, float *Y, int LDY, float *AS, int LDAS, float * Q, int LDQ,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LWORK, int *INFO) {
        mepack_single_gelyap_refine(TRANS, GUESS, M, A, LDA, X, LDX, Y, LDY, AS, LDAS, Q, LDQ, MAXIT, TAU, CONVLOG, WORK, LWORK, INFO);
    }

    static void blocksize_set(int nb) {
        mepack_single_trlyap_blocksize_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_trlyap_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_trlyap_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_single_trlyap_blocksize_2stage_set(nb);
    };
};


template <typename T>
static std::tuple<ArrayType<T>, ArrayType<T>>
LyapDriver(ArrayType<T> A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X, std::string opA = std::string("N"), std::string guess = std::string("N"),
            MepackOptions optset = MepackOptions())
{
    if (A.rows != Y.rows || A.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix and the right hand side have to be of the same dimension size(A) = [ %d %d ] size (Y) = [ %d %d ]", (int) A.rows, (int) A.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (S.rows != Y.rows || S.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix and the right hand side have to be of the same dimension size(S) = [ %d %d ] size (Y) = [ %d %d ]", (int) S.rows, (int) S.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (X.rows != Y.rows || X.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The initial guess matrix and the right hand side have to be of the same dimension size(X0) = [ %d %d ] size (Y) = [ %d %d ]", (int) X.rows, (int) X.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (Q.rows != Y.rows || Q.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The unitary transformation matrix and the right hand side have to be of the same dimension size(Q) = [ %d %d ] size (Y) = [ %d %d ]", (int) Q.rows, (int) Q.columns,
                (int) Y.rows, (int) Y.columns);
    }


    if ( optset.getMB() >= 0 ) {
        lyap_refine_solver<T>::blocksize_set(optset.getMB());
    }
    if ( optset.getISolver() >= 0 ) {
        lyap_refine_solver<T>::isolver_set(optset.getISolver());
    }
    if ( optset.getFSolver() >= 0 ) {
        lyap_refine_solver<T>::fsolver_set(optset.getFSolver());
    }

    if ( optset.getMaxit() < 1) {
        mexoct_error("MEPACK:MaxitInvalid", "The maximum number of iterations is invalid. maxit = %d", (int) optset.getMaxit());
    }

    if ( optset.getTau() < 0 ) {
        mexoct_error("MEPACK:TolInvalid", "The tolerance is invalid. tol = %20.15e", (double) optset.getTau());
    }

    int maxit = (mexoct_int) optset.getMaxit();
    T  tol = (T) optset.getTau();


    ssize_t auxmem;
    auxmem = mepack_memory_frontend(lyap_refine_solver<T>::fn, (const char *){"F"}, (const char *){"F"}, A.rows, A.columns);
    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    ArrayType<T> convlog_internal(maxit, 1, false);

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    int info = 0;

    lyap_refine_solver<T>::solve(opA.c_str(), guess.c_str(), A.rows, A.ptr(), A.ld, X.ptr(), X.ld, Y.ptr(), Y.ld,
                                 S.ptr(), S.ld, Q.ptr(), Q.ld, &maxit, &tol, convlog_internal.ptr(), auxbuf, auxmem, &info);

    ArrayType<T> convlog( maxit, 1, false);
    for ( auto i = 0; i < maxit; i++) convlog(i, 0) = convlog_internal(i, 0);

    if ( info != 0) {
        mexoct_error("MEPACK:GELYAP_REFINE", "GELYAP_REFINE failed with info = %d", (int) info);
    }
    return std::make_tuple(X, convlog) ;
}


/* Master routines */
template <typename T>
std::tuple<ArrayType<T>>
call_gelyap_op_optset(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{

    bool chk = mepack_mexoct_check_qtriangular(S);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix S is not (quasi) upper triangular.");
    }

    ArrayType<T> X(A.rows, A.rows);
    auto ret = LyapDriver(A, S, Q, Y, X, opA, "N", optset);
    return std::make_tuple(std::get<0>(ret));
}

template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gelyap_op_optset_convlog(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{

    bool chk = mepack_mexoct_check_qtriangular(S);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix S is not (quasi) upper triangular.");
    }

    ArrayType<T> X(A.rows, A.rows);
    auto ret = LyapDriver(A, S, Q, Y, X, opA, "N", optset);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret));
}

template <typename T>
std::tuple<ArrayType<T>>
call_gelyap_x0_op_optset(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{

    bool chk = mepack_mexoct_check_qtriangular(S);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix S is not (quasi) upper triangular.");
    }

    auto ret = LyapDriver(A, S, Q, Y, X, opA, "I", optset);
    return std::make_tuple(std::get<0>(ret));
}

template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gelyap_x0_op_optset_convlog(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{

    bool chk = mepack_mexoct_check_qtriangular(S);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix S is not (quasi) upper triangular.");
    }

    auto ret = LyapDriver(A, S, Q, Y, X, opA, "I", optset);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret));
}



/* Shorter Interfaces */

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y)
{
    return call_gelyap_op_optset(A, S, Q, Y);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gelyap_convlog(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y)
{
    return call_gelyap_op_optset_convlog(A, S, Q, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_x0(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X)
{
    return call_gelyap_x0_op_optset(A, S, Q, Y, X);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gelyap_x0_convlog(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X)
{
    return call_gelyap_x0_op_optset_convlog(A, S, Q, Y, X);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_op(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, std::string op)
{
    return call_gelyap_op_optset(A, S, Q, Y, op);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gelyap_op_convlog(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, std::string op)
{
    return call_gelyap_op_optset_convlog(A, S, Q, Y, op);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_x0_op(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X, std::string op)
{
    return call_gelyap_x0_op_optset(A, S, Q, Y, X, op);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gelyap_x0_op_convlog(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X, std::string op)
{
    return call_gelyap_x0_op_optset_convlog(A, S, Q, Y, X, op);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_opt(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, MepackOptions opt)
{
    return call_gelyap_op_optset(A, S, Q, Y, "N",  opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gelyap_opt_convlog(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, MepackOptions opt)
{
    return call_gelyap_op_optset_convlog(A, S, Q, Y, "N", opt);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_x0_opt(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X, MepackOptions opt)
{
    return call_gelyap_x0_op_optset(A, S, Q, Y, X, "N", opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gelyap_x0_opt_convlog(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X, MepackOptions opt)
{
    return call_gelyap_x0_op_optset_convlog(A, S, Q, Y, X, "N", opt);
}


template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_op_opt(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, std::string op, MepackOptions opt)
{
    return call_gelyap_op_optset(A, S, Q, Y, op, opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gelyap_op_opt_convlog(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, std::string op, MepackOptions opt)
{
    return call_gelyap_op_optset_convlog(A, S, Q, Y, op, opt);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_x0_op_opt(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X, std::string op, MepackOptions opt)
{
    return call_gelyap_x0_op_optset(A, S, Q, Y, X, op, opt );
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gelyap_x0_op_opt_convlog(ArrayType<T>  A, ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y, ArrayType<T> X, std::string op, MepackOptions opt)
{
    return call_gelyap_x0_op_optset_convlog(A, S, Q, Y, X, op, opt);
}








MEXOCT_ENTRY(mepack_lyap_refine,
R"===(
X = mepack_lyap_refine(A, S, Q, Y)
X = mepack_lyap_refine(A, S, Q, Y, X0)
X = mepack_lyap_refine(A, S, Q, Y, Op)
X = mepack_lyap_refine(A, S, Q, Y, X0, Op)
X = mepack_lyap_refine(A, S, Q, Y, OptSet)
X = mepack_lyap_refine(A, S, Q, Y, X0, OptSet)
X = mepack_lyap_refine(A, S, Q, Y, Op, OptSet)
X = mepack_lyap_refine(A, S, Q, Y, X0, Op, OptSet)

Input:
    A      -- coefficient matrix
    S      -- Schur factor of the coefficient matrix
    Q      -- Unitary transformation of the Schur decomposition
    Y      -- right hand side
    X0     -- Initial Guess
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.
 Return:
    X      -- Solution of the Lyapunov equation

[X, convlog] = mepack_lyap_refine(A, S, Q, Y)
[X, convlog] = mepack_lyap_refine(A, S, Q, Y, X0)
[X, convlog] = mepack_lyap_refine(A, S, Q, Y, Op)
[X, convlog] = mepack_lyap_refine(A, S, Q, Y, X0, Op)
[X, convlog] = mepack_lyap_refine(A, S, Q, Y, OptSet)
[X, convlog] = mepack_lyap_refine(A, S, Q, Y, X0, OptSet)
[X, convlog] = mepack_lyap_refine(A, S, Q, Y, Op, OptSet)
[X, convlog] = mepack_lyap_refine(A, S, Q, Y, X0, Op, OptSet)

Input:
    A      -- coefficient matrix
    S      -- Schur factor of the coefficient matrix
    Q      -- Unitary transformation of the Schur decomposition
    Y      -- right hand side
    X0     -- Initial Guess
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.
 Return:
    X        -- Solution of the Lyapunov equation
    convlog  -- convergence history


The mepack_lyap_refine function solves the Lyapunov equation

   op(A) * X + X * op(A)' = Y,

where A is a general matrix employing the iterative refinement
strategy. The operation op either returns the matrix or
if Op == 'T' its transpose. In this way, the transposed equation
can be solved without transposing the coefficient matrix before.
The function requires the matrix A and its Schur decomposition

   [ Q, S ] = schur (A, 'real')

in order to work. The Schur decomposition can also be obtain from
the output of a previous run to mepack_lyap.
The iterative refinement stops once

  || op(A) * X_k + X_k * op(A)' - Y ||_F
  --------------------------------------  < 2 * m * ||A||_F * u * tau
             || Y ||_F

is fulfilled. Thereby m is the order of the matrix A, u is the
machine precision in double or single precision, and tau is a
security factor. Tau can be changed in the OptSet.

The matrices A and Y must be of the same size and Y must be symmetric.

The function uses the level-3 solver D/SLA_GELYAP_REFINE of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_lyap_refine:

   - mb      -- block size of the level-3 algorithm (mb >=2)
   - bignb   -- block size of two-stage solver (bignb >=2)
   - isolver -- internal level-2 solver in the algorithm
       = 0   -- Default/Automatic selection
       = 1   -- Local Copies, with alignments (if supported by the compiler)
       = 2   -- Local Copies, without alignments
       = 3   -- Reordering and Fortran Vectorization
       = 4   -- Classic level-2 solver
       = 5   -- Naive level-2 solver
       = 6   -- Recursive Blocking
   - fsolver -- determine the internal level 3 solver
       = 0   -- automatic or default selection
       = 1   -- standard level 3 solver
       = 2   -- standard level 2 solver
       = 3   -- DAG OpenMP 4 solver
       = 4   -- Two-Stage level 3 solver
       = 5   -- Recursive Blocking Solver
   - maxit   -- maximum number of iterations in the iterative
                refinement, default: 15
   - tau     -- security in the iterative refinement, default: 1.0
   - openmp  -- enable/disable OpenMP support. By default
                the support is enabled GNU Octave and
                disabled in MATLAB


This file file is part of MEPACK, the Matrix Equation Package
by Martin Koehler.
)===")
{

    MEXOCT_INIT();
    mepack_init();

    auto ParamOpString = Parameter<std::string>("Op", "Transpose operation on the coefficient matrix");
    auto ParamOptSet   = Parameter<MepackOptions>("optset", "Options for MEPACK");

    auto ParamDoubleA  = Parameter<ArrayType<double>,Traits::SquareMatrix>("A", "double precision coefficient matrix");
    auto ParamDoubleS  = Parameter<ArrayType<double>,Traits::SquareMatrix>("S", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleQ  = Parameter<ArrayType<double>,Traits::SquareMatrix>("Q", "double precision unitary matrix corresponding to S");
    auto ParamDoubleY  = Parameter<ArrayType<double>,Traits::SquareMatrix>("Y", "double precision right hand side");
    auto ParamDoubleX0 = Parameter<ArrayType<double>,Traits::SquareMatrix>("X0", "double precision initial guess");

    auto ParamFloatA  = Parameter<ArrayType<float>,Traits::SquareMatrix>("As",  "single precision coefficient matrix");
    auto ParamFloatS  = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ss",  "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatQ  = Parameter<ArrayType<float>,Traits::SquareMatrix>("Qs",  "single precision unitary matrix corresponding to Ss");
    auto ParamFloatY  = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ys",  "single precision right hand side");
    auto ParamFloatX0 = Parameter<ArrayType<float>,Traits::SquareMatrix>("X0s", "single precision initial guess");


    auto gelyap_double                         = makeFunctionExt<1>(call_gelyap<double>, "X", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY);
    auto gelyap_double_convlog                 = makeFunctionExt<2>(call_gelyap_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY);

    auto gelyap_double_x0                      = makeFunctionExt<1>(call_gelyap_x0<double>, "X", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamDoubleX0);
    auto gelyap_double_x0_convlog              = makeFunctionExt<2>(call_gelyap_x0_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamDoubleX0);

    auto gelyap_double_op                      = makeFunctionExt<1>(call_gelyap_op<double>, "X", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOpString);
    auto gelyap_double_op_convlog              = makeFunctionExt<2>(call_gelyap_op_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOpString);

    auto gelyap_double_x0_op                   = makeFunctionExt<1>(call_gelyap_x0_op<double>, "X", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamDoubleX0, ParamOpString);
    auto gelyap_double_x0_op_convlog           = makeFunctionExt<2>(call_gelyap_x0_op_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamDoubleX0, ParamOpString);

    auto gelyap_double_opt                      = makeFunctionExt<1>(call_gelyap_opt<double>, "X", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOptSet);
    auto gelyap_double_opt_convlog              = makeFunctionExt<2>(call_gelyap_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOptSet);

    auto gelyap_double_x0_opt                   = makeFunctionExt<1>(call_gelyap_x0_opt<double>, "X", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamDoubleX0, ParamOptSet);
    auto gelyap_double_x0_opt_convlog           = makeFunctionExt<2>(call_gelyap_x0_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamDoubleX0, ParamOptSet);

    auto gelyap_double_op_opt                   = makeFunctionExt<1>(call_gelyap_op_opt<double>, "X", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOpString, ParamOptSet);
    auto gelyap_double_op_opt_convlog              = makeFunctionExt<2>(call_gelyap_op_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOpString, ParamOptSet);

    auto gelyap_double_x0_op_opt                   = makeFunctionExt<1>(call_gelyap_x0_op_opt<double>, "X", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamDoubleX0, ParamOpString, ParamOptSet);
    auto gelyap_double_x0_op_opt_convlog           = makeFunctionExt<2>(call_gelyap_x0_op_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamDoubleX0, ParamOpString, ParamOptSet);

    auto gelyap_float                         = makeFunctionExt<1>(call_gelyap<float>, "X", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY);
    auto gelyap_float_convlog                 = makeFunctionExt<2>(call_gelyap_convlog<float>, "X, convlog", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY);

    auto gelyap_float_x0                      = makeFunctionExt<1>(call_gelyap_x0<float>, "X", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamFloatX0);
    auto gelyap_float_x0_convlog              = makeFunctionExt<2>(call_gelyap_x0_convlog<float>, "X, convlog", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamFloatX0);

    auto gelyap_float_op                      = makeFunctionExt<1>(call_gelyap_op<float>, "X", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamOpString);
    auto gelyap_float_op_convlog              = makeFunctionExt<2>(call_gelyap_op_convlog<float>, "X, convlog", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamOpString);

    auto gelyap_float_x0_op                   = makeFunctionExt<1>(call_gelyap_x0_op<float>, "X", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamFloatX0, ParamOpString);
    auto gelyap_float_x0_op_convlog           = makeFunctionExt<2>(call_gelyap_x0_op_convlog<float>, "X, convlog", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamFloatX0, ParamOpString);

    auto gelyap_float_opt                      = makeFunctionExt<1>(call_gelyap_opt<float>, "X", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamOptSet);
    auto gelyap_float_opt_convlog              = makeFunctionExt<2>(call_gelyap_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamOptSet);

    auto gelyap_float_x0_opt                   = makeFunctionExt<1>(call_gelyap_x0_opt<float>, "X", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamFloatX0, ParamOptSet);
    auto gelyap_float_x0_opt_convlog           = makeFunctionExt<2>(call_gelyap_x0_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamFloatX0, ParamOptSet);

    auto gelyap_float_op_opt                   = makeFunctionExt<1>(call_gelyap_op_opt<float>, "X", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamOpString, ParamOptSet);
    auto gelyap_float_op_opt_convlog              = makeFunctionExt<2>(call_gelyap_op_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamOpString, ParamOptSet);

    auto gelyap_float_x0_op_opt                   = makeFunctionExt<1>(call_gelyap_x0_op_opt<float>, "X", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamFloatX0, ParamOpString, ParamOptSet);
    auto gelyap_float_x0_op_opt_convlog           = makeFunctionExt<2>(call_gelyap_x0_op_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatS, ParamFloatQ, ParamFloatY, ParamFloatX0, ParamOpString, ParamOptSet);

	auto parser = makeArgParser("mepack_lyap_refine", gelyap_double,  gelyap_double_convlog,
                                                      gelyap_double_x0, gelyap_double_x0_convlog,
                                                      gelyap_double_op,  gelyap_double_op_convlog,
                                                      gelyap_double_x0_op, gelyap_double_x0_op_convlog,
                                                      gelyap_double_opt,  gelyap_double_opt_convlog,
                                                      gelyap_double_x0_opt, gelyap_double_x0_opt_convlog,
                                                      gelyap_double_op_opt,  gelyap_double_op_opt_convlog,
                                                      gelyap_double_x0_op_opt, gelyap_double_x0_op_opt_convlog,
                                                      gelyap_float,  gelyap_float_convlog,
                                                      gelyap_float_x0, gelyap_float_x0_convlog,
                                                      gelyap_float_op,  gelyap_float_op_convlog,
                                                      gelyap_float_x0_op, gelyap_float_x0_op_convlog,
                                                      gelyap_float_opt,  gelyap_float_opt_convlog,
                                                      gelyap_float_x0_opt, gelyap_float_x0_opt_convlog,
                                                      gelyap_float_op_opt,  gelyap_float_op_opt_convlog,
                                                      gelyap_float_x0_op_opt, gelyap_float_x0_op_opt_convlog
                                                      );

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


