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
struct stein_solver {
    constexpr static const char *fn = NULL;
    virtual void solve(const char * FACT, const char *TRANSA, int M, T * A,  int LDA, T * Q, int LDQ, T *X, int LDX, T * SCALE, T *WORK, size_t LDWORK, int *INFO) = 0;
    virtual void blocksize_set(int nb) = 0;
    virtual void blocksize_big_set(int nb) = 0;
    virtual void isolver_set(int is) = 0;
    virtual void fsolver_set(int fs) = 0;
};

template<>
struct stein_solver<double> {
    constexpr static const char *fn = "DLA_GESTEIN";

    static void solve(const char * fact, const char *op, int M, double * Aptr, int lda, double *Qptr, int ldq, double *Yptr, int ldy, double * scale, double * buffer, size_t ldwork, int *info) {
        mepack_double_gestein(fact, op, M, Aptr, lda, Qptr, ldq, Yptr, ldy, scale, buffer, ldwork, info);
    };
    static void blocksize_set(int nb) {
        mepack_double_trstein_blocksize_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_trstein_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_trstein_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_double_trstein_blocksize_2stage_set(nb);
    };

};

template<>
struct stein_solver<float> {
    constexpr static const char *fn = "SLA_GESTEIN";

    static void solve(const char * fact, const char *op, int M, float * Aptr, int lda, float *Qptr, int ldq, float *Yptr, int ldy, float * scale, float * buffer, size_t ldwork, int *info) {
        mepack_single_gestein(fact, op, M, Aptr, lda, Qptr, ldq, Yptr, ldy, scale, buffer, ldwork, info);
    };
    static void blocksize_set(int nb) {
        mepack_single_trstein_blocksize_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_trstein_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_trstein_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_single_trstein_blocksize_2stage_set(nb);
    };
};


template <typename T>
static std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>>
SteinDriver( ArrayType<T> A, ArrayType<T> Y, ArrayType<T> Q, std::string opA = std::string("N"), MepackOptions optset = MepackOptions(), bool factorize = false)
{
    if (A.rows != Y.rows || A.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix and the right hand side have to be of the same dimension size(A) = [ %d %d ] size (Y) = [ %d %d ]", (int) A.rows, (int) A.columns,
                (int) Y.rows, (int) Y.columns);
    }
    if (Q.rows != Y.rows || Q.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The unitary transformation matrix and the right hand side have to be of the same dimension size(Q) = [ %d %d ] size (Y) = [ %d %d ]", (int) Q.rows, (int) Q.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if ( optset.getMB() >= 0 ) {
        stein_solver<T>::blocksize_set(optset.getMB());
    }
    if ( optset.getISolver() >= 0 ) {
        stein_solver<T>::isolver_set(optset.getISolver());
    }
    if ( optset.getFSolver() >= 0 ) {
        stein_solver<T>::fsolver_set(optset.getFSolver());
    }

    ssize_t auxmem;

    if ( factorize ) {
        auxmem = mepack_memory_frontend(stein_solver<T>::fn, (const char *){"N"}, (const char *){"N"}, A.rows, A.columns);
    } else {
        auxmem = mepack_memory_frontend(stein_solver<T>::fn, (const char *){"F"}, (const char *){"F"}, A.rows, A.columns);
    }
    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    T scale = 1;
    int info = 0;

    if ( factorize ) {
        stein_solver<T>::solve("N", opA.c_str(), A.rows, A.ptr(), A.ld, Q.ptr(), Q.ld, Y.ptr(), Y.ld, &scale, auxbuf, auxmem, &info);
    } else {
        stein_solver<T>::solve("F", opA.c_str(), A.rows, A.ptr(), A.ld, Q.ptr(), Q.ld, Y.ptr(), Y.ld, &scale, auxbuf, auxmem, &info);
    }

    if ( info != 0) {
        mexoct_error("MEPACK:GESTEIN", "GESTEIN failed with info = %d", (int) info);
    }
    return std::make_tuple(Y, A, Q);
}


template <typename T>
std::tuple<ArrayType<T>>
call_gestein_op_optset(ArrayType<T>  A, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    ArrayType<T> Q(A.rows, A.rows);
    auto ret = SteinDriver(A, Y, Q, opA, optset, true);
    return std::make_tuple(std::get<0>(ret));
}

template<typename T>
std::tuple<ArrayType<T>>
call_gestein(ArrayType<T>  A, ArrayType<T> Y)
{
    return call_gestein_op_optset(A, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gestein_op(ArrayType<T>  A, ArrayType<T> Y, std::string opA)
{
    return call_gestein_op_optset(A, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gestein_optset(ArrayType<T>  A, ArrayType<T> Y, MepackOptions optset)
{
   return call_gestein_op_optset(A, Y, "N", optset);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_3ret_gestein_op_optset(ArrayType<T>  A, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    ArrayType<T> Q(A.rows, A.rows);
    return SteinDriver(A, Y, Q, opA, optset, true);
}



template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>,ArrayType<T>>
call_3ret_gestein(ArrayType<T>  A, ArrayType<T> Y)
{
    return call_3ret_gestein_op_optset(A, Y);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_3ret_gestein_op(ArrayType<T>  A, ArrayType<T> Y, std::string opA)
{
    return call_3ret_gestein_op_optset(A, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>, ArrayType<T>>
call_3ret_gestein_optset(ArrayType<T>  A, ArrayType<T> Y, MepackOptions optset)
{
   return call_3ret_gestein_op_optset(A, Y, "N", optset);
}

template <typename T>
std::tuple<ArrayType<T>>
call_gestein_prefact_op_optset(ArrayType<T>  S, ArrayType<T> Q, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    bool chk = mepack_mexoct_check_qtriangular(S);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix is not (quasi) upper triangular.");
    }
    auto ret = SteinDriver(S, Y, Q, opA, optset, false);
    return std::make_tuple(std::get<0>(ret));
}

template<typename T>
std::tuple<ArrayType<T>>
call_gestein_prefact(ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y)
{
    return call_gestein_prefact_op_optset(S, Q, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gestein_prefact_op(ArrayType<T>  S, ArrayType<T> Q, ArrayType<T> Y, std::string opA)
{
    return call_gestein_prefact_op_optset(S, Q, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gestein_prefact_optset(ArrayType<T>  S, ArrayType<T> Q, ArrayType<T> Y, MepackOptions optset)
{
   return call_gestein_prefact_op_optset(S, Q, Y, "N", optset);
}



MEXOCT_ENTRY(mepack_stein,
R"===(
X = mepack_stein(A, Y)
X = mepack_stein(A, Y, Op)
X = mepack_stein(A, Y, OptSet)
X = mepack_stein(A, Y, Op, OptSet)
 Input:
    A      -- coefficient matrix
    Y      -- right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.
 Return:
    X      -- Solution of the Stein equation


[X, S, Q] = mepack_stein(A, Y)
[X, S, Q] = mepack_stein(A, Y, Op)
[X, S, Q] = mepack_stein(A, Y, OptSet)
[X, S, Q] = mepack_stein(A, Y, Op, OptSet)
  Input:
    A      -- coefficient matrix
    Y      -- right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

  Return:
    X      -- Solution of the Stein equation
    S, Q   -- Schur decomposition of A, A = Q*S*Q'

X = mepack_stein(S, Q, Y)
X = mepack_stein(S, Q, Y, Op)
X = mepack_stein(S, Q, Y, OptSet)
X = mepack_stein(S, Q, Y, Op, OptSet)
  Input:
    S      -- (quasi) upper triangular coefficient matrix, Schur form of A
    Q      -- unitary transformation matrix, A = Q*S*Q'
    Y      -- right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

  Return:
    X      -- Solution of the Stein equation

The mepack_stein function solves the Stein equation

   op(A) * X * op(A)' - X = Y,

where A is a general matrix. The operation op either returns the matrix or
if Op == 'T' its transpose. In this way, the transposed equation
can be solved without transposing the coefficient matrix before. The function
can also be called with A decomposed in to A = Q*S*Q', where S is the
Schur form of A. If the function is called with three return values,
the Schur form of A and its unitary transformation matrix is returned.

The matrices A and Y must be of the same size and Y must be symmetric.

The function uses the level-3 solver D/SLA_GESTEIN of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_stein:

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

    auto ParamDoubleA = Parameter<ArrayType<double>,Traits::SquareMatrix>("A", "double precision coefficient matrix");
    auto ParamDoubleS = Parameter<ArrayType<double>,Traits::SquareMatrix>("S", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleQ = Parameter<ArrayType<double>,Traits::SquareMatrix>("Q", "double precision unitary matrix corresponding to S");
    auto ParamDoubleY = Parameter<ArrayType<double>,Traits::SquareMatrix>("Y", "double precision right hand side");

    auto ParamFloatA = Parameter<ArrayType<float>,Traits::SquareMatrix>("As", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatY = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ys", "single precision right hand side");
    auto ParamFloatS = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ss", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatQ = Parameter<ArrayType<float>,Traits::SquareMatrix>("Qs", "single precision unitary matrix corresponding to Ss");


    auto gestein_double            = makeFunctionExt<1>(call_gestein<double>, "X", ParamDoubleA, ParamDoubleY);
    auto gestein_double_optset     = makeFunctionExt<1>(call_gestein_optset<double>, "X", ParamDoubleA, ParamDoubleY, ParamOptSet);
    auto gestein_double_op         = makeFunctionExt<1>(call_gestein_op<double>, "X", ParamDoubleA, ParamDoubleY, ParamOpString);
    auto gestein_double_op_optset  = makeFunctionExt<1>(call_gestein_op_optset<double>, "X", ParamDoubleA, ParamDoubleY, ParamOpString, ParamOptSet);
    auto gestein_3ret_double            = makeFunctionExt<3>(call_3ret_gestein<double>, "[X, S, Q]", ParamDoubleA, ParamDoubleY);
    auto gestein_3ret_double_optset     = makeFunctionExt<3>(call_3ret_gestein_optset<double>, "[X, S, Q]", ParamDoubleA, ParamDoubleY, ParamOptSet);
    auto gestein_3ret_double_op         = makeFunctionExt<3>(call_3ret_gestein_op<double>, "[X, S, Q]", ParamDoubleA, ParamDoubleY, ParamOpString);
    auto gestein_3ret_double_op_optset  = makeFunctionExt<3>(call_3ret_gestein_op_optset<double>, "[X, S, Q]", ParamDoubleA, ParamDoubleY, ParamOpString, ParamOptSet);
    auto gestein_prefact_double            = makeFunctionExt<1>(call_gestein_prefact<double>, "X",           ParamDoubleS, ParamDoubleQ, ParamDoubleY);
    auto gestein_prefact_double_optset     = makeFunctionExt<1>(call_gestein_prefact_optset<double>, "X",    ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOptSet);
    auto gestein_prefact_double_op         = makeFunctionExt<1>(call_gestein_prefact_op<double>, "X",        ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOpString);
    auto gestein_prefact_double_op_optset  = makeFunctionExt<1>(call_gestein_prefact_op_optset<double>, "X", ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOpString, ParamOptSet);




    auto gestein_float            = makeFunctionExt<1>(call_gestein<float>, "Xs", ParamFloatA, ParamFloatY);
    auto gestein_float_optset     = makeFunctionExt<1>(call_gestein_optset<float>, "Xs", ParamFloatA, ParamFloatY, ParamOptSet);
    auto gestein_float_op         = makeFunctionExt<1>(call_gestein_op<float>, "Xs", ParamFloatA, ParamFloatY, ParamOpString);
    auto gestein_float_op_optset  = makeFunctionExt<1>(call_gestein_op_optset<float>, "Xs", ParamFloatA, ParamFloatY, ParamOpString, ParamOptSet);
    auto gestein_3ret_float            = makeFunctionExt<3>(call_3ret_gestein<float>, "[Xs, Ss, Qs]", ParamFloatA, ParamFloatY);
    auto gestein_3ret_float_optset     = makeFunctionExt<3>(call_3ret_gestein_optset<float>, "[Xs, Ss, Qs]", ParamFloatA, ParamFloatY, ParamOptSet);
    auto gestein_3ret_float_op         = makeFunctionExt<3>(call_3ret_gestein_op<float>, "[Xs, Ss, Qs]", ParamFloatA, ParamFloatY, ParamOpString);
    auto gestein_3ret_float_op_optset  = makeFunctionExt<3>(call_3ret_gestein_op_optset<float>, "[Xs, Ss, Qs]", ParamFloatA, ParamFloatY, ParamOpString, ParamOptSet);
    auto gestein_prefact_float            = makeFunctionExt<1>(call_gestein_prefact<float>, "Xs",           ParamFloatS, ParamFloatQ, ParamFloatY);
    auto gestein_prefact_float_optset     = makeFunctionExt<1>(call_gestein_prefact_optset<float>, "Xs",    ParamFloatS, ParamFloatQ, ParamFloatY, ParamOptSet);
    auto gestein_prefact_float_op         = makeFunctionExt<1>(call_gestein_prefact_op<float>, "Xs",        ParamFloatS, ParamFloatQ, ParamFloatY, ParamOpString);
    auto gestein_prefact_float_op_optset  = makeFunctionExt<1>(call_gestein_prefact_op_optset<float>, "Xs", ParamFloatS, ParamFloatQ, ParamFloatY, ParamOpString, ParamOptSet);





	auto parser = makeArgParser("mepack_stein", gestein_double,  gestein_double_optset, gestein_double_op, gestein_double_op_optset,
                                                 gestein_3ret_double,  gestein_3ret_double_optset, gestein_3ret_double_op, gestein_3ret_double_op_optset,
                                                 gestein_prefact_double, gestein_prefact_double_optset, gestein_prefact_double_op, gestein_prefact_double_op_optset,
                                                 gestein_3ret_float,  gestein_3ret_float_optset, gestein_3ret_float_op, gestein_3ret_float_op_optset,
                                                 gestein_prefact_float, gestein_prefact_float_optset, gestein_prefact_float_op, gestein_prefact_float_op_optset,
                                                 gestein_float, gestein_float_op, gestein_float_optset, gestein_float_op_optset);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


