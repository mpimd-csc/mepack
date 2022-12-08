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
 * Copyright (C) Martin Koehler, 2017-2022
 */

#define MEXOCT_MATLAB_DQ_STRING

#include "mepack.h"
#include "mexoct/interface.hpp"
#include "mepackoptions.hpp"

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
struct tglyap_solver {
    constexpr static const char *fn = NULL;
    constexpr static const char *fn_openmp = NULL;

    virtual void solve(const char *op, int M, T * Aptr, int lda, T* Bptr, int ldb, T *Yptr, int ldy, T * scale, T * buffer, int *info) = 0;
    virtual void solve_openmp(const char *op, int M, T * Aptr, int lda, T* Bptr, int ldb, T *Yptr, int ldy, T * scale, T * buffer, int *info) = 0;
    virtual void blocksize_set(int nb) = 0;
    virtual void isolver_set(int) = 0;
};

template<>
struct tglyap_solver<double> {
    constexpr static const char *fn = "DLA_TGLYAP_L3";
    constexpr static const char *fn_openmp = "DLA_TGLYAP_DAG";

    static void solve(const char *op, int M, double * Aptr, int lda, double *Bptr, int ldb, double *Yptr, int ldy, double * scale, double * buffer, int *info) {
        mepack_double_tglyap_level3(op, M, Aptr, lda, Bptr, ldb, Yptr, ldy, scale, buffer, info);
    };
    static void solve_openmp(const char *op, int M, double * Aptr, int lda, double *Bptr, int ldb, double *Yptr, int ldy, double * scale, double * buffer, int *info) {
        mepack_double_tglyap_dag(op, M, Aptr, lda, Bptr, ldb, Yptr, ldy, scale, buffer, info);
    };
    static void blocksize_set(int nb) {
        mepack_double_tglyap_blocksize_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_tglyap_isolver_set(isolver);
    };
};

template<>
struct tglyap_solver<float> {
    constexpr static const char *fn = "SLA_TGLYAP_L3";
    constexpr static const char *fn_openmp = "SLA_TGLYAP_DAG";

    static void solve(const char *op, int M, float * Aptr, int lda, float *Bptr, int ldb, float *Yptr, int ldy, float * scale, float * buffer, int *info) {
        mepack_single_tglyap_level3(op, M, Aptr, lda, Bptr, ldb, Yptr, ldy, scale, buffer, info);
    };
    static void solve_openmp(const char *op, int M, float * Aptr, int lda, float *Bptr, int ldb, float *Yptr, int ldy, float * scale, float * buffer, int *info) {
        mepack_single_tglyap_dag(op, M, Aptr, lda, Bptr, ldb, Yptr, ldy, scale, buffer, info);
    };
    static void blocksize_set(int nb) {
        mepack_single_tglyap_blocksize_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_tglyap_isolver_set(isolver);
    };
};



template <typename T>
std::tuple<ArrayType<T>>
call_tglyap_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    if ( A.rows != Y.rows && A.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments one and three need to be of the same dimension size(A) = [ %d %d ] size (Y) = [ %d %d ]", (int) A.rows, (int) A.columns, (int) Y.rows, (int) Y.columns);
    }
    if ( B.rows != Y.rows && B.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments two and three need to be of the same dimension size(B) = [ %d %d ] size (Y) = [ %d %d ]", (int) B.rows, (int) B.columns, (int) Y.rows, (int) Y.columns);
    }

    if ( optset.getMB() >= 0 ) {
        tglyap_solver<T>::blocksize_set(optset.getMB());
    }
    if ( optset.getISolver() >= 0 ) {
        tglyap_solver<T>::isolver_set(optset.getISolver());
    }
    ssize_t auxmem;
    if ( optset.getOpenMP()) {
        auxmem = mepack_memory(tglyap_solver<T>::fn_openmp, A.rows, A.columns);
    } else {
        auxmem = mepack_memory(tglyap_solver<T>::fn, A.rows, A.columns);
    }

    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    T scale = 1;
    int info = 0;

    if ( optset.getOpenMP()) {
        tglyap_solver<T>::solve_openmp(opA.c_str(), A.rows, A.ptr(), A.ld, B.ptr(), B.ld, Y.ptr(), Y.ld, &scale, auxbuf, &info);
    } else {
        tglyap_solver<T>::solve(opA.c_str(), A.rows, A.ptr(), A.ld, B.ptr(), B.ld, Y.ptr(), Y.ld, &scale, auxbuf, &info);
    }
    if ( info != 0) {
        mexoct_error("MEPACK:TGLYAP_LEVEL3", "TGLYAP_LEVEL3 failed with info = %d", (int) info);
    }

    return std::make_tuple(Y);


}

template<typename T>
std::tuple<ArrayType<T>>
call_tglyap(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y)
{
    return call_tglyap_op_optset(A, B, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_tglyap_op(ArrayType<T>  A,ArrayType<T> B,  ArrayType<T> Y, std::string opA)
{
    return call_tglyap_op_optset(A, B, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>>
call_tglyap_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, MepackOptions optset)
{
// #define DEBUG_ALL
#ifdef DEBUG_ALL
    std::cout << "Optset Function \n";
    std::cout << " - MB " << optset.getMB() << "\n";
    std::cout << " - NB " << optset.getNB() << "\n";
    std::cout << " - BIGNB " << optset.getBigNB() << "\n";
    std::cout << " - ISOLVER " << optset.getISolver() << "\n";
    std::cout << " - FSOLVER " << optset.getFSolver() << "\n";
#endif
    return call_tglyap_op_optset(A, B, Y, "N", optset);
}



MEXOCT_ENTRY(mepack_tglyap,
R"===(
X = mepack_tglyap(A, B, Y)
X = mepack_tglyap(A, B, Y, Op)
X = mepack_tglyap(A, B, Y, OptSet)
X = mepack_tglyap(A, B, Y, Op, OptSet)
    A      -- double precision (quasi) upper triangular coefficient matrix
    B      -- double precision upper triangular coefficient matrix
    Y      -- double precision right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

Xs = mepack_tglyap(As, Bs, Ys)
Xs = mepack_tglyap(As, Bs, Ys, Op)
Xs = mepack_tglyap(As, Bs, Ys, OptSet)
Xs = mepack_tglyap(As, Bs, Ys, Op, OptSet)
    As     -- single precision (quasi) upper triangular coefficient matrix
    BS     -- single precision upper triangular coefficient matrix
    Ys     -- single precision right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

The mepack_tglyap function solves the Lyapunov equation

   op(A) * X * op(B)' + op(B) * X * op(A)' = Y,

where A is a (quasi) upper triangular matrix and B is a upper triangular
matrix coming from a generalized Schur decomposition. The operation op
either returns the matrix or if Op == 'T' its transpose. In this way, the
transposed equation can be solved without transposing the coefficient
matrix before.

The matrices A, B, and Y must be of the same size and Y must be symmetric.

The function uses the level-3 solver D/SLA_TGLYAP_LEVEL3 of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_tglyap:

   - mb      -- block size of the level-3 algorithm (mb >=2)
   - isolver -- internal level-2 solver in the algorithm
       = 0   -- Default/Automatic selection
       = 1   -- Local Copies, with alignments (if supported by the compiler)
       = 2   -- Local Copies, without alignments
       = 3   -- Reordering and Fortran Vectorization
       = 4   -- Classic level-2 solver
       = 5   -- Naive level-2 solver
       = 6   -- Recursive Blocking
   - openmp  -- OpenMP selection
       = 1   -- OpenMP - DAG solver enabled and OpenMP accelerated operations
       = 0   -- Standard level-3 solver without OpenMP acceleration



This file file is part of MEPACK, the Matrix Equation Package
by Martin Koehler.
)===")
{

    MEXOCT_INIT();
    mepack_init();

    auto ParamOpString = Parameter<std::string>("Op", "Transpose operation on the coefficient matrix");
    auto ParamOptSet   = Parameter<MepackOptions>("optset", "Options for MEPACK");

    auto ParamDoubleA = Parameter<ArrayType<double>,Traits::SquareMatrix>("A", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleB = Parameter<ArrayType<double>,Traits::SquareMatrix>("B", "double precision upper triangular coefficient matrix");
    auto ParamDoubleY = Parameter<ArrayType<double>,Traits::SquareMatrix>("Y", "double precision right hand side");

    auto ParamFloatA = Parameter<ArrayType<float>,Traits::SquareMatrix>("As", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatB = Parameter<ArrayType<float>,Traits::SquareMatrix>("Bs", "single precision upper triangular coefficient matrix");
    auto ParamFloatY = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ys", "single precision right hand side");


    auto tglyap_double            = makeFunctionExt<1>(call_tglyap<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY);
    auto tglyap_double_optset     = makeFunctionExt<1>(call_tglyap_optset<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOptSet);
    auto tglyap_double_op         = makeFunctionExt<1>(call_tglyap_op<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpString);
    auto tglyap_double_op_optset  = makeFunctionExt<1>(call_tglyap_op_optset<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpString, ParamOptSet);



    auto tglyap_float            = makeFunctionExt<1>(call_tglyap<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY);
    auto tglyap_float_optset     = makeFunctionExt<1>(call_tglyap_optset<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOptSet);
    auto tglyap_float_op         = makeFunctionExt<1>(call_tglyap_op<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpString);
    auto tglyap_float_op_optset  = makeFunctionExt<1>(call_tglyap_op_optset<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpString, ParamOptSet);




	auto parser = makeArgParser("mepack_tglyap", tglyap_double,  tglyap_double_optset, tglyap_double_op, tglyap_double_op_optset,
                                                 tglyap_float, tglyap_float_op, tglyap_float_optset, tglyap_float_op_optset);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


