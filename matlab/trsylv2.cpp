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
struct tr2_solver {
    constexpr static const char *fn = NULL;
    constexpr static const char *fn_openmp = NULL;

    virtual void solve(const char *opA, const char * opB, T sign, int M, int N, T * Aptr, int lda, T* Bptr, int ldb, T *Yptr, int ldy, T * scale, T * buffer, int *info) = 0;
    virtual void solve_openmp(const char *opA, const char *opB, T sign, int M, int N, T * Aptr, int lda, T* Bptr, int ldb, T *Yptr, int ldy, T * scale, T * buffer, int *info) = 0;
    virtual void blocksize_mb_set(int mb) = 0;
    virtual void blocksize_nb_set(int nb) = 0;
    virtual void isolver_set(int) = 0;
};

template<>
struct tr2_solver<double> {
    constexpr static const char *fn = "DLA_TRSYLV2_L3";
    constexpr static const char *fn_openmp = "DLA_TRSYLV2_DAG";

    static void solve(const char *opA, const char *opB, double sign, int M, int N, double * Aptr, int lda, double *Bptr, int ldb, double *Yptr, int ldy, double * scale, double * buffer, int *info) {
        mepack_double_trsylv2_level3(opA, opB, sign, M, N, Aptr, lda, Bptr, ldb, Yptr, ldy, scale, buffer, info);
    };

    static void solve_openmp(const char *opA, const char *opB, double sign, int M, int N, double * Aptr, int lda, double *Bptr, int ldb, double *Yptr, int ldy, double * scale, double * buffer, int *info) {
        mepack_double_trsylv2_dag(opA, opB, sign, M, N, Aptr, lda, Bptr, ldb, Yptr, ldy, scale, buffer, info);
    };

    static void blocksize_mb_set(int mb) {
        mepack_double_trsylv2_blocksize_mb_set(mb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_double_trsylv2_blocksize_nb_set(nb);
    };

    static void isolver_set(int isolver) {
        mepack_trsylv2_isolver_set(isolver);
    };
};

template<>
struct tr2_solver<float> {
    constexpr static const char *fn = "SLA_TRSYLV2_L3";
    constexpr static const char *fn_openmp = "SLA_TRSYLV2_DAG";

    static void solve(const char *opA, const char *opB, float sign, int M, int N, float * Aptr, int lda, float *Bptr, int ldb, float *Yptr, int ldy, float * scale, float * buffer, int *info) {
        mepack_single_trsylv2_level3(opA, opB, sign, M, N, Aptr, lda, Bptr, ldb, Yptr, ldy, scale, buffer, info);
    };

    static void solve_openmp(const char *opA, const char *opB, float sign, int M, int N, float * Aptr, int lda, float *Bptr, int ldb, float *Yptr, int ldy, float * scale, float * buffer, int *info) {
        mepack_single_trsylv2_dag(opA, opB, sign, M, N, Aptr, lda, Bptr, ldb, Yptr, ldy, scale, buffer, info);
    };

    static void blocksize_mb_set(int mb) {
        mepack_single_trsylv2_blocksize_mb_set(mb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_single_trsylv2_blocksize_nb_set(nb);
    };

    static void isolver_set(int isolver) {
        mepack_trsylv2_isolver_set(isolver);
    };
};




template <typename T>
static std::tuple<ArrayType<T>>
call_solver_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y,
                      std::string opA = std::string("N"), std::string opB = std::string("N"),
                      MepackOptions optset = MepackOptions())
{
    if ( A.rows != Y.rows) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments one and three need to have the same number of rows size(A) = [ %d %d ] size (Y) = [ %d %d ]", (int) A.rows, (int) A.columns, (int) Y.rows, (int) Y.columns);
    }
    if ( B.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments two and three need to have the same number of columns size(B) = [ %d %d ] size (Y) = [ %d %d ]", (int) B.rows, (int) B.columns, (int) Y.rows, (int) Y.columns);
    }

    if ( optset.getMB() >= 0 ) {
        tr2_solver<T>::blocksize_mb_set(optset.getMB());
    }
    if ( optset.getNB() >= 0 ) {
        tr2_solver<T>::blocksize_nb_set(optset.getNB());
    }

    if ( optset.getISolver() >= 0 ) {
        tr2_solver<T>::isolver_set(optset.getISolver());
    }
    ssize_t auxmem;
    if ( optset.getOpenMP()) {
        auxmem = mepack_memory(tr2_solver<T>::fn_openmp, A.rows, B.rows);
    } else {
        auxmem = mepack_memory(tr2_solver<T>::fn, A.rows, B.rows);
    }

    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    T scale = 1;
    T sign = 1;
    int info = 0;

    if ( optset.getSign() < 0 ) {
        sign = -1;
    }

    if ( optset.getOpenMP()) {
        tr2_solver<T>::solve_openmp(opA.c_str(), opB.c_str(), sign, A.rows, B.rows, A.ptr(), A.ld, B.ptr(), B.ld, Y.ptr(), Y.ld, &scale, auxbuf, &info);
    } else {
        tr2_solver<T>::solve(opA.c_str(), opB.c_str(), sign, A.rows, B.rows, A.ptr(), A.ld, B.ptr(), B.ld, Y.ptr(), Y.ld, &scale, auxbuf, &info);
    }
    if ( info != 0) {
        mexoct_error("MEPACK:TRSYLV2_LEVEL3", "TRSYLV2_LEVEL3 failed with info = %d", (int) info);
    }

    return std::make_tuple(Y);


}

template<typename T>
static std::tuple<ArrayType<T>>
call_solver(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y)
{
    return call_solver_op_optset(A, B, Y);
}

template<typename T>
static std::tuple<ArrayType<T>>
call_solver_op(ArrayType<T>  A,ArrayType<T> B,  ArrayType<T> Y, std::string opA, std::string opB)
{
    return call_solver_op_optset(A, B, Y, opA, opB);
}

template<typename T>
static std::tuple<ArrayType<T>>
call_solver_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, MepackOptions optset)
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
    return call_solver_op_optset(A, B, Y, "N", "N",optset);
}



MEXOCT_ENTRY(mepack_trsylv2,
R"===(
X = mepack_trsylv2(A, B, Y)
X = mepack_trsylv2(A, B, Y, OpA, OpB)
X = mepack_trsylv2(A, B, Y, OptSet)
X = mepack_trsylv2(A, B, Y, OpA, OpB, OptSet)
    A      -- double precision (quasi) upper triangular coefficient matrix
    B      -- double precision (quasi) upper triangular coefficient matrix
    Y      -- double precision right hand side
    OpA    -- Transpose operation on the left coefficient matrix
    OpB    -- Transpose operation on the right coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

Xs = mepack_trsylv2(As, Bs, Ys)
Xs = mepack_trsylv2(As, Bs, Ys, Op)
Xs = mepack_trsylv2(As, Bs, Ys, OptSet)
Xs = mepack_trsylv2(As, Bs, Ys, Op, OptSet)
    As     -- single precision (quasi) upper triangular coefficient matrix
    BS     -- single precision (quasi) upper triangular coefficient matrix
    Ys     -- single precision right hand side
    OpA    -- Transpose operation on the left coefficient matrix
    OpB    -- Transpose operation on the right coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

The mepack_trsylv function solves the Sylvester equation

   opA(A) * X * opB(B) + X = Y,

or

   opA(A) * X * op(B) - X  = Y,

where A and B are (quasi) upper triangular matrices, both coming from
a generalized Schur decomposition. The operation op
either returns the matrix or if Op == 'T' its transpose. In this way, the
transposed equation can be solved without transposing the coefficient
matrix before.

The function uses the level-3 solver D/SLA_TRSYLV2_* of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_trsylv2:

   - mb      -- row block size of the level-3 algorithm (mb >=2)
   - nb      -- column block size of the level-3 algorithm (nb >=2)
   - isolver -- internal level-2 solver in the algorithm
       = 0   -- Default/Automatic selection
       = 1   -- Local Copies, with alignments (if supported by the compiler)
       = 2   -- Local Copies, without alignments
       = 3   -- Reordering and Fortran Vectorization
       = 4   -- Classic level-2 solver
       = 5   -- Naive level-2 solver
       = 6   -- Recursive Blocking
   - sign    -- swap the sign between both terms
       = 1   -- solve the equation with '+'
       = -1  -- solve the equation with '-'
   - openmp  -- OpenMP selection
       = 1   -- OpenMP - DAG solver enabled and OpenMP accelerated operations
       = 0   -- Standard level-3 solver without OpenMP acceleration


This file file is part of MEPACK, the Matrix Equation Package
by Martin Koehler.
)===")
{

    MEXOCT_INIT();
    mepack_init();

    auto ParamOpAString = Parameter<std::string>("OpA", "Transpose operation on the left coefficient matrix");
    auto ParamOpBString = Parameter<std::string>("OpB", "Transpose operation on the right coefficient matrix");
    auto ParamOptSet   = Parameter<MepackOptions>("optset", "Options for MEPACK");

    auto ParamDoubleA = Parameter<ArrayType<double>, Traits::SquareMatrix>("A", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleB = Parameter<ArrayType<double>, Traits::SquareMatrix>("B", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleY = Parameter<ArrayType<double>>("Y", "double precision right hand side");

    auto ParamFloatA = Parameter<ArrayType<float>, Traits::SquareMatrix>("As", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatB = Parameter<ArrayType<float>, Traits::SquareMatrix>("Bs", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatY = Parameter<ArrayType<float>>("Ys", "single precision right hand side");


    auto trsylv_double            = makeFunctionExt<1>(call_solver<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY);
    auto trsylv_double_optset     = makeFunctionExt<1>(call_solver_optset<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOptSet);
    auto trsylv_double_op         = makeFunctionExt<1>(call_solver_op<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpAString, ParamOpBString);
    auto trsylv_double_op_optset  = makeFunctionExt<1>(call_solver_op_optset<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpAString,ParamOpBString,ParamOptSet);



    auto trsylv_float            = makeFunctionExt<1>(call_solver<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY);
    auto trsylv_float_optset     = makeFunctionExt<1>(call_solver_optset<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOptSet);
    auto trsylv_float_op         = makeFunctionExt<1>(call_solver_op<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpAString, ParamOpBString);
    auto trsylv_float_op_optset  = makeFunctionExt<1>(call_solver_op_optset<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpAString, ParamOpBString, ParamOptSet);


	auto parser = makeArgParser("mepack_trsylv2", trsylv_double,  trsylv_double_optset, trsylv_double_op, trsylv_double_op_optset,
                                                 trsylv_float, trsylv_float_op, trsylv_float_optset, trsylv_float_op_optset);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


