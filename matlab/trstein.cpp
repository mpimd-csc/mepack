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
struct trstein_solver {
    constexpr static const char *fn = NULL;
    constexpr static const char *fn_openmp = NULL;
    virtual void solve(const char *op, int M, T * Aptr, int lda, T *Yptr, int ldy, T * scale, T * buffer, int *info) = 0;
    virtual void solve_openmp(const char *op, int M, T * Aptr, int lda, T *Yptr, int ldy, T * scale, T * buffer, int *info) = 0;
    virtual void blocksize_set(int nb) = 0;
    virtual void isolver_set(int) = 0;
};


template<>
struct trstein_solver<double> {
    constexpr static const char *fn = "DLA_TRSTEIN_L3";
    constexpr static const char *fn_openmp = "DLA_TRSTEIN_DAG";


    static void solve(const char *op, int M, double * Aptr, int lda, double *Yptr, int ldy, double * scale, double * buffer, int *info) {
        mepack_double_trstein_level3(op, M, Aptr, lda, Yptr, ldy, scale, buffer, info);
    };
    static void solve_openmp(const char *op, int M, double * Aptr, int lda, double *Yptr, int ldy, double * scale, double * buffer, int *info) {
        mepack_double_trstein_dag(op, M, Aptr, lda, Yptr, ldy, scale, buffer, info);
    };

    static void blocksize_set(int nb) {
        mepack_double_trstein_blocksize_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_trstein_isolver_set(isolver);
    };
};

template<>
struct trstein_solver<float> {
    constexpr static const char *fn = "SLA_TRSTEIN_L3";
    constexpr static const char *fn_openmp = "SLA_TRSTEIN_DAG";

    static void solve(const char *op, int M, float * Aptr, int lda, float *Yptr, int ldy, float * scale, float * buffer, int *info) {
        mepack_single_trstein_level3(op, M, Aptr, lda, Yptr, ldy, scale, buffer, info);
    };
    static void solve_openmp(const char *op, int M, float * Aptr, int lda, float *Yptr, int ldy, float * scale, float * buffer, int *info) {
        mepack_single_trstein_dag(op, M, Aptr, lda, Yptr, ldy, scale, buffer, info);
    };
    static void blocksize_set(int nb) {
        mepack_single_trstein_blocksize_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_trstein_isolver_set(isolver);
    };
};



template <typename T>
std::tuple<ArrayType<T>>
call_trstein_op_optset(ArrayType<T>  A, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    if ( A.rows != Y.rows && A.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments one and two need to be of the same dimension size(A) = [ %d %d ] size (Y) = [ %d %d ]", (int) A.rows, (int) A.columns, (int) Y.rows, (int) Y.columns);
    }

    if ( optset.getMB() >= 0 ) {
        trstein_solver<T>::blocksize_set(optset.getMB());
    }
    if ( optset.getISolver() >= 0 ) {
        trstein_solver<T>::isolver_set(optset.getISolver());
    }
    ssize_t auxmem;
    if ( optset.getOpenMP()) {
        auxmem = mepack_memory(trstein_solver<T>::fn_openmp, A.rows, A.columns);
    } else {
        auxmem = mepack_memory(trstein_solver<T>::fn, A.rows, A.columns);
    }
    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    T scale = 1;
    int info = 0;

    if (optset.getOpenMP()) {
        trstein_solver<T>::solve_openmp(opA.c_str(), A.rows, A.ptr(), A.ld, Y.ptr(), Y.ld, &scale, auxbuf, &info);
    } else {
        trstein_solver<T>::solve(opA.c_str(), A.rows, A.ptr(), A.ld, Y.ptr(), Y.ld, &scale, auxbuf, &info);
    }

    if ( info != 0) {
        mexoct_error("MEPACK:TRSTEIN_LEVEL3", "TRSTEIN_LEVEL3 failed with info = %d", (int) info);
    }

    return std::make_tuple(Y);


}

template<typename T>
std::tuple<ArrayType<T>>
call_trstein(ArrayType<T>  A, ArrayType<T> Y)
{
    return call_trstein_op_optset(A, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_trstein_op(ArrayType<T>  A, ArrayType<T> Y, std::string opA)
{
    return call_trstein_op_optset(A, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>>
call_trstein_optset(ArrayType<T>  A, ArrayType<T> Y, MepackOptions optset)
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
    return call_trstein_op_optset(A, Y, "N", optset);
}





MEXOCT_ENTRY(mepack_trstein,
R"===(
X = mepack_trstein(A, Y)
X = mepack_trstein(A, Y, Op)
X = mepack_trstein(A, Y, OptSet)
X = mepack_trstein(A, Y, Op, OptSet)
    A      -- double precision coefficient matrix
    Y      -- double precision right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

Xs = mepack_trstein(As, Ys)
Xs = mepack_trstein(As, Ys, Op)
Xs = mepack_trstein(As, Ys, OptSet)
Xs = mepack_trstein(As, Ys, Op, OptSet)
    As     -- single precision coefficient matrix
    Ys     -- single precision right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

The mepack_trstein function solves the Stein equation

   op(A) * X op(A)' - X = Y,

where A is a (quasi) upper triangular matrix coming from a Schur
decomposition. The operation op either returns the matrix or
if Op == 'T' its transpose. In this way, the transposed equation
can be solved without transposing the coefficient matrix before.

The matrices A and Y must be of the same size and Y must be symmetric.

The function uses the level-3 solver D/SLA_TRSTEIN_LEVEL3 of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_trstein:

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

    auto ParamDoubleA = Parameter<ArrayType<double>,Traits::SquareMatrix>("A", "double precision coefficient matrix");
    auto ParamDoubleY = Parameter<ArrayType<double>,Traits::SquareMatrix>("Y", "double precision right hand side");

    auto ParamFloatA = Parameter<ArrayType<float>,Traits::SquareMatrix>("As", "single precision coefficient matrix");
    auto ParamFloatY = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ys", "single precision right hand side");


    auto trstein_double            = makeFunctionExt<1>(call_trstein<double>, "X", ParamDoubleA, ParamDoubleY);
    auto trstein_double_optset     = makeFunctionExt<1>(call_trstein_optset<double>, "X", ParamDoubleA, ParamDoubleY, ParamOptSet);
    auto trstein_double_op         = makeFunctionExt<1>(call_trstein_op<double>, "X", ParamDoubleA, ParamDoubleY, ParamOpString);
    auto trstein_double_op_optset  = makeFunctionExt<1>(call_trstein_op_optset<double>, "X", ParamDoubleA, ParamDoubleY, ParamOpString, ParamOptSet);



    auto trstein_float            = makeFunctionExt<1>(call_trstein<float>, "Xs", ParamFloatA, ParamFloatY);
    auto trstein_float_optset     = makeFunctionExt<1>(call_trstein_optset<float>, "Xs", ParamFloatA, ParamFloatY, ParamOptSet);
    auto trstein_float_op         = makeFunctionExt<1>(call_trstein_op<float>, "Xs", ParamFloatA, ParamFloatY, ParamOpString);
    auto trstein_float_op_optset  = makeFunctionExt<1>(call_trstein_op_optset<float>, "Xs", ParamFloatA, ParamFloatY, ParamOpString, ParamOptSet);




	auto parser = makeArgParser("mepack_trstein", trstein_double,  trstein_double_optset, trstein_double_op, trstein_double_op_optset,
                                                 trstein_float, trstein_float_op, trstein_float_optset, trstein_float_op_optset);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


