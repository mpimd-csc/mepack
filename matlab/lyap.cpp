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
#include "hess_test.hpp"
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
struct lyap_solver {
    constexpr static const char *fn = NULL;
    virtual void solve(const char * FACT, const char *TRANSA, int M, T * A,  int LDA, T * Q, int LDQ, T *X, int LDX, T * SCALE, T *WORK, size_t LDWORK, int *INFO) = 0;
    virtual void blocksize_set(int nb) = 0;
    virtual void blocksize_big_set(int nb) = 0;
    virtual void isolver_set(int is) = 0;
    virtual void fsolver_set(int fs) = 0;
};

template<>
struct lyap_solver<double> {
    constexpr static const char *fn = "DLA_GELYAP";

    static void solve(const char * fact, const char *op, int M, double * Aptr, int lda, double *Qptr, int ldq, double *Yptr, int ldy, double * scale, double * buffer, size_t ldwork, int *info) {
        mepack_double_gelyap(fact, op, M, Aptr, lda, Qptr, ldq, Yptr, ldy, scale, buffer, ldwork, info);
    };
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
struct lyap_solver<float> {
    constexpr static const char *fn = "SLA_GELYAP";

    static void solve(const char * fact, const char *op, int M, float * Aptr, int lda, float *Qptr, int ldq, float *Yptr, int ldy, float * scale, float * buffer, size_t ldwork, int *info) {
        mepack_single_gelyap(fact, op, M, Aptr, lda, Qptr, ldq, Yptr, ldy, scale, buffer, ldwork, info);
    };
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
static std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>>
LyapDriver( ArrayType<T> A, ArrayType<T> Y, ArrayType<T> Q, std::string opA = std::string("N"), MepackOptions optset = MepackOptions(), bool factorize = false)
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
        lyap_solver<T>::blocksize_set(optset.getMB());
    }
    if ( optset.getISolver() >= 0 ) {
        lyap_solver<T>::isolver_set(optset.getISolver());
    }
    if ( optset.getFSolver() >= 0 ) {
        lyap_solver<T>::fsolver_set(optset.getFSolver());
    }

    ssize_t auxmem;

    std::string factA = std::string("N");

    if ( factorize ) {
        /* Check if the matrices are Hessenberg */
        bool a_hess = mepack_mexoct_check_hess(A);
        if (a_hess ) factA = std::string("H");
    } else {
        factA = std::string("F");
    }

    auxmem = mepack_memory_frontend(lyap_solver<T>::fn, factA.c_str(), factA.c_str(), A.rows, A.columns);

    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    T scale = 1;
    int info = 0;

    lyap_solver<T>::solve(factA.c_str() , opA.c_str(), A.rows, A.ptr(), A.ld, Q.ptr(), Q.ld, Y.ptr(), Y.ld, &scale, auxbuf, auxmem, &info);

    if ( info != 0) {
        mexoct_error("MEPACK:GELYAP", "GELYAP failed with info = %d", (int) info);
    }
    return std::make_tuple(Y, A, Q);
}


template <typename T>
std::tuple<ArrayType<T>>
call_gelyap_op_optset(ArrayType<T>  A, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    ArrayType<T> Q(A.rows, A.rows);
    auto ret = LyapDriver(A, Y, Q, opA, optset, true);
    return std::make_tuple(std::get<0>(ret));
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap(ArrayType<T>  A, ArrayType<T> Y)
{
    return call_gelyap_op_optset(A, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_op(ArrayType<T>  A, ArrayType<T> Y, std::string opA)
{
    return call_gelyap_op_optset(A, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_optset(ArrayType<T>  A, ArrayType<T> Y, MepackOptions optset)
{
   return call_gelyap_op_optset(A, Y, "N", optset);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_3ret_gelyap_op_optset(ArrayType<T>  A, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    ArrayType<T> Q(A.rows, A.rows);
    return LyapDriver(A, Y, Q, opA, optset, true);
}



template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>,ArrayType<T>>
call_3ret_gelyap(ArrayType<T>  A, ArrayType<T> Y)
{
    return call_3ret_gelyap_op_optset(A, Y);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_3ret_gelyap_op(ArrayType<T>  A, ArrayType<T> Y, std::string opA)
{
    return call_3ret_gelyap_op_optset(A, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>, ArrayType<T>>
call_3ret_gelyap_optset(ArrayType<T>  A, ArrayType<T> Y, MepackOptions optset)
{
   return call_3ret_gelyap_op_optset(A, Y, "N", optset);
}

template <typename T>
std::tuple<ArrayType<T>>
call_gelyap_prefact_op_optset(ArrayType<T>  S, ArrayType<T> Q, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    bool chk = mepack_mexoct_check_qtriangular(S);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix is not (quasi) upper triangular.");
    }
    auto ret = LyapDriver(S, Y, Q, opA, optset, false);
    return std::make_tuple(std::get<0>(ret));
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_prefact(ArrayType<T> S, ArrayType<T> Q, ArrayType<T> Y)
{
    return call_gelyap_prefact_op_optset(S, Q, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_prefact_op(ArrayType<T>  S, ArrayType<T> Q, ArrayType<T> Y, std::string opA)
{
    return call_gelyap_prefact_op_optset(S, Q, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gelyap_prefact_optset(ArrayType<T>  S, ArrayType<T> Q, ArrayType<T> Y, MepackOptions optset)
{
   return call_gelyap_prefact_op_optset(S, Q, Y, "N", optset);
}



MEXOCT_ENTRY(mepack_lyap,
R"===(
X = mepack_lyap(A, Y)
X = mepack_lyap(A, Y, Op)
X = mepack_lyap(A, Y, OptSet)
X = mepack_lyap(A, Y, Op, OptSet)
 Input:
    A      -- coefficient matrix
    Y      -- right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.
 Return:
    X      -- Solution of the Lyapunov equation


[X, S, Q] = mepack_lyap(A, Y)
[X, S, Q] = mepack_lyap(A, Y, Op)
[X, S, Q] = mepack_lyap(A, Y, OptSet)
[X, S, Q] = mepack_lyap(A, Y, Op, OptSet)
  Input:
    A      -- coefficient matrix
    Y      -- right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

  Return:
    X      -- Solution of the Lyapunov equation
    S, Q   -- Schur decomposition of A, A = Q*S*Q'

X = mepack_lyap(S, Q, Y)
X = mepack_lyap(S, Q, Y, Op)
X = mepack_lyap(S, Q, Y, OptSet)
X = mepack_lyap(S, Q, Y, Op, OptSet)
  Input:
    S      -- (quasi) upper triangular coefficient matrix, Schur form of A
    Q      -- unitary transformation matrix, A = Q*S*Q'
    Y      -- right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

  Return:
    X      -- Solution of the Lyapunov equation

The mepack_lyap function solves the Lyapunov equation

   op(A) * X + X * op(A)' = Y,

where A is a general matrix. The operation op either returns the matrix or
if Op == 'T' its transpose. In this way, the transposed equation
can be solved without transposing the coefficient matrix before. The function
can also be called with A decomposed in to A = Q*S*Q', where S is the
Schur form of A. If the function is called with three return values,
the Schur form of A and its unitary transformation matrix is returned.

The matrices A and Y must be of the same size and Y must be symmetric.

The function uses the level-3 solver D/SLA_GELYAP of MEPACK.

The involved Schur decomposition is aware of matrices in Hessenberg form.
In this case the initial Hessenberg reduction of the Schur decomposition is skipped.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_lyap:

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

    auto ParamFloatA = Parameter<ArrayType<float>,Traits::SquareMatrix>("As", "single precision coefficient matrix");
    auto ParamFloatY = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ys", "single precision right hand side");
    auto ParamFloatS = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ss", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatQ = Parameter<ArrayType<float>,Traits::SquareMatrix>("Qs", "single precision unitary matrix corresponding to Ss");


    auto gelyap_double            = makeFunctionExt<1>(call_gelyap<double>, "X", ParamDoubleA, ParamDoubleY);
    auto gelyap_double_optset     = makeFunctionExt<1>(call_gelyap_optset<double>, "X", ParamDoubleA, ParamDoubleY, ParamOptSet);
    auto gelyap_double_op         = makeFunctionExt<1>(call_gelyap_op<double>, "X", ParamDoubleA, ParamDoubleY, ParamOpString);
    auto gelyap_double_op_optset  = makeFunctionExt<1>(call_gelyap_op_optset<double>, "X", ParamDoubleA, ParamDoubleY, ParamOpString, ParamOptSet);
    auto gelyap_3ret_double            = makeFunctionExt<3>(call_3ret_gelyap<double>, "[X, S, Q]", ParamDoubleA, ParamDoubleY);
    auto gelyap_3ret_double_optset     = makeFunctionExt<3>(call_3ret_gelyap_optset<double>, "[X, S, Q]", ParamDoubleA, ParamDoubleY, ParamOptSet);
    auto gelyap_3ret_double_op         = makeFunctionExt<3>(call_3ret_gelyap_op<double>, "[X, S, Q]", ParamDoubleA, ParamDoubleY, ParamOpString);
    auto gelyap_3ret_double_op_optset  = makeFunctionExt<3>(call_3ret_gelyap_op_optset<double>, "[X, S, Q]", ParamDoubleA, ParamDoubleY, ParamOpString, ParamOptSet);
    auto gelyap_prefact_double            = makeFunctionExt<1>(call_gelyap_prefact<double>, "X",           ParamDoubleS, ParamDoubleQ, ParamDoubleY);
    auto gelyap_prefact_double_optset     = makeFunctionExt<1>(call_gelyap_prefact_optset<double>, "X",    ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOptSet);
    auto gelyap_prefact_double_op         = makeFunctionExt<1>(call_gelyap_prefact_op<double>, "X",        ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOpString);
    auto gelyap_prefact_double_op_optset  = makeFunctionExt<1>(call_gelyap_prefact_op_optset<double>, "X", ParamDoubleS, ParamDoubleQ, ParamDoubleY, ParamOpString, ParamOptSet);




    auto gelyap_float            = makeFunctionExt<1>(call_gelyap<float>, "Xs", ParamFloatA, ParamFloatY);
    auto gelyap_float_optset     = makeFunctionExt<1>(call_gelyap_optset<float>, "Xs", ParamFloatA, ParamFloatY, ParamOptSet);
    auto gelyap_float_op         = makeFunctionExt<1>(call_gelyap_op<float>, "Xs", ParamFloatA, ParamFloatY, ParamOpString);
    auto gelyap_float_op_optset  = makeFunctionExt<1>(call_gelyap_op_optset<float>, "Xs", ParamFloatA, ParamFloatY, ParamOpString, ParamOptSet);
    auto gelyap_3ret_float            = makeFunctionExt<3>(call_3ret_gelyap<float>, "[Xs, Ss, Qs]", ParamFloatA, ParamFloatY);
    auto gelyap_3ret_float_optset     = makeFunctionExt<3>(call_3ret_gelyap_optset<float>, "[Xs, Ss, Qs]", ParamFloatA, ParamFloatY, ParamOptSet);
    auto gelyap_3ret_float_op         = makeFunctionExt<3>(call_3ret_gelyap_op<float>, "[Xs, Ss, Qs]", ParamFloatA, ParamFloatY, ParamOpString);
    auto gelyap_3ret_float_op_optset  = makeFunctionExt<3>(call_3ret_gelyap_op_optset<float>, "[Xs, Ss, Qs]", ParamFloatA, ParamFloatY, ParamOpString, ParamOptSet);
    auto gelyap_prefact_float            = makeFunctionExt<1>(call_gelyap_prefact<float>, "Xs",           ParamFloatS, ParamFloatQ, ParamFloatY);
    auto gelyap_prefact_float_optset     = makeFunctionExt<1>(call_gelyap_prefact_optset<float>, "Xs",    ParamFloatS, ParamFloatQ, ParamFloatY, ParamOptSet);
    auto gelyap_prefact_float_op         = makeFunctionExt<1>(call_gelyap_prefact_op<float>, "Xs",        ParamFloatS, ParamFloatQ, ParamFloatY, ParamOpString);
    auto gelyap_prefact_float_op_optset  = makeFunctionExt<1>(call_gelyap_prefact_op_optset<float>, "Xs", ParamFloatS, ParamFloatQ, ParamFloatY, ParamOpString, ParamOptSet);





	auto parser = makeArgParser("mepack_lyap", gelyap_double,  gelyap_double_optset, gelyap_double_op, gelyap_double_op_optset,
                                                 gelyap_3ret_double,  gelyap_3ret_double_optset, gelyap_3ret_double_op, gelyap_3ret_double_op_optset,
                                                 gelyap_prefact_double, gelyap_prefact_double_optset, gelyap_prefact_double_op, gelyap_prefact_double_op_optset,
                                                 gelyap_3ret_float,  gelyap_3ret_float_optset, gelyap_3ret_float_op, gelyap_3ret_float_op_optset,
                                                 gelyap_prefact_float, gelyap_prefact_float_optset, gelyap_prefact_float_op, gelyap_prefact_float_op_optset,
                                                 gelyap_float, gelyap_float_op, gelyap_float_optset, gelyap_float_op_optset);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


