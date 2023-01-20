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
#include "triangular_test.hpp"


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
struct gstein_solver {
    constexpr static const char *fn = NULL;
    virtual void solve(const char * FACT, const char *TRANSA, int M, T * A,  int LDA, T * B, int ldb, T * Q, int LDQ, T* Z, int ldz,
                       T *X, int LDX, T * SCALE, T *WORK, size_t LDWORK, int *INFO) = 0;
    virtual void blocksize_set(int nb) = 0;
    virtual void blocksize_big_set(int nb) = 0;
    virtual void isolver_set(int is) = 0;
    virtual void fsolver_set(int fs) = 0;
};

template<>
struct gstein_solver<double> {
    constexpr static const char *fn = "DLA_GGSTEIN";

    static void solve(const char * fact, const char *op, int M, double * Aptr, int lda, double *Bptr, int ldb, double *Qptr, int ldq, double *Zptr, int ldz, double *Yptr, int ldy, double * scale, double * buffer, size_t ldwork, int *info) {
        mepack_double_ggstein(fact, op, M, Aptr, lda, Bptr, ldb, Qptr, ldq, Zptr, ldz, Yptr, ldy, scale, buffer, ldwork, info);
    };
    static void blocksize_set(int nb) {
        mepack_double_tgstein_blocksize_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_tgstein_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_tgstein_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_double_tgstein_blocksize_2stage_set(nb);
    };

};

template<>
struct gstein_solver<float> {
    constexpr static const char *fn = "SLA_GGSTEIN";

    static void solve(const char * fact, const char *op, int M, float * Aptr, int lda, float *Bptr, int ldb, float *Qptr, int ldq, float * Zptr, int ldz,
                      float *Yptr, int ldy, float * scale, float * buffer, size_t ldwork, int *info) {
        mepack_single_ggstein(fact, op, M, Aptr, lda, Bptr, ldb, Qptr, ldq, Zptr, ldz, Yptr, ldy, scale, buffer, ldwork, info);
    };
    static void blocksize_set(int nb) {
        mepack_single_tgstein_blocksize_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_tgstein_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_tgstein_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_single_tgstein_blocksize_2stage_set(nb);
    };
};


template <typename T>
static std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
GSteinDriver( ArrayType<T>& A, ArrayType<T>& B, ArrayType<T>& Y, ArrayType<T>& Q, ArrayType<T>& Z, std::string opA = std::string("N"), MepackOptions optset = MepackOptions(), bool factorize = false)
{

    if (A.rows != Y.rows || A.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The first coefficient matrix and the right hand side have to be of the same dimension size(A) = [ %d %d ] size (Y) = [ %d %d ]", (int) A.rows, (int) A.columns,
                (int) Y.rows, (int) Y.columns);
    }
    if (B.rows != Y.rows || B.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The second coefficient matrix and the right hand side have to be of the same dimension size(B) = [ %d %d ] size (Y) = [ %d %d ]", (int) B.rows, (int) B.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (Q.rows != Y.rows || Q.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The first unitary transformation matrix and the right hand side have to be of the same dimension size(Q) = [ %d %d ] size (Y) = [ %d %d ]", (int) Q.rows, (int) Q.columns,
                (int) Y.rows, (int) Y.columns);
    }
    if (Z.rows != Y.rows || Z.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The second unitary transformation matrix and the right hand side have to be of the same dimension size(Z) = [ %d %d ] size (Y) = [ %d %d ]", (int) Z.rows, (int) Z.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if ( optset.getMB() >= 0 ) {
        gstein_solver<T>::blocksize_set(optset.getMB());
    }
    if ( optset.getISolver() >= 0 ) {
        gstein_solver<T>::isolver_set(optset.getISolver());
    }
    if ( optset.getFSolver() >= 0 ) {
        gstein_solver<T>::fsolver_set(optset.getFSolver());
    }

    ssize_t auxmem;

    if ( factorize ) {
        auxmem = mepack_memory_frontend(gstein_solver<T>::fn, (const char *){"N"}, (const char *){"N"}, A.rows, A.columns);
    } else {
        auxmem = mepack_memory_frontend(gstein_solver<T>::fn, (const char *){"F"}, (const char *){"F"}, A.rows, A.columns);
    }

    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    T scale = 1;
    int info = 0;

    if ( factorize ) {
        gstein_solver<T>::solve("N", opA.c_str(), A.rows, A.ptr(), A.ld, B.ptr(), B.ld, Q.ptr(), Q.ld, Z.ptr(), Z.ld, Y.ptr(), Y.ld, &scale, auxbuf, auxmem, &info);
    } else {
        gstein_solver<T>::solve("F", opA.c_str(), A.rows, A.ptr(), A.ld, B.ptr(), B.ld, Q.ptr(), Q.ld, Z.ptr(), Z.ld, Y.ptr(), Y.ld, &scale, auxbuf, auxmem, &info);
    }

    if ( info != 0) {
        mexoct_error("MEPACK:GGSTEIN", "GGLYAP failed with info = %d", (int) info);
    }

    return std::make_tuple(std::move(Y), std::move(A), std::move(B), std::move(Q), std::move(Z));
}


template <typename T>
std::tuple<ArrayType<T>>
call_ggstein_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    ArrayType<T> Q(A.rows, A.rows);
    ArrayType<T> Z(A.rows, A.rows);

    auto ret = GSteinDriver(A, B, Y, Q, Z, opA, optset, true);
    return std::make_tuple(std::get<0>(ret));
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggstein(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y)
{
    return call_ggstein_op_optset(A, B, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggstein_op(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, std::string opA)
{
    return call_ggstein_op_optset(A, B, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggstein_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, MepackOptions optset)
{
   return call_ggstein_op_optset(A, B, Y, "N", optset);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_5ret_ggstein_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    ArrayType<T> Q(A.rows, A.rows);
    ArrayType<T> Z(A.rows, A.rows);
    return GSteinDriver(A, B, Y, Q, Z, opA, optset, true);
}



template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>,ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_5ret_ggstein(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y)
{

    return call_5ret_ggstein_op_optset(A, B, Y);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_5ret_ggstein_op(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, std::string opA)
{
    return call_5ret_ggstein_op_optset(A, B, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_5ret_ggstein_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, MepackOptions optset)
{
   return call_5ret_ggstein_op_optset(A, B, Y, "N", optset);
}

template <typename T>
std::tuple<ArrayType<T>>
call_ggstein_prefact_op_optset(ArrayType<T>  S, ArrayType<T> Ts, ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> Y, std::string opA = std::string("N"), MepackOptions optset = MepackOptions())
{
    bool chk = mepack_mexoct_check_qtriangular(S);
    bool chk2 = mepack_mexoct_check_triangular(Ts);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The first coefficient matrix is not (quasi) upper triangular.");
    }
    if ( !chk2 ) {
        mexoct_error("MEPACK:UpperTriangular", "The second coefficient matrix is not upper triangular.");
    }


    auto ret = GSteinDriver(S, Ts, Y, Q, Z, opA, optset, false);
    return std::make_tuple(std::get<0>(ret));
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggstein_prefact(ArrayType<T> S, ArrayType<T> Ts, ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> Y)
{
    return call_ggstein_prefact_op_optset(S, Ts, Q, Z, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggstein_prefact_op(ArrayType<T>  S, ArrayType<T> Ts, ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> Y, std::string opA)
{
    return call_ggstein_prefact_op_optset(S, Ts, Q, Z, Y, opA);
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggstein_prefact_optset(ArrayType<T>  S, ArrayType<T> Ts, ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> Y, MepackOptions optset)
{
   return call_ggstein_prefact_op_optset(S, Ts, Q, Z, Y, "N", optset);
}



MEXOCT_ENTRY(mepack_gstein,
R"===(
X = mepack_Stein(A, B, Y)
X = mepack_Stein(A, B, Y, Op)
X = mepack_gstein(A, B, Y, OptSet)
X = mepack_gstein(A, B, Y, Op, OptSet)
 Input:
    A      -- first coefficient matrix
    B      -- second coefficient matrix
    Y      -- right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.
 Return:
    X      -- Solution of the generalized Stein equation


[X, S, T, Q, Z] = mepack_gstein(A, B, Y)
[X, S, T, Q, Z] = mepack_gstein(A, B, Y, Op)
[X, S, T, Q, Z] = mepack_gstein(A, B, Y, OptSet)
[X, S, T, Q, Z] = mepack_gstein(A, B, Y, Op, OptSet)
  Input:
    A      -- first coefficient matrix
    B      -- second coefficient matrix
    Y      -- right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

  Return:
    X      -- Solution of the generalized Stein equation
    S, T,
    Q, Z   -- generalized Schur decomposition of (A,B),
              A = Q*S*Z', B = Q*T*Z'

X = mepack_gstein(S, T, Q, Z, Y)
X = mepack_gstein(S, T, Q, Z, Y, Op)
X = mepack_gstein(S, T, Q, Z, Y, OptSet)
X = mepack_gstein(S, T, Q, Z, Y, Op, OptSet)
  Input:
    S      -- (quasi) upper triangular coefficient matrix
    T      -- upper triangular coefficient matrix
    Q      -- unitary transformation matrix, A = Q*S*Z'
    Z      -- unitary transformation matrix, B = Q*T*Z'
    Y      -- right hand side
    Op     -- Transpose operation on the coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

  Return:
    X      -- Solution of the generalized Stein equation

The mepack_gstein function solves the generalized Lyapunov equation

   op(A) * X * op(A)' - op(B) * X * op(A)' = Y,

where A and B are general matrices. The operation op either returns the matrix or
if Op == 'T' its transpose. In this way, the transposed equation
can be solved without transposing the coefficient matrix before. The function
can also be called with (A,B) decomposed in to A = Q*S*Z' and B = Q*T*Z',
where (S,T) is the generalized Schur form of (A, B).
If the function is called with five return values,
the generalized Schur form of (A,B) and its unitary transformation matrices are returned.

The matrices A, B, S, T, Q, Z  and Y must be of the same size and Y must be symmetric.

The function uses the level-3 solver D/SLA_GGSTEIN of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_gstein:

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
    auto ParamDoubleB = Parameter<ArrayType<double>,Traits::SquareMatrix>("B", "double precision coefficient matrix");
    auto ParamDoubleS = Parameter<ArrayType<double>,Traits::SquareMatrix>("S", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleT = Parameter<ArrayType<double>,Traits::SquareMatrix>("T", "double precision upper triangular coefficient matrix");
    auto ParamDoubleQ = Parameter<ArrayType<double>,Traits::SquareMatrix>("Q", "double precision unitary matrix corresponding to (S,T)");
    auto ParamDoubleZ = Parameter<ArrayType<double>,Traits::SquareMatrix>("Z", "double precision unitary matrix corresponding to (S,T)");
    auto ParamDoubleY = Parameter<ArrayType<double>,Traits::SquareMatrix>("Y", "double precision right hand side");

    auto ParamFloatA = Parameter<ArrayType<float>,Traits::SquareMatrix>("As", "single precision coefficient matrix");
    auto ParamFloatB = Parameter<ArrayType<float>,Traits::SquareMatrix>("Bs", "single precision coefficient matrix");
    auto ParamFloatY = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ys", "single precision right hand side");
    auto ParamFloatS = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ss", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatT = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ts", "single precision upper triangular coefficient matrix");
    auto ParamFloatQ = Parameter<ArrayType<float>,Traits::SquareMatrix>("Qs", "single precision unitary matrix corresponding to (Ss,Ts)");
    auto ParamFloatZ = Parameter<ArrayType<float>,Traits::SquareMatrix>("Zs", "single precision unitary matrix corresponding to (Ss,Ts)");



    auto ggstein_double                    = makeFunctionExt<1>(call_ggstein<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY);
    auto ggstein_double_optset             = makeFunctionExt<1>(call_ggstein_optset<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOptSet);
    auto ggstein_double_op                 = makeFunctionExt<1>(call_ggstein_op<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpString);
    auto ggstein_double_op_optset          = makeFunctionExt<1>(call_ggstein_op_optset<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpString, ParamOptSet);
    auto ggstein_3ret_double               = makeFunctionExt<5>(call_5ret_ggstein<double>, "[X, S, T, Q, Z]", ParamDoubleA, ParamDoubleB, ParamDoubleY);
    auto ggstein_3ret_double_optset        = makeFunctionExt<5>(call_5ret_ggstein_optset<double>, "[X, S, T, Q, Z]", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOptSet);
    auto ggstein_3ret_double_op            = makeFunctionExt<5>(call_5ret_ggstein_op<double>, "[X, S, T, Q, Z]", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpString);
    auto ggstein_3ret_double_op_optset     = makeFunctionExt<5>(call_5ret_ggstein_op_optset<double>, "[X, S, T, Q, Z]", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpString, ParamOptSet);
    auto ggstein_prefact_double            = makeFunctionExt<1>(call_ggstein_prefact<double>, "X",           ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleZ, ParamDoubleY);
    auto ggstein_prefact_double_optset     = makeFunctionExt<1>(call_ggstein_prefact_optset<double>, "X",    ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleZ, ParamDoubleY, ParamOptSet);
    auto ggstein_prefact_double_op         = makeFunctionExt<1>(call_ggstein_prefact_op<double>, "X",        ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleZ, ParamDoubleY, ParamOpString);
    auto ggstein_prefact_double_op_optset  = makeFunctionExt<1>(call_ggstein_prefact_op_optset<double>, "X", ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleZ, ParamDoubleY, ParamOpString, ParamOptSet);




    auto ggstein_float                    = makeFunctionExt<1>(call_ggstein<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY);
    auto ggstein_float_optset             = makeFunctionExt<1>(call_ggstein_optset<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOptSet);
    auto ggstein_float_op                 = makeFunctionExt<1>(call_ggstein_op<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpString);
    auto ggstein_float_op_optset          = makeFunctionExt<1>(call_ggstein_op_optset<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpString, ParamOptSet);
    auto ggstein_3ret_float               = makeFunctionExt<5>(call_5ret_ggstein<float>, "[Xs, Ss, Ts, Qs, Zs]", ParamFloatA, ParamFloatB,  ParamFloatY);
    auto ggstein_3ret_float_optset        = makeFunctionExt<5>(call_5ret_ggstein_optset<float>, "[Xs, Ss, Ts, Qs, Zs]", ParamFloatA, ParamFloatB, ParamFloatY, ParamOptSet);
    auto ggstein_3ret_float_op            = makeFunctionExt<5>(call_5ret_ggstein_op<float>, "[Xs, Ss, Ts, Qs, Zs]", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpString);
    auto ggstein_3ret_float_op_optset     = makeFunctionExt<5>(call_5ret_ggstein_op_optset<float>, "[Xs, Ss, Ts, Qs, Zs]", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpString, ParamOptSet);
    auto ggstein_prefact_float            = makeFunctionExt<1>(call_ggstein_prefact<float>, "Xs",           ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatZ, ParamFloatY);
    auto ggstein_prefact_float_optset     = makeFunctionExt<1>(call_ggstein_prefact_optset<float>, "Xs",    ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatZ, ParamFloatY, ParamOptSet);
    auto ggstein_prefact_float_op         = makeFunctionExt<1>(call_ggstein_prefact_op<float>, "Xs",        ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatZ, ParamFloatY, ParamOpString);
    auto ggstein_prefact_float_op_optset  = makeFunctionExt<1>(call_ggstein_prefact_op_optset<float>, "Xs", ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatZ, ParamFloatY, ParamOpString, ParamOptSet);





	auto parser = makeArgParser("mepack_gstein", ggstein_double,  ggstein_double_optset, ggstein_double_op, ggstein_double_op_optset,
                                                 ggstein_3ret_double,  ggstein_3ret_double_optset, ggstein_3ret_double_op, ggstein_3ret_double_op_optset,
                                                 ggstein_prefact_double, ggstein_prefact_double_optset, ggstein_prefact_double_op, ggstein_prefact_double_op_optset,
                                                 ggstein_3ret_float,  ggstein_3ret_float_optset, ggstein_3ret_float_op, ggstein_3ret_float_op_optset,
                                                 ggstein_prefact_float, ggstein_prefact_float_optset, ggstein_prefact_float_op, ggstein_prefact_float_op_optset,
                                                 ggstein_float, ggstein_float_op, ggstein_float_optset, ggstein_float_op_optset);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


