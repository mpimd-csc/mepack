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
struct sylv2_solver {
    constexpr static const char *fn = NULL;
    virtual void solve(const char * FACTA, const char *FACT,  const char *TRANSA, const char *TRANSB, T sgn,
                       int M, int N, T * A,  int LDA, T * B, int ldb, T * QA, int LDQA, T* QB, int ldza,
                       T *X, int LDX, T * SCALE, T *WORK, size_t LDWORK, int *INFO) = 0;
    virtual void blocksize_set(int nb) = 0;
    virtual void blocksize_big_set(int nb) = 0;
    virtual void isolver_set(int is) = 0;
    virtual void fsolver_set(int fs) = 0;
};

template<>
struct sylv2_solver<double> {
    constexpr static const char *fn = "DLA_GESYLV2";

    static void solve(const char * facta, const char *factb, const char *opA, const char * opB, double sgn, int M, int N ,
            double * Aptr, int lda, double *Bptr, int ldb, double *Qptr, int ldq, double *Zptr, int ldz, double *Yptr, int ldy,
            double * scale, double * buffer, size_t ldwork, int *info) {
        mepack_double_gesylv2(facta, factb, opA, opB, sgn, M, N, Aptr, lda, Bptr, ldb, Qptr, ldq, Zptr, ldz, Yptr, ldy, scale, buffer, ldwork, info);
    };
    static void blocksize_mb_set(int nb) {
        mepack_double_trsylv2_blocksize_mb_set(nb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_double_trsylv2_blocksize_nb_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_trsylv2_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_trsylv2_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_double_trsylv2_blocksize_2stage_set(nb);
    };

};

template<>
struct sylv2_solver<float> {
    constexpr static const char *fn = "SLA_GESYLV2";

    static void solve(const char * facta, const char *factb, const char *opA, const char * opB, float sgn, int M, int N ,
            float * Aptr, int lda, float *Bptr, int ldb, float *Qptr, int ldq, float *Zptr, int ldz, float *Yptr, int ldy,
            float * scale, float * buffer, size_t ldwork, int *info) {
        mepack_single_gesylv2(facta, factb, opA, opB, sgn, M, N, Aptr, lda, Bptr, ldb, Qptr, ldq, Zptr, ldz, Yptr, ldy, scale, buffer, ldwork, info);
    };
    static void blocksize_mb_set(int nb) {
        mepack_single_trsylv2_blocksize_mb_set(nb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_single_trsylv2_blocksize_nb_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_trsylv2_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_trsylv2_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_single_trsylv2_blocksize_2stage_set(nb);
    };

};


template <typename T>
static std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
SylvDriver( ArrayType<T>& A, ArrayType<T>& B, ArrayType<T>& Y, ArrayType<T>& Q, ArrayType<T>& Z, std::string opA = std::string("N"), std::string opB = std::string("N"),
        MepackOptions optset = MepackOptions(), bool factorize = false)
{

    if ( A.rows != Y.rows) {
        mexoct_error("MEPACK:DimMissMatch",
                "The first coefficient matrix and the right hand side need to have the same number of rows size(A) = [ %d %d ] size (Y) = [ %d %d ]",
                (int) A.rows, (int) A.columns, (int) Y.rows, (int) Y.columns);
    }
    if ( B.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The second coefficient matrix and the right hand side need to have the same number of columns size(B) = [ %d %d ] size (Y) = [ %d %d ]",
                (int) B.rows, (int) B.columns, (int) Y.rows, (int) Y.columns);
    }
    if ( Q.rows != A.rows || Q.columns != A.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The dimensions of the first coefficient matrix and its transformation matrix must fit. size(A) = [ %d %d ] size (Q) = [ %d %d ]",
                (int) A.rows, (int) A.columns, (int) Q.rows, (int) Q.columns);
    }
    if ( Z.rows != B.rows || Z.columns != B.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The dimensions of the second coefficient matrix and its transformation matrix must fit. size(B) = [ %d %d ] size (Z) = [ %d %d ]",
                (int) B.rows, (int) B.columns, (int) Z.rows, (int) Z.columns);
    }



    if ( optset.getMB() >= 0 ) {
        sylv2_solver<T>::blocksize_mb_set(optset.getMB());
    }
    if ( optset.getNB() >= 0 ) {
        sylv2_solver<T>::blocksize_nb_set(optset.getNB());
    }

    if ( optset.getISolver() >= 0 ) {
        sylv2_solver<T>::isolver_set(optset.getISolver());
    }

    if ( optset.getFSolver() >= 0 ) {
        sylv2_solver<T>::fsolver_set(optset.getFSolver());
    }


    if ( optset.getBigNB() >= 0 ) {
        sylv2_solver<T>::blocksize_big_set(optset.getBigNB());
    }


    ssize_t auxmem;

    if ( factorize ) {
        auxmem = mepack_memory_frontend(sylv2_solver<T>::fn, (const char *){"N"}, (const char *){"N"}, A.rows, B.rows);
    } else {
        auxmem = mepack_memory_frontend(sylv2_solver<T>::fn, (const char *){"F"}, (const char *){"F"}, A.rows, B.rows);
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


    if ( factorize ) {
        sylv2_solver<T>::solve("N", "N", opA.c_str(), opB.c_str(), sign, A.rows, B.rows, A.ptr(), A.ld, B.ptr(), B.ld, Q.ptr(), Q.ld, Z.ptr(), Z.ld, Y.ptr(), Y.ld, &scale, auxbuf, auxmem, &info);
    } else {
        sylv2_solver<T>::solve("F", "F", opA.c_str(), opB.c_str(), sign, A.rows, B.rows, A.ptr(), A.ld, B.ptr(), B.ld, Q.ptr(), Q.ld, Z.ptr(), Z.ld, Y.ptr(), Y.ld, &scale, auxbuf, auxmem, &info);
    }

    if ( info != 0) {
        mexoct_error("MEPACK:GESTEIN", "GESTEIN failed with info = %d", (int) info);
    }

    return std::make_tuple(std::move(Y), std::move(A), std::move(B), std::move(Q), std::move(Z));
}


template <typename T>
std::tuple<ArrayType<T>>
call_sylv2_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    ArrayType<T> Q(A.rows, A.rows);
    ArrayType<T> Z(B.rows, B.rows);

    auto ret = SylvDriver(A, B, Y, Q, Z, opA, opB, optset, true);
    return std::make_tuple(std::get<0>(ret));
}

template<typename T>
std::tuple<ArrayType<T>>
call_sylv2(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y)
{
    return call_sylv2_op_optset(A, B, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_sylv2_op(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, std::string opA, std::string opB)
{
    return call_sylv2_op_optset(A, B, Y, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>>
call_sylv2_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, MepackOptions optset)
{
   return call_sylv2_op_optset(A, B, Y, "N", "N", optset);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_5ret_sylv2_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, std::string opA = std::string("N"), std::string opB = std::string("N"),
                          MepackOptions optset = MepackOptions())
{
    ArrayType<T> Q(A.rows, A.rows);
    ArrayType<T> Z(B.rows, B.rows);
    return SylvDriver(A, B, Y, Q, Z, opA, opB, optset, true);
}



template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>,ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_5ret_sylv2(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y)
{

    return call_5ret_sylv2_op_optset(A, B, Y);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_5ret_sylv2_op(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, std::string opA, std::string opB)
{
    return call_5ret_sylv2_op_optset(A, B, Y, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_5ret_sylv2_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> Y, MepackOptions optset)
{
   return call_5ret_sylv2_op_optset(A, B, Y, "N", "N", optset);
}

template <typename T>
std::tuple<ArrayType<T>>
call_sylv2_prefact_op_optset(ArrayType<T>  S, ArrayType<T> Ts, ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> Y,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    bool chk = mepack_mexoct_check_qtriangular(S);
    bool chk2 = mepack_mexoct_check_qtriangular(Ts);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The first coefficient matrix is not (quasi) upper triangular.");
    }
    if ( !chk2 ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The second coefficient matrix is not (quasi) upper triangular.");
    }


    auto ret = SylvDriver(S, Ts, Y, Q, Z, opA, opB, optset, false);
    return std::make_tuple(std::get<0>(ret));
}

template<typename T>
std::tuple<ArrayType<T>>
call_sylv2_prefact(ArrayType<T> S, ArrayType<T> Ts, ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> Y)
{
    return call_sylv2_prefact_op_optset(S, Ts, Q, Z, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_sylv2_prefact_op(ArrayType<T>  S, ArrayType<T> Ts, ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> Y, std::string opA, std::string opB)
{
    return call_sylv2_prefact_op_optset(S, Ts, Q, Z, Y, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>>
call_sylv2_prefact_optset(ArrayType<T>  S, ArrayType<T> Ts, ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> Y, MepackOptions optset)
{
   return call_sylv2_prefact_op_optset(S, Ts, Q, Z, Y, "N", "N", optset);
}



MEXOCT_ENTRY(mepack_sylv2,
R"===(
X = mepack_sylv2(A, B, Y)
X = mepack_sylv2(A, B, Y, OpA, OpB)
X = mepack_sylv2(A, B, Y, OptSet)
X = mepack_sylv2(A, B, Y, OpA, OpB, OptSet)
 Input:
    A      -- first coefficient matrix
    B      -- second coefficient matrix
    Y      -- right hand side
    OpA    -- Transpose operation on the left coefficient matrix
    OpB    -- Transpose operation on the right coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.
 Return:
    X      -- Solution of the Sylvester equation


[X, S, T, Q, Z] = mepack_sylv2(A, B, Y)
[X, S, T, Q, Z] = mepack_sylv2(A, B, Y, OpA, OpB)
[X, S, T, Q, Z] = mepack_sylv2(A, B, Y, OptSet)
[X, S, T, Q, Z] = mepack_sylv2(A, B, Y, OpA, OpB, OptSet)
  Input:
    A      -- first coefficient matrix
    B      -- second coefficient matrix
    Y      -- right hand side
    OpA    -- Transpose operation on the left coefficient matrix
    OpB    -- Transpose operation on the right coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

  Return:
    X      -- Solution of the Sylvester  equation
    S, Q   -- Schur decomposition of A, A = Q * S * Q'
    T, Z   -- Schur decomposition of B, B = Z * T * B'

X = mepack_sylv2(S, T, Q, Z, Y)
X = mepack_sylv2(S, T, Q, Z, Y, OpA, OpB)
X = mepack_sylv2(S, T, Q, Z, Y, OptSet)
X = mepack_sylv2(S, T, Q, Z, Y, OpA, OpB, OptSet)
  Input:
    S      -- (quasi) upper triangular coefficient matrix
    T      -- (quasi) upper triangular coefficient matrix
    Q      -- unitary transformation matrix, A = Q*S*Q'
    Z      -- unitary transformation matrix, B = Z*T*Z'
    Y      -- right hand side
    OpA    -- Transpose operation on the left coefficient matrix
    OpB    -- Transpose operation on the right coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

  Return:
    X      -- Solution of the Sylvester equation

The mepack_sylv2 function solves the Sylvester equation

   opA(A) * X * opB(B) + X = Y,

or

   opA(A) * X * opB(B) - X = Y,

where A and B are general matrices. The operation opA (or opB) either returns the matrix or
if OpA == 'T' (or opB == 'T') its transpose. In this way, the transposed equation
can be solved without transposing the coefficient matrix before. The function
can also be called with A and B decomposed into their Schur decompositions,
A = Q*S*Q' and B = Q*T*Q', where S is the Schur form of A and T is the Schur form
of B. If the function is called with five return values,
the Schur forms of A and B and their unitary transformation matrices are returned.

The function uses the level-3 solver D/SLA_GESYLV2 of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_sylv2

   - mb      -- row block size of the level-3 algorithm (mb >=2)
   - bignb   -- block size of two-stage solver (bignb >=2)
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

    auto ParamOpAString = Parameter<std::string>("OpA", "Transpose operation on the first coefficient matrix");
    auto ParamOpBString = Parameter<std::string>("OpB", "Transpose operation on the second coefficient matrix");
    auto ParamOptSet   = Parameter<MepackOptions>("optset", "Options for MEPACK");

    auto ParamDoubleA = Parameter<ArrayType<double>,Traits::SquareMatrix>("A", "double precision coefficient matrix");
    auto ParamDoubleB = Parameter<ArrayType<double>,Traits::SquareMatrix>("B", "double precision coefficient matrix");
    auto ParamDoubleS = Parameter<ArrayType<double>,Traits::SquareMatrix>("S", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleT = Parameter<ArrayType<double>,Traits::SquareMatrix>("T", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleQ = Parameter<ArrayType<double>,Traits::SquareMatrix>("Q", "double precision unitary matrix corresponding to A and S");
    auto ParamDoubleZ = Parameter<ArrayType<double>,Traits::SquareMatrix>("Z", "double precision unitary matrix corresponding to B and T)");
    auto ParamDoubleY = Parameter<ArrayType<double>>("Y", "double precision right hand side");

    auto ParamFloatA = Parameter<ArrayType<float>, Traits::SquareMatrix>("As", "single precision coefficient matrix");
    auto ParamFloatB = Parameter<ArrayType<float>, Traits::SquareMatrix>("Bs", "single precision coefficient matrix");
    auto ParamFloatS = Parameter<ArrayType<float>, Traits::SquareMatrix>("Ss", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatT = Parameter<ArrayType<float>, Traits::SquareMatrix>("Ts", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatQ = Parameter<ArrayType<float>, Traits::SquareMatrix>("Qs", "single precision unitary matrix corresponding to As and Ss");
    auto ParamFloatZ = Parameter<ArrayType<float>, Traits::SquareMatrix>("Zs", "single precision unitary matrix corresponding to Bs and Ts");
    auto ParamFloatY = Parameter<ArrayType<float>>("Ys", "single precision right hand side");


    auto sylv2_double                    = makeFunctionExt<1>(call_sylv2<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY);
    auto sylv2_double_optset             = makeFunctionExt<1>(call_sylv2_optset<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOptSet);
    auto sylv2_double_op                 = makeFunctionExt<1>(call_sylv2_op<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpAString, ParamOpBString);
    auto sylv2_double_op_optset          = makeFunctionExt<1>(call_sylv2_op_optset<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpAString, ParamOpBString, ParamOptSet);
    auto sylv2_3ret_double               = makeFunctionExt<5>(call_5ret_sylv2<double>, "[X, S, T, Q, Z]", ParamDoubleA, ParamDoubleB, ParamDoubleY);
    auto sylv2_3ret_double_optset        = makeFunctionExt<5>(call_5ret_sylv2_optset<double>, "[X, S, T, Q, Z]", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOptSet);
    auto sylv2_3ret_double_op            = makeFunctionExt<5>(call_5ret_sylv2_op<double>, "[X, S, T, Q, Z]", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpAString, ParamOpBString);
    auto sylv2_3ret_double_op_optset     = makeFunctionExt<5>(call_5ret_sylv2_op_optset<double>, "[X, S, T, Q, Z]", ParamDoubleA, ParamDoubleB, ParamDoubleY, ParamOpAString,
                                                                                                                    ParamOpBString, ParamOptSet);
    auto sylv2_prefact_double            = makeFunctionExt<1>(call_sylv2_prefact<double>, "X",           ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleZ, ParamDoubleY);
    auto sylv2_prefact_double_optset     = makeFunctionExt<1>(call_sylv2_prefact_optset<double>, "X",    ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleZ, ParamDoubleY, ParamOptSet);
    auto sylv2_prefact_double_op         = makeFunctionExt<1>(call_sylv2_prefact_op<double>, "X",        ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleZ, ParamDoubleY,
                                                                                                         ParamOpAString, ParamOpBString);
    auto sylv2_prefact_double_op_optset  = makeFunctionExt<1>(call_sylv2_prefact_op_optset<double>, "X", ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleZ, ParamDoubleY,
                                                                                                         ParamOpAString, ParamOpBString, ParamOptSet);

    auto sylv2_float                    = makeFunctionExt<1>(call_sylv2<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY);
    auto sylv2_float_optset             = makeFunctionExt<1>(call_sylv2_optset<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOptSet);
    auto sylv2_float_op                 = makeFunctionExt<1>(call_sylv2_op<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpAString, ParamOpBString);
    auto sylv2_float_op_optset          = makeFunctionExt<1>(call_sylv2_op_optset<float>, "Xs", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpAString, ParamOpBString, ParamOptSet);
    auto sylv2_3ret_float               = makeFunctionExt<5>(call_5ret_sylv2<float>, "[Xs, Ss, Ts, Qs, Zs]", ParamFloatA, ParamFloatB,  ParamFloatY);
    auto sylv2_3ret_float_optset        = makeFunctionExt<5>(call_5ret_sylv2_optset<float>, "[Xs, Ss, Ts, Qs, Zs]", ParamFloatA, ParamFloatB, ParamFloatY, ParamOptSet);
    auto sylv2_3ret_float_op            = makeFunctionExt<5>(call_5ret_sylv2_op<float>, "[Xs, Ss, Ts, Qs, Zs]", ParamFloatA, ParamFloatB, ParamFloatY, ParamOpAString, ParamOpBString);
    auto sylv2_3ret_float_op_optset     = makeFunctionExt<5>(call_5ret_sylv2_op_optset<float>, "[Xs, Ss, Ts, Qs, Zs]", ParamFloatA, ParamFloatB, ParamFloatY,
                                                                                                      ParamOpAString, ParamOpBString, ParamOptSet);
    auto sylv2_prefact_float            = makeFunctionExt<1>(call_sylv2_prefact<float>, "Xs",           ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatZ, ParamFloatY);
    auto sylv2_prefact_float_optset     = makeFunctionExt<1>(call_sylv2_prefact_optset<float>, "Xs",    ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatZ, ParamFloatY, ParamOptSet);
    auto sylv2_prefact_float_op         = makeFunctionExt<1>(call_sylv2_prefact_op<float>, "Xs",        ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatZ, ParamFloatY, ParamOpAString,
                                                                                                      ParamOpBString);
    auto sylv2_prefact_float_op_optset  = makeFunctionExt<1>(call_sylv2_prefact_op_optset<float>, "Xs", ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatZ, ParamFloatY,
                                                                                                      ParamOpAString, ParamOpBString, ParamOptSet);





	auto parser = makeArgParser("mepack_sylv2", sylv2_double,  sylv2_double_optset, sylv2_double_op, sylv2_double_op_optset,
                                                 sylv2_3ret_double,  sylv2_3ret_double_optset, sylv2_3ret_double_op, sylv2_3ret_double_op_optset,
                                                 sylv2_prefact_double, sylv2_prefact_double_optset, sylv2_prefact_double_op, sylv2_prefact_double_op_optset,
                                                 sylv2_3ret_float,  sylv2_3ret_float_optset, sylv2_3ret_float_op, sylv2_3ret_float_op_optset,
                                                 sylv2_prefact_float, sylv2_prefact_float_optset, sylv2_prefact_float_op, sylv2_prefact_float_op_optset,
                                                 sylv2_float, sylv2_float_op, sylv2_float_optset, sylv2_float_op_optset);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


