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
struct sylv2_refine_solver {
    constexpr static const char *fn = NULL;
    virtual void solve (const char * TRANSA, const char * TRANSB, const char *GUESS, T SGN, int M, int N , T * A, int LDA,
        T *B, int LDB, T *X, int LDX, T *Y, int LDY, T *AS, int LDAS, T *BS, int LDBS,
        T * Q, int LDQ, T *U, int LDU,
        int *MAXIT, T *TAU, T *CONVLOG, T * WORK, size_t LWORK, int *INFO);
    virtual void blocksize_set(int nb) = 0;
    virtual void blocksize_big_set(int nb) = 0;
    virtual void isolver_set(int is) = 0;
    virtual void fsolver_set(int fs) = 0;
};

template<>
struct sylv2_refine_solver<double> {
    constexpr static const char *fn = "DLA_GESYLV2_REFINE";

    static void solve (const char * TRANSA, const char * TRANSB, const char *GUESS, double SGN, int M, int N , double * A, int LDA,
        double *B, int LDB, double *X, int LDX, double *Y, int LDY, double *AS, int LDAS, double *BS, int LDBS,
        double * Q, int LDQ, double *U, int LDU,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LWORK, int *INFO)
    {
        mepack_double_gesylv2_refine(TRANSA, TRANSB, GUESS, SGN, M, N, A, LDA, B, LDB, X, LDX, Y, LDY, AS, LDAS, BS, LDBS, Q, LDQ, U, LDU, MAXIT, TAU, CONVLOG, WORK, LWORK, INFO);
    }

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
struct sylv2_refine_solver<float> {
    constexpr static const char *fn = "SLA_GESYLV2_REFINE";

    static void solve (const char * TRANSA, const char * TRANSB, const char *GUESS, float SGN, int M, int N , float * A, int LDA,
        float *B, int LDB, float *X, int LDX, float *Y, int LDY, float *AS, int LDAS, float *BS, int LDBS,
        float * Q, int LDQ, float *U, int LDU,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LWORK, int *INFO)
    {
        mepack_single_gesylv2_refine(TRANSA, TRANSB, GUESS, SGN, M, N, A, LDA, B, LDB, X, LDX, Y, LDY, AS, LDAS, BS, LDBS, Q, LDQ, U, LDU, MAXIT, TAU, CONVLOG, WORK, LWORK, INFO);
    }

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
static std::tuple<ArrayType<T>, ArrayType<T>>
Sylv2Driver(ArrayType<T> A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT, ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X,
        std::string opA = std::string("N"), std::string opB = std::string("N"), std::string guess = std::string("N"),
        MepackOptions optset = MepackOptions())
{
    if (A.rows != Y.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix A and the right hand side must have the same number of rows size(A) = [ %d %d ] size (Y) = [ %d %d ]", (int) A.rows, (int) A.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (B.columns != Y.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix B and the right hand side must have the same number of columns size(B) = [ %d %d ] size (Y) = [ %d %d ]", (int) B.rows, (int) B.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (S.rows != Y.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix S and the right hand side must have the same number of rows size(S) = [ %d %d ] size (Y) = [ %d %d ]", (int) S.rows, (int) S.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (TT.columns != Y.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix T and the right hand side must have the same number of columns size(T) = [ %d %d ] size (Y) = [ %d %d ]", (int) TT.rows, (int) TT.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (X.rows != Y.rows || X.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The initial guess matrix and the right hand side have to be of the same dimension. size(X0) = [ %d %d ] size (Y) = [ %d %d ]", (int) X.rows, (int) X.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (Q.rows != Y.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The unitary transformation matrix Q and the right hand side must have the same number of rows. size(Q) = [ %d %d ] size (Y) = [ %d %d ]", (int) Q.rows, (int) Q.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (U.columns != Y.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The unitary transformation matrix U and the right hand side must have the same number of columns. size(U) = [ %d %d ] size (Y) = [ %d %d ]", (int) U.rows, (int) U.columns,
                (int) Y.rows, (int) Y.columns);
    }


    if ( optset.getMB() >= 0 ) {
        sylv2_refine_solver<T>::blocksize_mb_set(optset.getMB());
    }
    if ( optset.getNB() >= 0 ) {
        sylv2_refine_solver<T>::blocksize_nb_set(optset.getNB());
    }

    if ( optset.getISolver() >= 0 ) {
        sylv2_refine_solver<T>::isolver_set(optset.getISolver());
    }

    if ( optset.getBigNB() >= 0 ) {
        sylv2_refine_solver<T>::blocksize_big_set(optset.getBigNB());
    }

    if ( optset.getFSolver() >= 0 ) {
        sylv2_refine_solver<T>::fsolver_set(optset.getFSolver());
    }

    if ( optset.getMaxit() < 1) {
        mexoct_error("MEPACK:MaxitInvalid", "The maximum number of iterations is invalid. maxit = %d", (int) optset.getMaxit());
    }

    if ( optset.getTau() < 0 ) {
        mexoct_error("MEPACK:TolInvalid", "The tolerance is invalid. tol = %20.15e", (double) optset.getTau());
    }

    int maxit = (mexoct_int) optset.getMaxit();
    T  tol = (T) optset.getTau();

    T sign = 1;

    if ( optset.getSign() < 0 ) {
        sign = -1;
    }



    ssize_t auxmem;
    auxmem = mepack_memory_frontend(sylv2_refine_solver<T>::fn, (const char *){"F"}, (const char *){"F"}, A.rows, B.rows);
    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    ArrayType<T> convlog_internal(maxit, 1, false);

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    int info = 0;

    sylv2_refine_solver<T>::solve(opA.c_str(), opB.c_str(), guess.c_str(), sign, A.rows, B.rows, A.ptr(), A.ld, B.ptr(), B.ld,
                                 X.ptr(), X.ld, Y.ptr(), Y.ld,
                                 S.ptr(), S.ld, TT.ptr(), TT.ld, Q.ptr(), Q.ld, U.ptr(), U.ld,
                                 &maxit, &tol, convlog_internal.ptr(), auxbuf, auxmem, &info);

    ArrayType<T> convlog( maxit, 1, false);
    for ( auto i = 0; i < maxit; i++) convlog(i, 0) = convlog_internal(i, 0);

    if ( info != 0) {
        mexoct_error("MEPACK:GESYLV2_REFINE", "GESYLV2_REFINE failed with info = %d", (int) info);
    }
    return std::make_tuple(X, convlog) ;
}


/* Master routines */
template <typename T>
std::tuple<ArrayType<T>>
call_gesylv2_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT, ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{

    bool chk = mepack_mexoct_check_qtriangular(S);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix S is not (quasi) upper triangular.");
    }
    bool chk2 = mepack_mexoct_check_qtriangular(TT);
    if ( !chk2 ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix T is not (quasi) upper triangular.");
    }


    ArrayType<T> X(A.rows, B.rows);
    auto ret = Sylv2Driver(A, B, S, TT, Q, U, Y, X, opA, opB, "N", optset);
    return std::make_tuple(std::get<0>(ret));
}

template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gesylv2_op_optset_convlog(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT, ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    bool chk = mepack_mexoct_check_qtriangular(S);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix S is not (quasi) upper triangular.");
    }
    bool chk2 = mepack_mexoct_check_qtriangular(TT);
    if ( !chk2 ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix T is not (quasi) upper triangular.");
    }

    ArrayType<T> X(A.rows, B.rows);
    auto ret = Sylv2Driver(A, B, S, TT, Q, U, Y, X, opA, opB, "N", optset);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret));
}

template <typename T>
std::tuple<ArrayType<T>>
call_gesylv2_x0_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    bool chk = mepack_mexoct_check_qtriangular(S);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix S is not (quasi) upper triangular.");
    }
    bool chk2 = mepack_mexoct_check_qtriangular(TT);
    if ( !chk2 ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix T is not (quasi) upper triangular.");
    }

    auto ret = Sylv2Driver(A, B, S, TT, Q, U, Y, X, opA, opB, "I", optset);
    return std::make_tuple(std::get<0>(ret));
}

template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gesylv2_x0_op_optset_convlog(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{

    bool chk = mepack_mexoct_check_qtriangular(S);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix S is not (quasi) upper triangular.");
    }
    bool chk2 = mepack_mexoct_check_qtriangular(TT);
    if ( !chk2 ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix T is not (quasi) upper triangular.");
    }

    auto ret = Sylv2Driver(A, B, S, TT, Q, U, Y, X, opA, opB, "I", optset);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret));
}



/* Shorter Interfaces */

template<typename T>
std::tuple<ArrayType<T>>
call_gesylv2(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y)
{
    return call_gesylv2_op_optset(A, B, S, TT, Q, U, Y);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gesylv2_convlog(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S,  ArrayType<T> TT, ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y)
{
    return call_gesylv2_op_optset_convlog(A, B, S, TT, Q, U, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gesylv2_x0(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X)
{
    return call_gesylv2_x0_op_optset(A, B, S, TT, Q, U, Y, X);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gesylv2_x0_convlog(ArrayType<T>  A,  ArrayType<T> B,ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X)
{
    return call_gesylv2_x0_op_optset_convlog(A, B, S, TT, Q, U, Y, X);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gesylv2_op(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, std::string opA, std::string opB)
{
    return call_gesylv2_op_optset(A, B, S, TT, Q, U, Y, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gesylv2_op_convlog(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, std::string opA, std::string opB)
{
    return call_gesylv2_op_optset_convlog(A, B, S, TT, Q, U, Y, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gesylv2_x0_op(ArrayType<T>  A, ArrayType<T> B,ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X, std::string opA, std::string opB)
{
    return call_gesylv2_x0_op_optset(A, B, S, TT, Q, U, Y, X, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gesylv2_x0_op_convlog(ArrayType<T>  A,  ArrayType<T> B,ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X, std::string opA, std::string opB)
{
    return call_gesylv2_x0_op_optset_convlog(A, B, S, TT, Q, U,  Y, X, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gesylv2_opt(ArrayType<T>  A,  ArrayType<T> B,ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, MepackOptions opt)
{
    return call_gesylv2_op_optset(A, B, S, TT, Q, U, Y, "N", "N",  opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gesylv2_opt_convlog(ArrayType<T>  A,  ArrayType<T> B,ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, MepackOptions opt)
{
    return call_gesylv2_op_optset_convlog(A, B, S, TT, Q, U, Y, "N", "N", opt);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gesylv2_x0_opt(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X, MepackOptions opt)
{
    return call_gesylv2_x0_op_optset(A, B, S, TT, Q, U, Y, X, "N", "N", opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gesylv2_x0_opt_convlog(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X, MepackOptions opt)
{
    return call_gesylv2_x0_op_optset_convlog(A, B, S, TT, Q, U, Y, X, "N", "N",  opt);
}


template<typename T>
std::tuple<ArrayType<T>>
call_gesylv2_op_opt(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_gesylv2_op_optset(A, B, S, TT, Q, U, Y, opA, opB, opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gesylv2_op_opt_convlog(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_gesylv2_op_optset_convlog(A, B, S, TT, Q, U, Y, opA, opB, opt);
}

template<typename T>
std::tuple<ArrayType<T>>
call_gesylv2_x0_op_opt(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_gesylv2_x0_op_optset(A, B, S, TT, Q, U, Y, X, opA, opB, opt );
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_gesylv2_x0_op_opt_convlog(ArrayType<T>  A,  ArrayType<T> B,ArrayType<T> S, ArrayType<T> TT,  ArrayType<T> Q, ArrayType<T> U, ArrayType<T> Y, ArrayType<T> X,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_gesylv2_x0_op_optset_convlog(A, B, S, TT, Q, U, Y, X, opA, opB, opt);
}








MEXOCT_ENTRY(mepack_sylv2_refine,
R"===(
X = mepack_sylv2_refine(A, B, S, T, Q, U, Y)
X = mepack_sylv2_refine(A, B, S, T, Q, U, Y, X0)
X = mepack_sylv2_refine(A, B, S, T, Q, U, Y, OpA, OpB)
X = mepack_sylv2_refine(A, B, S, T, Q, U, Y, X0, OpA, OpB)
X = mepack_sylv2_refine(A, B, S, T, Q, U, Y, OptSet)
X = mepack_sylv2_refine(A, B, S, T, Q, U, Y, X0, OptSet)
X = mepack_sylv2_refine(A, B, S, T, Q, U, Y, Op, OptSet)
X = mepack_sylv2_refine(A, B, S, T, Q, U, Y, X0, OpA, OpB, OptSet)

Input:
    A      -- coefficient matrix
    B      -- coefficient matrix
    S      -- Schur factor of the coefficient matrix A
    T      -- Schur factor of the coefficient matrix B
    Q      -- Unitary transformation of the Schur decomposition of A
    U      -- Unitary transformation  of the Schur decomposition of B
    Y      -- right hand side
    X0     -- Initial Guess
    OpA    -- Transpose operation on the coefficient matrix A
    OpB    -- Transpose operation on the coefficient matrix B
    OptSet -- Options structure to override the default settings
              of MEPACK.
Return:
    X      -- Solution of the Lyapunov equation

[X, convlog] = mepack_sylv2_refine(A, B, S, T, Q, U, Y)
[X, convlog] = mepack_sylv2_refine(A, B, S, T, Q, U, Y, X0)
[X, convlog] = mepack_sylv2_refine(A, B, S, T, Q, U, Y, OpA, OpB)
[X, convlog] = mepack_sylv2_refine(A, B, S, T, Q, U, Y, X0, OpA, OpB)
[X, convlog] = mepack_sylv2_refine(A, B, S, T, Q, U, Y, OptSet)
[X, convlog] = mepack_sylv2_refine(A, B, S, T, Q, U, Y, X0, OptSet)
[X, convlog] = mepack_sylv2_refine(A, B, S, T, Q, U, Y, Op, OptSet)
[X, convlog] = mepack_sylv2_refine(A, B, S, T, Q, U, Y, X0, OpA, OpB, OptSet)

Input:
    A      -- coefficient matrix
    B      -- coefficient matrix
    S      -- Schur factor of the coefficient matrix A
    T      -- Schur factor of the coefficient matrix B
    Q      -- Unitary transformation of the Schur decomposition of A
    U      -- Unitary transformation  of the Schur decomposition of B
    Y      -- right hand side
    X0     -- Initial Guess
    OpA    -- Transpose operation on the coefficient matrix A
    OpB    -- Transpose operation on the coefficient matrix B
    OptSet -- Options structure to override the default settings
              of MEPACK.
Return:
    X        -- Solution of the Lyapunov equation
    convlog  -- convergence history


The mepack_sylv2_refine function solves the Sylvester equation

   opA(A) * X * opB(B) +/- X = Y,

where A and B are general matrices employing the iterative refinement
strategy. The operation opA (or opB) either returns the matrix or
if Op == 'T' its transpose. In this way, the transposed equation
can be solved without transposing the coefficient matrix before.
The function requires the matrices  A  and B and their Schur decompositions

   [ Q, S ] = schur (A, 'real')
   [ U, T ] = schur (B, 'real')

in order to work. The Schur decomposition can also be obtain from
the output of a previous run to mepack_lyap.
The iterative refinement stops once

  || opA(A) * X_k * opB(B) +/- X_k - Y ||_F
  -----------------------------------------  < sqrt(m*n) * (||A||_F * ||B||_F + 1) * u * tau
             || Y ||_F

is fulfilled. Thereby m is the order of the matrix A, u is the
machine precision in double or single precision, and tau is a
security factor. Tau can be changed in the OptSet.

The function uses the level-3 solver D/SLA_GESYLV_REFINE of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_sylv2_refine:

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

    auto ParamOpAString = Parameter<std::string>("OpA", "Transpose operation on the first coefficient matrix");
    auto ParamOpBString = Parameter<std::string>("OpB", "Transpose operation on the second coefficient matrix");


    auto ParamOptSet   = Parameter<MepackOptions>("optset", "Options for MEPACK");

    auto ParamDoubleA  = Parameter<ArrayType<double>,Traits::SquareMatrix>("A", "double precision coefficient matrix");
    auto ParamDoubleB  = Parameter<ArrayType<double>,Traits::SquareMatrix>("B", "double precision coefficient matrix");
    auto ParamDoubleS  = Parameter<ArrayType<double>,Traits::SquareMatrix>("S", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleT  = Parameter<ArrayType<double>,Traits::SquareMatrix>("T", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleQ  = Parameter<ArrayType<double>,Traits::SquareMatrix>("Q", "double precision unitary matrix corresponding to S");
    auto ParamDoubleU  = Parameter<ArrayType<double>,Traits::SquareMatrix>("U", "double precision unitary matrix corresponding to S");
    auto ParamDoubleY  = Parameter<ArrayType<double>>("Y", "double precision right hand side");
    auto ParamDoubleX0 = Parameter<ArrayType<double>>("X0", "double precision initial guess");

    auto ParamFloatA  = Parameter<ArrayType<float>,Traits::SquareMatrix>("As",  "single precision coefficient matrix");
    auto ParamFloatB  = Parameter<ArrayType<float>,Traits::SquareMatrix>("Bs",  "single precision coefficient matrix");
    auto ParamFloatS  = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ss",  "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatT  = Parameter<ArrayType<float>,Traits::SquareMatrix>("Ts",  "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatQ  = Parameter<ArrayType<float>,Traits::SquareMatrix>("Qs",  "single precision unitary matrix corresponding to Ss");
    auto ParamFloatU  = Parameter<ArrayType<float>,Traits::SquareMatrix>("Us",  "single precision unitary matrix corresponding to Ss");
    auto ParamFloatY  = Parameter<ArrayType<float>>("Ys",  "single precision right hand side");
    auto ParamFloatX0 = Parameter<ArrayType<float>>("X0s", "single precision initial guess");


    auto gesylv2_double                         = makeFunctionExt<1>(call_gesylv2<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY);
    auto gesylv2_double_convlog                 = makeFunctionExt<2>(call_gesylv2_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY);

    auto gesylv2_double_x0                      = makeFunctionExt<1>(call_gesylv2_x0<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamDoubleX0);
    auto gesylv2_double_x0_convlog              = makeFunctionExt<2>(call_gesylv2_x0_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamDoubleX0);

    auto gesylv2_double_op                      = makeFunctionExt<1>(call_gesylv2_op<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamOpAString, ParamOpBString);
    auto gesylv2_double_op_convlog              = makeFunctionExt<2>(call_gesylv2_op_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamOpAString, ParamOpBString);

    auto gesylv2_double_x0_op                   = makeFunctionExt<1>(call_gesylv2_x0_op<double>, "X", ParamDoubleA, ParamDoubleB,  ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU,
                                                                                                    ParamDoubleY, ParamDoubleX0, ParamOpAString, ParamOpBString);
    auto gesylv2_double_x0_op_convlog           = makeFunctionExt<2>(call_gesylv2_x0_op_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT,
                                                                                                                     ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamDoubleX0,
                                                                                                                     ParamOpAString, ParamOpBString);

    auto gesylv2_double_opt                      = makeFunctionExt<1>(call_gesylv2_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamOptSet);
    auto gesylv2_double_opt_convlog              = makeFunctionExt<2>(call_gesylv2_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamOptSet);

    auto gesylv2_double_x0_opt                   = makeFunctionExt<1>(call_gesylv2_x0_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamDoubleX0, ParamOptSet);
    auto gesylv2_double_x0_opt_convlog           = makeFunctionExt<2>(call_gesylv2_x0_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamDoubleX0, ParamOptSet);

    auto gesylv2_double_op_opt                   = makeFunctionExt<1>(call_gesylv2_op_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamOpAString, ParamOpBString, ParamOptSet);
    auto gesylv2_double_op_opt_convlog              = makeFunctionExt<2>(call_gesylv2_op_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamOpAString, ParamOpBString, ParamOptSet);

    auto gesylv2_double_x0_op_opt                   = makeFunctionExt<1>(call_gesylv2_x0_op_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamDoubleX0, ParamOpAString, ParamOpBString, ParamOptSet);
    auto gesylv2_double_x0_op_opt_convlog           = makeFunctionExt<2>(call_gesylv2_x0_op_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleS, ParamDoubleT, ParamDoubleQ, ParamDoubleU, ParamDoubleY, ParamDoubleX0, ParamOpAString, ParamOpBString, ParamOptSet);

    auto gesylv2_float                         = makeFunctionExt<1>(call_gesylv2<float>, "X", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY);
    auto gesylv2_float_convlog                 = makeFunctionExt<2>(call_gesylv2_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY);

    auto gesylv2_float_x0                      = makeFunctionExt<1>(call_gesylv2_x0<float>, "X", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamFloatX0);
    auto gesylv2_float_x0_convlog              = makeFunctionExt<2>(call_gesylv2_x0_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamFloatX0);

    auto gesylv2_float_op                      = makeFunctionExt<1>(call_gesylv2_op<float>, "X", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamOpAString, ParamOpBString);
    auto gesylv2_float_op_convlog              = makeFunctionExt<2>(call_gesylv2_op_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamOpAString, ParamOpBString);

    auto gesylv2_float_x0_op                   = makeFunctionExt<1>(call_gesylv2_x0_op<float>, "X", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamFloatX0, ParamOpAString, ParamOpBString);
    auto gesylv2_float_x0_op_convlog           = makeFunctionExt<2>(call_gesylv2_x0_op_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamFloatX0, ParamOpAString, ParamOpBString);

    auto gesylv2_float_opt                      = makeFunctionExt<1>(call_gesylv2_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamOptSet);
    auto gesylv2_float_opt_convlog              = makeFunctionExt<2>(call_gesylv2_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamOptSet);

    auto gesylv2_float_x0_opt                   = makeFunctionExt<1>(call_gesylv2_x0_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamFloatX0, ParamOptSet);
    auto gesylv2_float_x0_opt_convlog           = makeFunctionExt<2>(call_gesylv2_x0_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamFloatX0, ParamOptSet);

    auto gesylv2_float_op_opt                   = makeFunctionExt<1>(call_gesylv2_op_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamOpAString, ParamOpBString, ParamOptSet);
    auto gesylv2_float_op_opt_convlog              = makeFunctionExt<2>(call_gesylv2_op_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamOpAString, ParamOpBString, ParamOptSet);

    auto gesylv2_float_x0_op_opt                   = makeFunctionExt<1>(call_gesylv2_x0_op_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamFloatX0, ParamOpAString, ParamOpBString, ParamOptSet);
    auto gesylv2_float_x0_op_opt_convlog           = makeFunctionExt<2>(call_gesylv2_x0_op_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatS, ParamFloatT, ParamFloatQ, ParamFloatU, ParamFloatY, ParamFloatX0, ParamOpAString, ParamOpBString, ParamOptSet);

	auto parser = makeArgParser("mepack_sylv2_refine", gesylv2_double,  gesylv2_double_convlog
                                                      , gesylv2_double_x0, gesylv2_double_x0_convlog
                                                      , gesylv2_double_op,  gesylv2_double_op_convlog
                                                      , gesylv2_double_x0_op , gesylv2_double_x0_op_convlog
                                                      , gesylv2_double_opt,  gesylv2_double_opt_convlog
                                                      , gesylv2_double_x0_opt, gesylv2_double_x0_opt_convlog
                                                      , gesylv2_double_op_opt,  gesylv2_double_op_opt_convlog
                                                      , gesylv2_double_x0_op_opt, gesylv2_double_x0_op_opt_convlog
                                                      , gesylv2_float,  gesylv2_float_convlog
                                                      , gesylv2_float_x0, gesylv2_float_x0_convlog
                                                      , gesylv2_float_op,  gesylv2_float_op_convlog
                                                      , gesylv2_float_x0_op, gesylv2_float_x0_op_convlog
                                                      , gesylv2_float_opt,  gesylv2_float_opt_convlog
                                                      , gesylv2_float_x0_opt, gesylv2_float_x0_opt_convlog
                                                      , gesylv2_float_op_opt,  gesylv2_float_op_opt_convlog
                                                      , gesylv2_float_x0_op_opt, gesylv2_float_x0_op_opt_convlog
                                                      );

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


