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
struct gsylv_refine_solver {
    constexpr static const char *fn = NULL;
    virtual void solve (const char * TRANSA, const char * TRANSB, const char *GUESS, T SGN, int M, int N , T * A, int LDA,
        T *B, int LDB, T * C, int LDC, T * D, int LDD, T *X, int LDX, T *Y, int LDY, T *AS, int LDAS, T *BS, int LDBS,
        T *CS, int LDCS, T *DS, int LDDS,
        T * Q, int LDQ, T * Z, int LDZ, T *U, int LDU, T *V, int LDV,
        int *MAXIT, T *TAU, T *CONVLOG, T * WORK, size_t LWORK, int *INFO);
    virtual void blocksize_set(int nb) = 0;
    virtual void blocksize_big_set(int nb) = 0;
    virtual void isolver_set(int is) = 0;
    virtual void fsolver_set(int fs) = 0;
};

template<>
struct gsylv_refine_solver<double> {
    constexpr static const char *fn = "DLA_GGSYLV_REFINE";

    static void solve (const char * TRANSA, const char * TRANSB, const char *GUESS, double SGN, int M, int N , double * A, int LDA,
        double *B, int LDB, double *C, int LDC, double *D, int LDD, double *X, int LDX, double *Y, int LDY, double *AS, int LDAS, double *BS, int LDBS,
        double *CS, int LDCS, double *DS, int LDDS, double * Q, int LDQ, double *Z, int LDZ, double *U, int LDU, double *V, int LDV,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LWORK, int *INFO)
    {
        mepack_double_ggsylv_refine(TRANSA, TRANSB, GUESS, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, Y, LDY,
                                     AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, Q, LDQ, Z, LDZ, U, LDU, V, LDV,
                                     MAXIT, TAU, CONVLOG, WORK, LWORK, INFO);
    }

    static void blocksize_mb_set(int nb) {
        mepack_double_tgsylv_blocksize_mb_set(nb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_double_tgsylv_blocksize_nb_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_tgsylv_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_tgsylv_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_double_tgsylv_blocksize_2stage_set(nb);
    };

};

template<>
struct gsylv_refine_solver<float> {
    constexpr static const char *fn = "SLA_GGSYLV_REFINE";

    static void solve (const char * TRANSA, const char * TRANSB, const char *GUESS, float SGN, int M, int N , float * A, int LDA,
        float *B, int LDB, float *C, int LDC, float *D, int LDD, float *X, int LDX, float *Y, int LDY, float *AS, int LDAS, float *BS, int LDBS,
        float *CS, int LDCS, float *DS, int LDDS, float * Q, int LDQ, float *Z, int LDZ, float *U, int LDU, float *V, int LDV,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LWORK, int *INFO)
    {
        mepack_single_ggsylv_refine(TRANSA, TRANSB, GUESS, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD,
                                    X, LDX, Y, LDY, AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, Q, LDQ, Z, LDZ, U, LDU, V, LDV,
                                    MAXIT, TAU, CONVLOG, WORK, LWORK, INFO);
    }

    static void blocksize_mb_set(int nb) {
        mepack_single_tgsylv_blocksize_mb_set(nb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_single_tgsylv_blocksize_nb_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_tgsylv_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_tgsylv_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_single_tgsylv_blocksize_2stage_set(nb);
    };

};


template <typename T>
static std::tuple<ArrayType<T>, ArrayType<T>>
GSylvDriver(ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
            ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
            ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
            ArrayType<T> Y, ArrayType<T> X,
            std::string opA = std::string("N"), std::string opB = std::string("N"), std::string guess = std::string("N"),
            MepackOptions optset = MepackOptions())
{

    if ( !mepack_mexoct_check_qtriangular(AS) ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix AS is not (quasi) upper triangular.");
    }
    if ( !mepack_mexoct_check_qtriangular(BS) ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix BS is not (quasi) upper triangular.");
    }
    if ( !mepack_mexoct_check_qtriangular(CS) ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix CS is not (quasi) upper triangular.");
    }
    if ( !mepack_mexoct_check_qtriangular(DS) ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The coefficient matrix DS is not (quasi) upper triangular.");
    }

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

    if (C.rows != Y.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix C and the right hand side must have the same number of rows size(C) = [ %d %d ] size (Y) = [ %d %d ]", (int) C.rows, (int) C.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (D.columns != Y.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix D and the right hand side must have the same number of columns size(D) = [ %d %d ] size (Y) = [ %d %d ]", (int) D.rows, (int) D.columns,
                (int) Y.rows, (int) Y.columns);
    }


    if (AS.rows != Y.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix AS and the right hand side must have the same number of rows size(AS) = [ %d %d ] size (Y) = [ %d %d ]", (int) AS.rows, (int) AS.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (BS.columns != Y.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix BS and the right hand side must have the same number of columns size(BS) = [ %d %d ] size (Y) = [ %d %d ]", (int) BS.rows, (int) BS.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (CS.rows != Y.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix CS and the right hand side must have the same number of rows size(CS) = [ %d %d ] size (Y) = [ %d %d ]", (int) CS.rows, (int) CS.columns,
                (int) Y.rows, (int) Y.columns);
    }

    if (DS.columns != Y.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix DS and the right hand side must have the same number of columns size(DS) = [ %d %d ] size (Y) = [ %d %d ]", (int) DS.rows, (int) DS.columns,
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
        gsylv_refine_solver<T>::blocksize_mb_set(optset.getMB());
    }
    if ( optset.getNB() >= 0 ) {
        gsylv_refine_solver<T>::blocksize_nb_set(optset.getNB());
    }

    if ( optset.getISolver() >= 0 ) {
        gsylv_refine_solver<T>::isolver_set(optset.getISolver());
    }

    if ( optset.getBigNB() >= 0 ) {
        gsylv_refine_solver<T>::blocksize_big_set(optset.getBigNB());
    }

    if ( optset.getFSolver() >= 0 ) {
        gsylv_refine_solver<T>::fsolver_set(optset.getFSolver());
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
    auxmem = mepack_memory_frontend(gsylv_refine_solver<T>::fn, (const char *){"F"}, (const char *){"F"}, A.rows, B.rows);
    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    ArrayType<T> convlog_internal(maxit, 1, false);

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    int info = 0;

    gsylv_refine_solver<T>::solve(opA.c_str(), opB.c_str(), guess.c_str(), sign, A.rows, B.rows,
                                 A.ptr(), A.ld, B.ptr(), B.ld, C.ptr(), C.ld, D.ptr(), D.ld,
                                 X.ptr(), X.ld, Y.ptr(), Y.ld,
                                 AS.ptr(), AS.ld, BS.ptr(), BS.ld, CS.ptr(), CS.ld, DS.ptr(), DS.ld,
                                 Q.ptr(), Q.ld, Z.ptr(), Z.ld, U.ptr(), U.ld, V.ptr(), V.ld,
                                 &maxit, &tol, convlog_internal.ptr(), auxbuf, auxmem, &info);

    ArrayType<T> convlog( maxit, 1, false);
    for ( auto i = 0; i < maxit; i++) convlog(i, 0) = convlog_internal(i, 0);

    if ( info != 0) {
        mexoct_error("MEPACK:GGSYLV_REFINE", "GGSYLV_REFINE failed with info = %d", (int) info);
    }
    return std::make_tuple(X, convlog) ;
}


/* Master routines */
template <typename T>
std::tuple<ArrayType<T>>
call_ggsylv_op_optset(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{

    ArrayType<T> X(A.rows, B.rows);
    auto ret = GSylvDriver(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X, opA, opB, "N", optset);
    return std::make_tuple(std::get<0>(ret));
}

template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggsylv_op_optset_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    ArrayType<T> X(A.rows, B.rows);
    auto ret = GSylvDriver(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X, opA, opB, "N", optset);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret));
}

template <typename T>
std::tuple<ArrayType<T>>
call_ggsylv_x0_op_optset(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y, ArrayType<T> X,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    auto ret = GSylvDriver(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X, opA, opB, "I", optset);
    return std::make_tuple(std::get<0>(ret));
}

template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggsylv_x0_op_optset_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y, ArrayType<T> X,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    auto ret = GSylvDriver(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X, opA, opB, "I", optset);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret));
}



/* Shorter Interfaces */

template<typename T>
std::tuple<ArrayType<T>>
call_ggsylv(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y)
{
    return call_ggsylv_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggsylv_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y)
{
    return call_ggsylv_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y);
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggsylv_x0(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y, ArrayType<T> X)
{
    return call_ggsylv_x0_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggsylv_x0_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y, ArrayType<T> X)
{
    return call_ggsylv_x0_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X);
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggsylv_op(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        std::string opA, std::string opB)
{
    return call_ggsylv_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggsylv_op_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        std::string opA, std::string opB)
{
    return call_ggsylv_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggsylv_x0_op(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        ArrayType<T> X, std::string opA, std::string opB)
{
    return call_ggsylv_x0_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggsylv_x0_op_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        ArrayType<T> X, std::string opA, std::string opB)
{
    return call_ggsylv_x0_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggsylv_opt(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        MepackOptions opt)
{
    return call_ggsylv_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, "N", "N",  opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggsylv_opt_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        MepackOptions opt)
{
    return call_ggsylv_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, "N", "N", opt);
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggsylv_x0_opt(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        ArrayType<T> X, MepackOptions opt)
{
    return call_ggsylv_x0_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X, "N", "N", opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggsylv_x0_opt_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        ArrayType<T> X, MepackOptions opt)
{
    return call_ggsylv_x0_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V , Y, X, "N", "N",  opt);
}


template<typename T>
std::tuple<ArrayType<T>>
call_ggsylv_op_opt(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_ggsylv_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, opA, opB, opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggsylv_op_opt_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_ggsylv_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, opA, opB, opt);
}

template<typename T>
std::tuple<ArrayType<T>>
call_ggsylv_x0_op_opt(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        ArrayType<T> X,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_ggsylv_x0_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X, opA, opB, opt );
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggsylv_x0_op_opt_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> Y,
        ArrayType<T> X,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_ggsylv_x0_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X, opA, opB, opt);
}

MEXOCT_ENTRY(mepack_gsylv_refine,
R"===(
X = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y)
X = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0)
X = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, OpA, OpB)
X = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0, OpA, OpB)
X = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, OptSet)
X = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0, OptSet)
X = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, Op, OptSet)
X = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0, OpA, OpB, OptSet)

Input:
    A      -- coefficient matrix
    B      -- coefficient matrix
    C      -- coefficient matrix
    B      -- coefficient matrix
    AS     -- Schur factor of the matrix pair (A, C)
    BS     -- Schur factor of the matrix pair (B, D)
    CS     -- Schur factor of the matrix pair (A, C)
    DS     -- Schur factor of the matrix pair (B, D)
    Q      -- Unitary transformation of the Schur decomposition of (A,C)
    Z      -- Unitary transformation of the Schur decomposition of (A,C)
    U      -- Unitary transformation  of the Schur decomposition of (B,D)
    V      -- Unitary transformation  of the Schur decomposition of (B,D)
    Y      -- right hand side
    X0     -- Initial Guess
    OpA    -- Transpose operation on the coefficient matrix A
    OpB    -- Transpose operation on the coefficient matrix B
    OptSet -- Options structure to override the default settings
              of MEPACK.
Return:
    X      -- Solution of the Sylvester equation

[X, convlog] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y)
[X, convlog] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0)
[X, convlog] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, OpA, OpB)
[X, convlog] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0, OpA, OpB)
[X, convlog] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, OptSet)
[X, convlog] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0, OptSet)
[X, convlog] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, Op, OptSet)
[X, convlog] = mepack_gsylv_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, Y, X0, OpA, OpB, OptSet)

Input:
    A      -- coefficient matrix
    B      -- coefficient matrix
    C      -- coefficient matrix
    B      -- coefficient matrix
    AS     -- Schur factor of the matrix pair (A, C)
    BS     -- Schur factor of the matrix pair (B, D)
    CS     -- Schur factor of the matrix pair (A, C)
    DS     -- Schur factor of the matrix pair (B, D)
    Q      -- Unitary transformation of the Schur decomposition of (A,C)
    Z      -- Unitary transformation of the Schur decomposition of (A,C)
    U      -- Unitary transformation  of the Schur decomposition of (B,D)
    V      -- Unitary transformation  of the Schur decomposition of (B,D)
    Y      -- right hand side
    X0     -- Initial Guess
    OpA    -- Transpose operation on the coefficient matrix A
    OpB    -- Transpose operation on the coefficient matrix B
    OptSet -- Options structure to override the default settings
              of MEPACK.
Return:
    X       -- Solution of the Sylvester equation
    convlog -- convergence log

The mepack_gsylv_refine function solves the generalized Sylvester equation

   opA(A) * X * opB(B) +/- opA(C) * X * opB(D) = Y,

where A, B, C, and D are general matrices employing the iterative refinement
strategy. The operation opA (or opB) either returns the matrix or
if Op == 'T' its transpose. In this way, the transposed equation
can be solved without transposing the coefficient matrix before.
The function requires the matrix pairs (A,C) and (B,D) and their Schur decompositions

    if (exist('OCTAVE_VERSION', 'builtin'))
        [AS, CS, Q, Z] = qz(A, C);
        [BS, DS, U, V] = qz(B, D);
    else
        [AS, CS, Q, Z] = qz(A, C, 'real');
        [BS, DS, U, V] = qz(B, D, 'real');
    end
    Q = Q';
    U = U';

in order to work. The Schur decomposition can also be obtain from
the output of a previous run to mepack_gsylv.
The iterative refinement stops once

  || opA(A) * X_k * opB(B) +/- opA(C) * X_k * opB(D) - Y ||_F
  -----------------------------------------------------------
                          || Y ||_F

                < sqrt(m * n) * (||A||_F*||B||_F+||C||_F*||D||_F) * u * tau

is fulfilled. Thereby m is the order of the matrix A, n the order
of the matrix B, u is the machine precision in double or single precision,
and tau is a security factor. Tau can be changed in the OptSet.

The function uses the level-3 solver D/SLA_GGSYLV_REFINE of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_gsylv_refine:

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
    auto ParamDoubleC  = Parameter<ArrayType<double>,Traits::SquareMatrix>("C", "double precision coefficient matrix");
    auto ParamDoubleD  = Parameter<ArrayType<double>,Traits::SquareMatrix>("D", "double precision coefficient matrix");
    auto ParamDoubleAS  = Parameter<ArrayType<double>,Traits::SquareMatrix>("AS", "double precision coefficient matrix");
    auto ParamDoubleBS  = Parameter<ArrayType<double>,Traits::SquareMatrix>("BS", "double precision coefficient matrix");
    auto ParamDoubleCS  = Parameter<ArrayType<double>,Traits::SquareMatrix>("CS", "double precision coefficient matrix");
    auto ParamDoubleDS  = Parameter<ArrayType<double>,Traits::SquareMatrix>("DS", "double precision coefficient matrix");
    auto ParamDoubleQ  = Parameter<ArrayType<double>,Traits::SquareMatrix>("Q", "double precision unitary matrix corresponding to (A, C)");
    auto ParamDoubleZ  = Parameter<ArrayType<double>,Traits::SquareMatrix>("Z", "double precision unitary matrix corresponding to (A, C)");
    auto ParamDoubleU  = Parameter<ArrayType<double>,Traits::SquareMatrix>("U", "double precision unitary matrix corresponding to (B, D)");
    auto ParamDoubleV  = Parameter<ArrayType<double>,Traits::SquareMatrix>("V", "double precision unitary matrix corresponding to (B, D)");
    auto ParamDoubleY  = Parameter<ArrayType<double>>("Y", "double precision right hand side");
    auto ParamDoubleX0 = Parameter<ArrayType<double>>("X0", "double precision initial guess");

    auto ParamFloatA  = Parameter<ArrayType<float>,Traits::SquareMatrix>("A", "single precision coefficient matrix");
    auto ParamFloatB  = Parameter<ArrayType<float>,Traits::SquareMatrix>("B", "single precision coefficient matrix");
    auto ParamFloatC  = Parameter<ArrayType<float>,Traits::SquareMatrix>("C", "single precision coefficient matrix");
    auto ParamFloatD  = Parameter<ArrayType<float>,Traits::SquareMatrix>("D", "single precision coefficient matrix");
    auto ParamFloatAS = Parameter<ArrayType<float>,Traits::SquareMatrix>("AS", "single precision coefficient matrix");
    auto ParamFloatBS = Parameter<ArrayType<float>,Traits::SquareMatrix>("BS", "single precision coefficient matrix");
    auto ParamFloatCS = Parameter<ArrayType<float>,Traits::SquareMatrix>("CS", "single precision coefficient matrix");
    auto ParamFloatDS = Parameter<ArrayType<float>,Traits::SquareMatrix>("DS", "single precision coefficient matrix");
    auto ParamFloatQ  = Parameter<ArrayType<float>,Traits::SquareMatrix>("Q", "single precision unitary matrix corresponding to (A, C)");
    auto ParamFloatZ  = Parameter<ArrayType<float>,Traits::SquareMatrix>("Z", "single precision unitary matrix corresponding to (A, C)");
    auto ParamFloatU  = Parameter<ArrayType<float>,Traits::SquareMatrix>("U", "single precision unitary matrix corresponding to (B, D)");
    auto ParamFloatV  = Parameter<ArrayType<float>,Traits::SquareMatrix>("V", "single precision unitary matrix corresponding to (B, D)");
    auto ParamFloatY  = Parameter<ArrayType<float>>("Y", "single precision right hand side");
    auto ParamFloatX0 = Parameter<ArrayType<float>>("X0", "single precision initial guess");

    auto ggsylv_double                         = makeFunctionExt<1>(call_ggsylv<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY);
    auto ggsylv_double_convlog                 = makeFunctionExt<2>(call_ggsylv_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY);

    auto ggsylv_double_x0                      = makeFunctionExt<1>(call_ggsylv_x0<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamDoubleX0);
    auto ggsylv_double_x0_convlog              = makeFunctionExt<2>(call_ggsylv_x0_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamDoubleX0);

    auto ggsylv_double_op                      = makeFunctionExt<1>(call_ggsylv_op<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamOpAString, ParamOpBString);
    auto ggsylv_double_op_convlog              = makeFunctionExt<2>(call_ggsylv_op_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamOpAString, ParamOpBString);

    auto ggsylv_double_x0_op                   = makeFunctionExt<1>(call_ggsylv_x0_op<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV,
                                                                                                    ParamDoubleY, ParamDoubleX0, ParamOpAString, ParamOpBString);
    auto ggsylv_double_x0_op_convlog           = makeFunctionExt<2>(call_ggsylv_x0_op_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD,
                                                                                                                     ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS,
                                                                                                                     ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV,
                                                                                                                     ParamDoubleY, ParamDoubleX0,
                                                                                                                     ParamOpAString, ParamOpBString);

    auto ggsylv_double_opt                      = makeFunctionExt<1>(call_ggsylv_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamOptSet);
    auto ggsylv_double_opt_convlog              = makeFunctionExt<2>(call_ggsylv_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamOptSet);

    auto ggsylv_double_x0_opt                   = makeFunctionExt<1>(call_ggsylv_x0_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamDoubleX0, ParamOptSet);
    auto ggsylv_double_x0_opt_convlog           = makeFunctionExt<2>(call_ggsylv_x0_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamDoubleX0, ParamOptSet);

    auto ggsylv_double_op_opt                   = makeFunctionExt<1>(call_ggsylv_op_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamOpAString, ParamOpBString, ParamOptSet);
    auto ggsylv_double_op_opt_convlog              = makeFunctionExt<2>(call_ggsylv_op_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamOpAString, ParamOpBString, ParamOptSet);

    auto ggsylv_double_x0_op_opt                   = makeFunctionExt<1>(call_ggsylv_x0_op_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamDoubleX0, ParamOpAString, ParamOpBString, ParamOptSet);
    auto ggsylv_double_x0_op_opt_convlog           = makeFunctionExt<2>(call_ggsylv_x0_op_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleY, ParamDoubleX0, ParamOpAString, ParamOpBString, ParamOptSet);

    auto ggsylv_float                         = makeFunctionExt<1>(call_ggsylv<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY);
    auto ggsylv_float_convlog                 = makeFunctionExt<2>(call_ggsylv_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY);

    auto ggsylv_float_x0                      = makeFunctionExt<1>(call_ggsylv_x0<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamFloatX0);
    auto ggsylv_float_x0_convlog              = makeFunctionExt<2>(call_ggsylv_x0_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamFloatX0);

    auto ggsylv_float_op                      = makeFunctionExt<1>(call_ggsylv_op<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamOpAString, ParamOpBString);
    auto ggsylv_float_op_convlog              = makeFunctionExt<2>(call_ggsylv_op_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamOpAString, ParamOpBString);

    auto ggsylv_float_x0_op                   = makeFunctionExt<1>(call_ggsylv_x0_op<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamFloatX0, ParamOpAString, ParamOpBString);
    auto ggsylv_float_x0_op_convlog           = makeFunctionExt<2>(call_ggsylv_x0_op_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamFloatX0, ParamOpAString, ParamOpBString);

    auto ggsylv_float_opt                      = makeFunctionExt<1>(call_ggsylv_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamOptSet);
    auto ggsylv_float_opt_convlog              = makeFunctionExt<2>(call_ggsylv_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamOptSet);

    auto ggsylv_float_x0_opt                   = makeFunctionExt<1>(call_ggsylv_x0_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamFloatX0, ParamOptSet);
    auto ggsylv_float_x0_opt_convlog           = makeFunctionExt<2>(call_ggsylv_x0_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamFloatX0, ParamOptSet);

    auto ggsylv_float_op_opt                   = makeFunctionExt<1>(call_ggsylv_op_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamOpAString, ParamOpBString, ParamOptSet);
    auto ggsylv_float_op_opt_convlog              = makeFunctionExt<2>(call_ggsylv_op_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamOpAString, ParamOpBString, ParamOptSet);

    auto ggsylv_float_x0_op_opt                   = makeFunctionExt<1>(call_ggsylv_x0_op_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamFloatX0, ParamOpAString, ParamOpBString, ParamOptSet);
    auto ggsylv_float_x0_op_opt_convlog           = makeFunctionExt<2>(call_ggsylv_x0_op_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatY, ParamFloatX0, ParamOpAString, ParamOpBString, ParamOptSet);

	auto parser = makeArgParser("mepack_gsylv_refine", ggsylv_double,  ggsylv_double_convlog
                                                      , ggsylv_double_x0, ggsylv_double_x0_convlog
                                                      , ggsylv_double_op,  ggsylv_double_op_convlog
                                                      , ggsylv_double_x0_op , ggsylv_double_x0_op_convlog
                                                      , ggsylv_double_opt,  ggsylv_double_opt_convlog
                                                      , ggsylv_double_x0_opt, ggsylv_double_x0_opt_convlog
                                                      , ggsylv_double_op_opt,  ggsylv_double_op_opt_convlog
                                                      , ggsylv_double_x0_op_opt, ggsylv_double_x0_op_opt_convlog
                                                      , ggsylv_float,  ggsylv_float_convlog
                                                      , ggsylv_float_x0, ggsylv_float_x0_convlog
                                                      , ggsylv_float_op,  ggsylv_float_op_convlog
                                                      , ggsylv_float_x0_op, ggsylv_float_x0_op_convlog
                                                      , ggsylv_float_opt,  ggsylv_float_opt_convlog
                                                      , ggsylv_float_x0_opt, ggsylv_float_x0_opt_convlog
                                                      , ggsylv_float_op_opt,  ggsylv_float_op_opt_convlog
                                                      , ggsylv_float_x0_op_opt, ggsylv_float_x0_op_opt_convlog
                                                      );

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


