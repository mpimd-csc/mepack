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
struct cgsylv_dual_refine_solver {
    constexpr static const char *fn = NULL;
    virtual void solve (const char * TRANSA, const char * TRANSB, const char *GUESS, T SGN1, T SGN2, int M, int N , T * A, int LDA,
        T *B, int LDB, T * C, int LDC, T * D, int LDD, T *R, int LDR, T* L, int LDL, T* E, int LDE, T * F, int LDF,
        T *AS, int LDAS, T *BS, int LDBS,
        T *CS, int LDCS, T *DS, int LDDS,
        T * Q, int LDQ, T * Z, int LDZ, T *U, int LDU, T *V, int LDV,
        int *MAXIT, T *TAU, T *CONVLOG, T * WORK, size_t LWORK, int *INFO);
    virtual void blocksize_set(int nb) = 0;
    virtual void blocksize_big_set(int nb) = 0;
    virtual void isolver_set(int is) = 0;
    virtual void fsolver_set(int fs) = 0;
};

template<>
struct cgsylv_dual_refine_solver<double> {
    constexpr static const char *fn = "DLA_GGCSYLV_DUAL_REFINE";

    static void solve (const char * TRANSA, const char * TRANSB, const char *GUESS, double SGN1, double SGN2, int M, int N , double * A, int LDA,
        double *B, int LDB, double *C, int LDC, double *D, int LDD,
        double *R, int LDR, double *L, int LDL, double *E, int LDE, double *F, int LDF,
        double *AS, int LDAS, double *BS, int LDBS,
        double *CS, int LDCS, double *DS, int LDDS, double * Q, int LDQ, double *Z, int LDZ, double *U, int LDU, double *V, int LDV,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LWORK, int *INFO)
    {
        mepack_double_ggcsylv_dual_refine(TRANSA, TRANSB, GUESS, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, R, LDR, L, LDL, E, LDE, F, LDF,
                                     AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, Q, LDQ, Z, LDZ, U, LDU, V, LDV,
                                     MAXIT, TAU, CONVLOG, WORK, LWORK, INFO);
    }

    static void blocksize_mb_set(int nb) {
        mepack_double_tgcsylv_blocksize_mb_set(nb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_double_tgcsylv_blocksize_nb_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_tgcsylv_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_tgcsylv_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_double_tgcsylv_blocksize_2stage_set(nb);
    };

};

template<>
struct cgsylv_dual_refine_solver<float> {
    constexpr static const char *fn = "SLA_GGCSYLV_DUAL_REFINE";

    static void solve (const char * TRANSA, const char * TRANSB, const char *GUESS, float SGN1, float SGN2, int M, int N , float * A, int LDA,
        float *B, int LDB, float *C, int LDC, float *D, int LDD,
        float *R, int LDR, float *L, int LDL, float *E, int LDE, float *F, int LDF,
        float *AS, int LDAS, float *BS, int LDBS,
        float *CS, int LDCS, float *DS, int LDDS, float * Q, int LDQ, float *Z, int LDZ, float *U, int LDU, float *V, int LDV,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LWORK, int *INFO)
    {
        mepack_single_ggcsylv_dual_refine(TRANSA, TRANSB, GUESS, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD,
                                    R, LDR, L, LDL,E, LDE, F, LDF, AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, Q, LDQ, Z, LDZ, U, LDU, V, LDV,
                                    MAXIT, TAU, CONVLOG, WORK, LWORK, INFO);
    }

    static void blocksize_mb_set(int nb) {
        mepack_single_tgcsylv_blocksize_mb_set(nb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_single_tgcsylv_blocksize_nb_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_tgcsylv_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_tgsylv_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_single_tgcsylv_blocksize_2stage_set(nb);
    };

};


template <typename T>
static std::tuple<ArrayType<T>, ArrayType<T>,ArrayType<T>>
CGSylvDualDriver(ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
            ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
            ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
            ArrayType<T> E, ArrayType<T> F, ArrayType<T> R, ArrayType<T> L,
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

    if ( E.rows != F.rows || E.columns != F.columns) {
         mexoct_error("MEPACK:DimMissMatch",
                "The right hand sides have different shapes. size(E) = [ %d %d ] size (F) = [ %d %d ]", (int) E.rows, (int) E.columns,
                (int) F.rows, (int) F.columns);
    }

    if ( R.rows != L.rows || R.columns != L.columns) {
         mexoct_error("MEPACK:DimMissMatch",
                "The initial guesses have different shapes. size(R) = [ %d %d ] size (L) = [ %d %d ]", (int) R.rows, (int) R.columns,
                (int) L.rows, (int) L.columns);
    }
    if (A.rows != E.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix A and the right hand side must have the same number of rows size(A) = [ %d %d ] size (E) = [ %d %d ]", (int) A.rows, (int) A.columns,
                (int) E.rows, (int) E.columns);
    }

    if (B.columns != E.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix B and the right hand side must have the same number of columns size(B) = [ %d %d ] size (E) = [ %d %d ]", (int) B.rows, (int) B.columns,
                (int) E.rows, (int) E.columns);
    }

    if (C.rows != E.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix C and the right hand side must have the same number of rows size(C) = [ %d %d ] size (E) = [ %d %d ]", (int) C.rows, (int) C.columns,
                (int) E.rows, (int) E.columns);
    }

    if (D.columns != E.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix D and the right hand side must have the same number of columns size(D) = [ %d %d ] size (E) = [ %d %d ]", (int) D.rows, (int) D.columns,
                (int) E.rows, (int) E.columns);
    }


    if (AS.rows != E.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix AS and the right hand side must have the same number of rows size(AS) = [ %d %d ] size (E) = [ %d %d ]", (int) AS.rows, (int) AS.columns,
                (int) E.rows, (int) E.columns);
    }

    if (BS.columns != E.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix BS and the right hand side must have the same number of columns size(BS) = [ %d %d ] size (E) = [ %d %d ]", (int) BS.rows, (int) BS.columns,
                (int) E.rows, (int) E.columns);
    }

    if (CS.rows != E.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix CS and the right hand side must have the same number of rows size(CS) = [ %d %d ] size (E) = [ %d %d ]", (int) CS.rows, (int) CS.columns,
                (int) E.rows, (int) E.columns);
    }

    if (DS.columns != E.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The coefficient matrix DS and the right hand side must have the same number of columns size(DS) = [ %d %d ] size (E) = [ %d %d ]", (int) DS.rows, (int) DS.columns,
                (int) E.rows, (int) E.columns);
    }

    if (E.rows != R.rows || E.columns != R.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The initial guess matrix and the right hand side have to be of the same dimension. size(R0) = [ %d %d ] size (E) = [ %d %d ]", (int) R.rows, (int) R.columns,
                (int) E.rows, (int) E.columns);
    }

    if (Q.rows != E.rows ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The unitary transformation matrix Q and the right hand side must have the same number of rows. size(Q) = [ %d %d ] size (E) = [ %d %d ]", (int) Q.rows, (int) Q.columns,
                (int) E.rows, (int) E.columns);
    }

    if (U.columns != E.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The unitary transformation matrix U and the right hand side must have the same number of columns. size(U) = [ %d %d ] size (E) = [ %d %d ]", (int) U.rows, (int) U.columns,
                (int) E.rows, (int) E.columns);
    }


    if ( optset.getMB() >= 0 ) {
        cgsylv_dual_refine_solver<T>::blocksize_mb_set(optset.getMB());
    }
    if ( optset.getNB() >= 0 ) {
        cgsylv_dual_refine_solver<T>::blocksize_nb_set(optset.getNB());
    }

    if ( optset.getISolver() >= 0 ) {
        cgsylv_dual_refine_solver<T>::isolver_set(optset.getISolver());
    }

    if ( optset.getBigNB() >= 0 ) {
        cgsylv_dual_refine_solver<T>::blocksize_big_set(optset.getBigNB());
    }

    if ( optset.getFSolver() >= 0 ) {
        cgsylv_dual_refine_solver<T>::fsolver_set(optset.getFSolver());
    }

    if ( optset.getMaxit() < 1) {
        mexoct_error("MEPACK:MaxitInvalid", "The maximum number of iterations is invalid. maxit = %d", (int) optset.getMaxit());
    }

    if ( optset.getTau() < 0 ) {
        mexoct_error("MEPACK:TolInvalid", "The tolerance is invalid. tol = %20.15e", (double) optset.getTau());
    }

    int maxit = (mexoct_int) optset.getMaxit();
    T  tol = (T) optset.getTau();

    T sign1 = 1;
    T sign2 = 1;

    if ( optset.getSign() < 0 ) {
        sign1 = -1;
    }

    if ( optset.getSign2() < 0 ) {
        sign2 = -1;
    }
    ssize_t auxmem;
    auxmem = mepack_memory_frontend(cgsylv_dual_refine_solver<T>::fn, (const char *){"F"}, (const char *){"F"}, A.rows, B.rows);
    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    ArrayType<T> convlog_internal(maxit, 1, false);

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    int info = 0;

    cgsylv_dual_refine_solver<T>::solve(opA.c_str(), opB.c_str(), guess.c_str(), sign1, sign2 , A.rows, B.rows,
                                 A.ptr(), A.ld, B.ptr(), B.ld, C.ptr(), C.ld, D.ptr(), D.ld,
                                 R.ptr(), R.ld, L.ptr(), L.ld, E.ptr(), E.ld, F.ptr(), F.ld,
                                 AS.ptr(), AS.ld, BS.ptr(), BS.ld, CS.ptr(), CS.ld, DS.ptr(), DS.ld,
                                 Q.ptr(), Q.ld, Z.ptr(), Z.ld, U.ptr(), U.ld, V.ptr(), V.ld,
                                 &maxit, &tol, convlog_internal.ptr(), auxbuf, auxmem, &info);

    ArrayType<T> convlog( maxit, 1, false);
    for ( auto i = 0; i < maxit; i++) convlog(i, 0) = convlog_internal(i, 0);

    if ( info != 0) {
        mexoct_error("MEPACK:CGSYLV_DUAL_REFINE", "CGSYLV_DUAL_REFINE failed with info = %d", (int) info);
    }
    return std::make_tuple(R, L, convlog) ;
}


/* Master routines */
template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggcsylv_dual_op_optset(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{

    ArrayType<T> R(A.rows, B.rows);
    ArrayType<T> L(A.rows, B.rows);
    auto ret = CGSylvDualDriver(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L, opA, opB, "N", optset);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret));
}

template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_ggcsylv_dual_op_optset_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    ArrayType<T> R(A.rows, B.rows);
    ArrayType<T> L(A.rows, B.rows);
    auto ret = CGSylvDualDriver(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L, opA, opB, "N", optset);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret), std::get<2>(ret));
}

template <typename T>
std::tuple<ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_x0_op_optset(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F, ArrayType<T> R, ArrayType<T> L,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    auto ret = CGSylvDualDriver(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L, opA, opB, "I", optset);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret));
}

template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_x0_op_optset_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F,  ArrayType<T> R, ArrayType<T> L,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    auto ret = CGSylvDualDriver(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L, opA, opB, "I", optset);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret), std::get<2>(ret));
}



/* Shorter Interfaces */

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F)
{
    return call_ggcsylv_dual_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F)
{
    return call_ggcsylv_dual_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_x0(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F, ArrayType<T> R, ArrayType<T> L)
{
    return call_ggcsylv_dual_x0_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_x0_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F, ArrayType<T> R, ArrayType<T> L)
{
    return call_ggcsylv_dual_x0_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_op(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F,
        std::string opA, std::string opB)
{
    return call_ggcsylv_dual_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_op_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F,
        std::string opA, std::string opB)
{
    return call_ggcsylv_dual_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_x0_op(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F, ArrayType<T> R, ArrayType<T> L,
        std::string opA, std::string opB)
{
    return call_ggcsylv_dual_x0_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_x0_op_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F, ArrayType<T> R, ArrayType<T> L,
        std::string opA, std::string opB)
{
    return call_ggcsylv_dual_x0_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_opt(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F,
        MepackOptions opt)
{
    return call_ggcsylv_dual_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, "N", "N",  opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_opt_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F,
        MepackOptions opt)
{
    return call_ggcsylv_dual_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, "N", "N", opt);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_x0_opt(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F, ArrayType<T> R, ArrayType<T> L,
        MepackOptions opt)
{
    return call_ggcsylv_dual_x0_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L, "N", "N", opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_x0_opt_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F, ArrayType<T> R, ArrayType<T> L,
        MepackOptions opt)
{
    return call_ggcsylv_dual_x0_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V , E, F, R, L, "N", "N",  opt);
}


template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_op_opt(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_ggcsylv_dual_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, opA, opB, opt);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_op_opt_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_ggcsylv_dual_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E,F, opA, opB, opt);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_x0_op_opt(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F, ArrayType<T> R, ArrayType<T> L,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_ggcsylv_dual_x0_op_optset(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L, opA, opB, opt );
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_x0_op_opt_convlog(
        ArrayType<T> A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
        ArrayType<T> AS, ArrayType<T> BS, ArrayType<T> CS, ArrayType<T> DS,
        ArrayType<T> Q, ArrayType<T> Z, ArrayType<T> U, ArrayType<T> V,
        ArrayType<T> E, ArrayType<T> F, ArrayType<T> R, ArrayType<T> L,
        std::string opA, std::string opB, MepackOptions opt)
{
    return call_ggcsylv_dual_x0_op_optset_convlog(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R, L, opA, opB, opt);
}

MEXOCT_ENTRY(mepack_csylv_dual_refine,
R"===(
[R, L] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F)
[R, L] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0)
[R, L] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, OpA, OpB)
[R, L] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OpA, OpB)
[R, L] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, OptSet)
[R, L] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OptSet)
[R, L] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, Op, OptSet)
[R, L] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OpA, OpB, OptSet)
[R, L, convlog] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F)
[R, L, convlog] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0)
[R, L, convlog] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, OpA, OpB)
[R, L, convlog] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OpA, OpB)
[R, L, convlog] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, OptSet)
[R, L, convlog] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OptSet)
[R, L, convlog] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, Op, OptSet)
[R, L, convlog] = mepack_csylv_dual_refine(A, B, C, D, AS, BS, CS, DS, Q, Z, U, V, E, F, R0, L0, OpA, OpB, OptSet)


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
    E      -- first right hand side
    F      -- second right hand side
    R0     -- Initial Guess for the first solution
    L0     -- Initial Guess for the second solution
    OpA    -- Transpose operation on the coefficient matrix A
    OpB    -- Transpose operation on the coefficient matrix B
    OptSet -- Options structure to override the default settings
              of MEPACK.
Return:
    R       -- first solution of the coupled Sylvester equation
    L       -- second solution of the coupled Sylvester equation
    convlog -- convergence history

The mepack_csylv_dual_refine function solves the coupled generalized Sylvester equation

    opA(A)'  * R  + opA(C)' * L = E
    +/- R * opB(B)' +/- L * opB(D)' =  F

which is the dual of


   opA(A) * R +/- L * opB(B) = E,
   opA(C) * R +/- L * opB(D) = F,

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
the output of a previous run to mepack_csylv_dual.
The iterative refinement stops once

  || opA(A)' * R_k + opA(C)' * L_k - E ||_F
  ----------------------------------------- < tol
                 || E ||_F

and

  || +/- R_k * opB(B)' +/- L_k * opB(D)' - F||_F
  ---------------------------------------------- < tol
                 || F ||_F

with

  tol = sqrt ( m * n) * (||A||_F + ||C||_F + ||B||_F + ||D||_F) * u * tau

are fulfilled. Thereby m is the order of the matrix A, n the order
of the matrix B, u is the machine precision in double or single precision,
and tau is a security factor. Tau can be changed in the OptSet.

The function uses the level-3 solver D/SLA_GGCSYLV_DUAL_REFINE of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_csylv_dual_refine:

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
   - sign    -- sign in the first of the two equations
       = 1   -- solve the equation with '+'
       = -1  -- solve the equation with '-'
   - sign2   -- sign in the second of the two equations
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
    auto ParamDoubleE  = Parameter<ArrayType<double>>("E", "double precision first right hand side");
    auto ParamDoubleF  = Parameter<ArrayType<double>>("F", "double precision second right hand side");
    auto ParamDoubleR0 = Parameter<ArrayType<double>>("R0", "double precision initial guess for R");
    auto ParamDoubleL0 = Parameter<ArrayType<double>>("L0", "double precision initial guess for L");

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

    auto ParamFloatE  = Parameter<ArrayType<float>>("E",  "single precision first right hand side");
    auto ParamFloatF  = Parameter<ArrayType<float>>("F",  "single precision second right hand side");
    auto ParamFloatR0 = Parameter<ArrayType<float>>("R0", "single precision initial guess for R");
    auto ParamFloatL0 = Parameter<ArrayType<float>>("L0", "single precision initial guess for L");


    auto cgsylv_dual_double                         = makeFunctionExt<2>(call_ggcsylv_dual<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF);
    auto cgsylv_dual_double_convlog                 = makeFunctionExt<3>(call_ggcsylv_dual_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF);

    auto cgsylv_dual_double_x0                      = makeFunctionExt<2>(call_ggcsylv_dual_x0<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamDoubleR0, ParamDoubleL0);
    auto cgsylv_dual_double_x0_convlog              = makeFunctionExt<3>(call_ggcsylv_dual_x0_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamDoubleR0, ParamDoubleL0);

    auto cgsylv_dual_double_op                      = makeFunctionExt<2>(call_ggcsylv_dual_op<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString);
    auto cgsylv_dual_double_op_convlog              = makeFunctionExt<3>(call_ggcsylv_dual_op_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString);

    auto cgsylv_dual_double_x0_op                   = makeFunctionExt<2>(call_ggcsylv_dual_x0_op<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV,
                                                                                                    ParamDoubleE, ParamDoubleF, ParamDoubleR0, ParamDoubleL0, ParamOpAString, ParamOpBString);
    auto cgsylv_dual_double_x0_op_convlog           = makeFunctionExt<3>(call_ggcsylv_dual_x0_op_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD,
                                                                                                                     ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS,
                                                                                                                     ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV,
                                                                                                                     ParamDoubleE, ParamDoubleF, ParamDoubleR0, ParamDoubleL0,
                                                                                                                     ParamOpAString, ParamOpBString);

    auto cgsylv_dual_double_opt                      = makeFunctionExt<2>(call_ggcsylv_dual_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamOptSet);
    auto cgsylv_dual_double_opt_convlog              = makeFunctionExt<3>(call_ggcsylv_dual_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamOptSet);

    auto cgsylv_dual_double_x0_opt                   = makeFunctionExt<2>(call_ggcsylv_dual_x0_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamDoubleR0, ParamDoubleL0, ParamOptSet);
    auto cgsylv_dual_double_x0_opt_convlog           = makeFunctionExt<3>(call_ggcsylv_dual_x0_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamDoubleR0, ParamDoubleL0, ParamOptSet);

    auto cgsylv_dual_double_op_opt                   = makeFunctionExt<2>(call_ggcsylv_dual_op_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString, ParamOptSet);
    auto cgsylv_dual_double_op_opt_convlog              = makeFunctionExt<3>(call_ggcsylv_dual_op_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString, ParamOptSet);

    auto cgsylv_dual_double_x0_op_opt                   = makeFunctionExt<2>(call_ggcsylv_dual_x0_op_opt<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamDoubleR0, ParamDoubleL0, ParamOpAString, ParamOpBString, ParamOptSet);
    auto cgsylv_dual_double_x0_op_opt_convlog           = makeFunctionExt<3>(call_ggcsylv_dual_x0_op_opt_convlog<double>, "X, convlog", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS, ParamDoubleQ, ParamDoubleZ, ParamDoubleU, ParamDoubleV, ParamDoubleE, ParamDoubleF, ParamDoubleR0, ParamDoubleL0, ParamOpAString, ParamOpBString, ParamOptSet);

    auto cgsylv_dual_float                         = makeFunctionExt<2>(call_ggcsylv_dual<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF);
    auto cgsylv_dual_float_convlog                 = makeFunctionExt<3>(call_ggcsylv_dual_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF);

    auto cgsylv_dual_float_x0                      = makeFunctionExt<2>(call_ggcsylv_dual_x0<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamFloatR0, ParamFloatL0);
    auto cgsylv_dual_float_x0_convlog              = makeFunctionExt<3>(call_ggcsylv_dual_x0_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamFloatR0, ParamFloatL0);

    auto cgsylv_dual_float_op                      = makeFunctionExt<2>(call_ggcsylv_dual_op<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString);
    auto cgsylv_dual_float_op_convlog              = makeFunctionExt<3>(call_ggcsylv_dual_op_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString);

    auto cgsylv_dual_float_x0_op                   = makeFunctionExt<2>(call_ggcsylv_dual_x0_op<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamFloatR0, ParamFloatL0, ParamOpAString, ParamOpBString);
    auto cgsylv_dual_float_x0_op_convlog           = makeFunctionExt<3>(call_ggcsylv_dual_x0_op_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamFloatR0, ParamFloatL0, ParamOpAString, ParamOpBString);

    auto cgsylv_dual_float_opt                      = makeFunctionExt<2>(call_ggcsylv_dual_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamOptSet);
    auto cgsylv_dual_float_opt_convlog              = makeFunctionExt<3>(call_ggcsylv_dual_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamOptSet);

    auto cgsylv_dual_float_x0_opt                   = makeFunctionExt<2>(call_ggcsylv_dual_x0_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamFloatR0, ParamFloatL0, ParamOptSet);
    auto cgsylv_dual_float_x0_opt_convlog           = makeFunctionExt<3>(call_ggcsylv_dual_x0_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamFloatR0, ParamFloatL0, ParamOptSet);

    auto cgsylv_dual_float_op_opt                   = makeFunctionExt<2>(call_ggcsylv_dual_op_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString, ParamOptSet);
    auto cgsylv_dual_float_op_opt_convlog              = makeFunctionExt<3>(call_ggcsylv_dual_op_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString, ParamOptSet);

    auto cgsylv_dual_float_x0_op_opt                   = makeFunctionExt<2>(call_ggcsylv_dual_x0_op_opt<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamFloatR0, ParamFloatL0, ParamOpAString, ParamOpBString, ParamOptSet);
    auto cgsylv_dual_float_x0_op_opt_convlog           = makeFunctionExt<3>(call_ggcsylv_dual_x0_op_opt_convlog<float>, "X, convlog", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS, ParamFloatQ, ParamFloatZ, ParamFloatU, ParamFloatV, ParamFloatE, ParamFloatF, ParamFloatR0, ParamFloatL0, ParamOpAString, ParamOpBString, ParamOptSet);

	auto parser = makeArgParser("mepack_ggcsylv_dual_refine", cgsylv_dual_double,  cgsylv_dual_double_convlog
                                                      , cgsylv_dual_double_x0, cgsylv_dual_double_x0_convlog
                                                      , cgsylv_dual_double_op,  cgsylv_dual_double_op_convlog
                                                      , cgsylv_dual_double_x0_op , cgsylv_dual_double_x0_op_convlog
                                                      , cgsylv_dual_double_opt,  cgsylv_dual_double_opt_convlog
                                                      , cgsylv_dual_double_x0_opt, cgsylv_dual_double_x0_opt_convlog
                                                      , cgsylv_dual_double_op_opt,  cgsylv_dual_double_op_opt_convlog
                                                      , cgsylv_dual_double_x0_op_opt, cgsylv_dual_double_x0_op_opt_convlog
                                                      , cgsylv_dual_float,  cgsylv_dual_float_convlog
                                                      , cgsylv_dual_float_x0, cgsylv_dual_float_x0_convlog
                                                      , cgsylv_dual_float_op,  cgsylv_dual_float_op_convlog
                                                      , cgsylv_dual_float_x0_op, cgsylv_dual_float_x0_op_convlog
                                                      , cgsylv_dual_float_opt,  cgsylv_dual_float_opt_convlog
                                                      , cgsylv_dual_float_x0_opt, cgsylv_dual_float_x0_opt_convlog
                                                      , cgsylv_dual_float_op_opt,  cgsylv_dual_float_op_opt_convlog
                                                      , cgsylv_dual_float_x0_op_opt, cgsylv_dual_float_x0_op_opt_convlog
                                                      );

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


