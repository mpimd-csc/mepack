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
#include "hess_triangular_test.hpp"

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
struct cgsylv_dual_solver {
    constexpr static const char *fn = NULL;
    virtual void solve(const char * FACTA, const char *FACT,  const char *TRANSA, const char *TRANSB, T sgn1, T sgn2,
                       int M, int N, T * A,  int LDA, T * B, int ldb, T * C, int ldc, T * D, int ldd,
                       T * QA, int LDQA, T* ZA, int ldza, T* QB, int ldqb, T * ZB, int ldzb,
                       T *X, int LDX, T *Y, int ldy,  T * SCALE, T *WORK, size_t LDWORK, int *INFO) = 0;
    virtual void blocksize_set(int nb) = 0;
    virtual void blocksize_big_set(int nb) = 0;
    virtual void isolver_set(int is) = 0;
    virtual void fsolver_set(int fs) = 0;
};

template<>
struct cgsylv_dual_solver<double> {
    constexpr static const char *fn = "DLA_GGCSYLV_DUAL";

    static void solve(const char * facta, const char *factb, const char *opA, const char * opB, double sgn1, double sgn2, int M, int N ,
            double * Aptr, int lda, double *Bptr, int ldb, double *Cptr, int ldc, double *Dptr, int ldd,
            double *Qptr, int ldq, double *Zptr, int ldz, double *QBptr, int ldqb, double *ZBptr, int ldzb,
            double *Xptr, int ldx, double *Yptr, int ldy,
            double * scale, double * buffer, size_t ldwork, int *info) {
        mepack_double_ggcsylv_dual(facta, factb, opA, opB, sgn1, sgn2, M, N, Aptr, lda, Bptr, ldb, Cptr, ldc, Dptr, ldd, Qptr, ldq, Zptr, ldz, QBptr, ldqb, ZBptr, ldzb,
                Xptr, ldx, Yptr, ldy, scale, buffer, ldwork, info);
    };
    static void blocksize_mb_set(int nb) {
        mepack_double_tgcsylv_dual_blocksize_mb_set(nb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_double_tgcsylv_dual_blocksize_nb_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_tgcsylv_dual_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_tgcsylv_dual_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_double_tgcsylv_dual_blocksize_2stage_set(nb);
    };

};

template<>
struct cgsylv_dual_solver<float> {
    constexpr static const char *fn = "SLA_GGCSYLV_DUAL";

    static void solve(const char * facta, const char *factb, const char *opA, const char * opB, float sgn1, float sgn2, int M, int N ,
            float * Aptr, int lda, float *Bptr, int ldb, float *Cptr, int ldc, float *Dptr, int ldd,
            float *Qptr, int ldq, float *Zptr, int ldz, float *QBptr, int ldqb, float *ZBptr, int ldzb,
            float *Xptr, int ldx, float *Yptr, int ldy,
            float * scale, float * buffer, size_t ldwork, int *info) {
        mepack_single_ggcsylv_dual(facta, factb, opA, opB, sgn1, sgn2, M, N, Aptr, lda, Bptr, ldb, Cptr, ldc, Dptr, ldd, Qptr, ldq, Zptr, ldz, QBptr, ldqb, ZBptr, ldzb,
                Xptr, ldx, Yptr, ldy, scale, buffer, ldwork, info);
    };
    static void blocksize_mb_set(int nb) {
        mepack_single_tgcsylv_dual_blocksize_mb_set(nb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_single_tgcsylv_dual_blocksize_nb_set(nb);
    };
    static void isolver_set(int isolver) {
        mepack_tgcsylv_dual_isolver_set(isolver);
    };
    static void fsolver_set(int fsolver) {
        mepack_tgcsylv_dual_frontend_solver_set(fsolver);
    };
    static void blocksize_big_set(int nb) {
        mepack_single_tgcsylv_dual_blocksize_2stage_set(nb);
    };

};


template <typename T>
/*                          X,    Y, A,              B           C              D           QA             ZA              QB          ZB */
static std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
GSylvDriver( ArrayType<T>& A, ArrayType<T>& B, ArrayType<T>& C, ArrayType<T>& D, ArrayType<T>&X, ArrayType<T>& Y,
            ArrayType<T>& QA, ArrayType<T>& ZA, ArrayType<T>& QB, ArrayType<T> &ZB,
            std::string opA = std::string("N"), std::string opB = std::string("N"),
        MepackOptions optset = MepackOptions(), bool factorize = false)
{
    if ( X.rows != Y.rows || X.columns != Y.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The dimension of both right hand sides are not the same. size(X) = [ %d  %d ] size(Y) = [ %d %d ]",
                (int) X.rows, (int) X.columns, (int) Y.rows, (int) Y.columns);
    }

    if ( A.rows != Y.rows) {
        mexoct_error("MEPACK:DimMissMatch",
                "The first coefficient matrix and the right hand side need to have the same number of rows size(A) = [ %d %d ] size (Y) = [ %d %d ]",
                (int) A.rows, (int) A.columns, (int) Y.rows, (int) Y.columns);
    }
    if ( C.rows != Y.rows) {
        mexoct_error("MEPACK:DimMissMatch",
                "The third coefficient matrix and the right hand side need to have the same number of rows size(C) = [ %d %d ] size (Y) = [ %d %d ]",
                (int) C.rows, (int) C.columns, (int) Y.rows, (int) Y.columns);
    }

    if ( B.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The second coefficient matrix and the right hand side need to have the same number of columns size(B) = [ %d %d ] size (Y) = [ %d %d ]",
                (int) B.rows, (int) B.columns, (int) Y.rows, (int) Y.columns);
    }

    if ( D.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The fourth coefficient matrix and the right hand side need to have the same number of columns size(D) = [ %d %d ] size (Y) = [ %d %d ]",
                (int) D.rows, (int) D.columns, (int) Y.rows, (int) Y.columns);
    }


    if ( QA.rows != A.rows || QA.columns != A.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The dimensions of the first coefficient matrix and its transformation matrix must fit. size(A) = [ %d %d ] size (QA) = [ %d %d ]",
                (int) A.rows, (int) A.columns, (int) QA.rows, (int) QA.columns);
    }
    if ( ZA.rows != A.rows || ZA.columns != A.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The dimensions of the first coefficient matrix and its transformation matrix must fit. size(A) = [ %d %d ] size (ZA) = [ %d %d ]",
                (int) A.rows, (int) A.columns, (int) ZA.rows, (int) ZA.columns);
    }



    if ( QB.rows != B.rows || QB.columns != B.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The dimensions of the second coefficient matrix and its transformation matrix must fit. size(B) = [ %d %d ] size (QB) = [ %d %d ]",
                (int) B.rows, (int) B.columns, (int) QB.rows, (int) QB.columns);
    }

    if ( ZB.rows != B.rows || ZB.columns != B.columns ) {
        mexoct_error("MEPACK:DimMissMatch",
                "The dimensions of the second coefficient matrix and its transformation matrix must fit. size(B) = [ %d %d ] size (ZB) = [ %d %d ]",
                (int) B.rows, (int) B.columns, (int) ZB.rows, (int) ZB.columns);
    }



    if ( optset.getMB() >= 0 ) {
        cgsylv_dual_solver<T>::blocksize_mb_set(optset.getMB());
    }
    if ( optset.getNB() >= 0 ) {
        cgsylv_dual_solver<T>::blocksize_nb_set(optset.getNB());
    }

    if ( optset.getISolver() >= 0 ) {
        cgsylv_dual_solver<T>::isolver_set(optset.getISolver());
    }

    if ( optset.getFSolver() >= 0 ) {
        cgsylv_dual_solver<T>::fsolver_set(optset.getFSolver());
    }


    if ( optset.getBigNB() >= 0 ) {
        cgsylv_dual_solver<T>::blocksize_big_set(optset.getBigNB());
    }


    ssize_t auxmem;
    std::string factA = std::string("N");
    std::string factB = std::string("N");

    if ( factorize ) {
        /* Check if the matrices are Hessenberg */
        bool ac_hess = mepack_mexoct_check_hess_triangular(A,C);
        bool bd_hess = mepack_mexoct_check_hess_triangular(B,D);

        if ( ac_hess ) factA = std::string("H");
        if ( bd_hess ) factB = std::string("H");
    } else {
        factA = std::string("F");
        factB = std::string("F");
    }

    auxmem = mepack_memory_frontend(cgsylv_dual_solver<T>::fn, factA.c_str(), factB.c_str(), A.rows, B.rows);

    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    T scale = 1;
    T sign1 = 1;
    T sign2 = 1;
    int info = 0;

    if ( optset.getSign() < 0 ) {
        sign1 = -1;
    }
    if ( optset.getSign2() < 0 ) {
        sign2 = -1;
    }


    cgsylv_dual_solver<T>::solve(factA.c_str(), factB.c_str(), opA.c_str(), opB.c_str(), sign1, sign2, A.rows, B.rows, A.ptr(), A.ld, B.ptr(), B.ld, C.ptr(), C.ld, D.ptr(), D.ld,
                QA.ptr(), QA.ld, ZA.ptr(), ZA.ld, QB.ptr(), QB.ld, ZB.ptr(), ZB.ld, X.ptr(), X.ld, Y.ptr(), Y.ld, &scale, auxbuf, auxmem, &info);

    if ( info != 0) {
        mexoct_error("MEPACK:GGSYLV_DUAL", "GGSYLV_DUAL failed with info = %d", (int) info);
    }

    return std::make_tuple(std::move(X), std::move(Y), std::move(A), std::move(B), std::move(C), std::move(D), std::move(QA), std::move(ZA), std::move(QB), std::move(ZB));
}


template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggcsylv_dual_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D,
                     ArrayType<T> X, ArrayType<T> Y, std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    ArrayType<T> QA(A.rows, A.rows);
    ArrayType<T> ZA(A.rows, A.rows);
    ArrayType<T> QB(B.rows, B.rows);
    ArrayType<T> ZB(B.rows, B.rows);

    auto ret = GSylvDriver(A, B, C, D, X, Y, QA, ZA, QB, ZB, opA, opB, optset, true);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret));
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggcsylv_dual(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D, ArrayType<T> X, ArrayType<T> Y)
{
    return call_ggcsylv_dual_op_optset(A, B, C, D, X, Y);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggcsylv_dual_op(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D, ArrayType<T> X, ArrayType<T> Y, std::string opA, std::string opB)
{
    return call_ggcsylv_dual_op_optset(A, B, C, D, X, Y, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>,ArrayType<T>>
call_ggcsylv_dual_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D, ArrayType<T> X, ArrayType<T> Y, MepackOptions optset)
{
   return call_ggcsylv_dual_op_optset(A, B, C, D, X, Y, "N", "N", optset);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_10ret_ggcsylv_dual_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D, ArrayType<T> X, ArrayType<T> Y, std::string opA = std::string("N"), std::string opB = std::string("N"),
                          MepackOptions optset = MepackOptions())
{
    ArrayType<T> QA(A.rows, A.rows);
    ArrayType<T> ZA(A.rows, A.rows);
    ArrayType<T> QB(B.rows, B.rows);
    ArrayType<T> ZB(B.rows, B.rows);
    return GSylvDriver(A, B, C, D, X, Y, QA, ZA, QB, ZB, opA, opB, optset, true);
}



template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_10ret_ggcsylv_dual(ArrayType<T>  A, ArrayType<T> B,ArrayType<T> C, ArrayType<T> D,  ArrayType<T> X, ArrayType<T> Y)
{

    return call_10ret_ggcsylv_dual_op_optset(A, B, C, D, X, Y);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_10ret_ggcsylv_dual_op(ArrayType<T>  A, ArrayType<T> B,ArrayType<T> C, ArrayType<T> D,  ArrayType<T> X, ArrayType<T> Y, std::string opA, std::string opB)
{
    return call_10ret_ggcsylv_dual_op_optset(A, B, C, D, X, Y, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>, ArrayType<T>>
call_10ret_ggcsylv_dual_optset(ArrayType<T>  A, ArrayType<T> B,ArrayType<T> C, ArrayType<T> D,  ArrayType<T> X, ArrayType<T> Y, MepackOptions optset)
{
   return call_10ret_ggcsylv_dual_op_optset(A, B, C, D, X, Y, "N", "N", optset);
}

template <typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggcsylv_dual_prefact_op_optset(ArrayType<T>  As, ArrayType<T> Bs, ArrayType<T> Cs, ArrayType<T> Ds,
        ArrayType<T> QA, ArrayType<T> ZA, ArrayType<T> QB, ArrayType<T> ZB, ArrayType<T> X, ArrayType<T> Y,
        std::string opA = std::string("N"), std::string opB = std::string("N"), MepackOptions optset = MepackOptions())
{
    bool chk = mepack_mexoct_check_qtriangular(As);
    bool chk2 = mepack_mexoct_check_qtriangular(Bs);
    bool chk3 = mepack_mexoct_check_triangular(Cs);
    bool chk4 = mepack_mexoct_check_triangular(Ds);
    if ( !chk ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The first coefficient matrix is not (quasi) upper triangular.");
    }
    if ( !chk2 ) {
        mexoct_error("MEPACK:QuasiUpperTriangular", "The second coefficient matrix is not (quasi) upper triangular.");
    }
    if ( !chk3 ) {
        mexoct_error("MEPACK:UpperTriangular", "The third coefficient matrix is not (quasi) upper triangular.");
    }
    if ( !chk4 ) {
        mexoct_error("MEPACK:UpperTriangular", "The fourth coefficient matrix is not (quasi) upper triangular.");
    }

    auto ret = GSylvDriver(As, Bs, Cs, Ds, X, Y, QA, ZA, QB, ZB, opA, opB, optset, false);
    return std::make_tuple(std::get<0>(ret), std::get<1>(ret));
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggcsylv_dual_prefact(ArrayType<T> As, ArrayType<T> Bs, ArrayType<T> Cs, ArrayType<T> Ds, ArrayType<T> QA, ArrayType<T> ZA, ArrayType<T> QB, ArrayType<T> ZB, ArrayType<T> X, ArrayType<T> Y)
{
    return call_ggcsylv_dual_prefact_op_optset(As, Bs, Cs, Ds, QA, ZA, QB, ZB, X, Y);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggcsylv_dual_prefact_op(ArrayType<T> As, ArrayType<T> Bs, ArrayType<T> Cs, ArrayType<T> Ds, ArrayType<T> QA, ArrayType<T> ZA, ArrayType<T> QB, ArrayType<T> ZB,
        ArrayType<T> X, ArrayType<T> Y, std::string opA, std::string opB)
{
    return call_ggcsylv_dual_prefact_op_optset(As, Bs, Cs, Ds, QA, ZA, QB, ZB, X, Y, opA, opB);
}

template<typename T>
std::tuple<ArrayType<T>, ArrayType<T>>
call_ggcsylv_dual_prefact_optset(ArrayType<T> As, ArrayType<T> Bs, ArrayType<T> Cs, ArrayType<T> Ds, ArrayType<T> QA, ArrayType<T> ZA, ArrayType<T> QB, ArrayType<T> ZB,
        ArrayType<T> X, ArrayType<T> Y, MepackOptions optset)
{
   return call_ggcsylv_dual_prefact_op_optset(As, Bs, Cs, Ds, QA, ZA, QB, ZB, X, Y, "N", "N", optset);
}



MEXOCT_ENTRY(mepack_csylv_dual,
R"===(
[R, L] = mepack_csylv_dual(A, B, C, D, E, F)
[R, L] = mepack_csylv_dual(A, B, C, D, E, F, OpA, OpB)
[R, L] = mepack_csylv_dual(A, B, C, D, E, F, OptSet)
[R, L] = mepack_csylv_dual(A, B, C, D, E, F, OpA, OpB, OptSet)
 Input:
    A      -- first coefficient matrix
    B      -- second coefficient matrix
    C      -- third coefficient matrix
    D      -- fourth coefficient matrix
    E      -- first right hand side
    F      -- second right hand side
    OpA    -- Transpose operation on the left coefficient matrix
    OpB    -- Transpose operation on the right coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.
 Return:
    R      -- first solution of the Sylvester equation
    L      -- second solution of the Sylvester equation

[R, L, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E, F)
[R, L, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E, F, OpA, OpB)
[R, L, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E, F, OptSet)
[R, L, AS, BS, CS, DS, QA, ZA, QB, ZB] = mepack_csylv_dual(A, B, C, D, E, F, OpA, OpB, OptSet)
  Input:
    A      -- first coefficient matrix
    B      -- second coefficient matrix
    C      -- third coefficient matrix
    D      -- fourth coefficient matrix
    E      -- first right hand side
    F      -- second right hand side
    OpA    -- Transpose operation on the left coefficient matrix
    OpB    -- Transpose operation on the right coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

  Return:
    R      -- first solution of the Sylvester equation
    L      -- second solution of the Sylvester equation
    AS, CS -- generalized Schur decomposition of (A, C)
    BS, DS -- generalized Schur decomposition of (B, D)
    QA, ZA -- transformation matrices for (A,C) = (QA*AS*ZA',QA*CS*ZA')
    QB, ZB -- transformation matrices for (B,D) = (QB*BS*ZB',QB*DS*ZB')


[R, L] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E, F)
[R, L] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E, F, OpA, OpB)
[R, L] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E, F, OptSet)
[R, L] = mepack_csylv_dual(AS, BS, CS, DS, QA, ZA, QB, ZB, E, F, OpA, OpB, OptSet)
  Input:
    AS, CS -- Schur form of (A,C)
    BS, DS -- Schur form of (B,D)
    QA, ZA -- transformation matrices for (A,C) = (QA*AS*ZA',QA*CS*ZA')
    QB, ZB -- transformation matrices for (B,D) = (QB*BS*ZB',QB*DS*ZB')
    E      -- first right hand side
    F      -- second right hand side
    OpA    -- Transpose operation on the left coefficient matrix
    OpB    -- Transpose operation on the right coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.
  Return:
    R      -- first solution of the Sylvester equation
    L      -- second solution of the Sylvester equation

The mepack_csylv_dual function solves the coupled Sylvester equation

   opA(A)**T * R + opA(C)**T * L = E,
   +/- R * opB(B)**T +/- L * opB(D)**T = F,

where A, B, C, D are general matrices. The operation opA (or opB) either returns the matrix or
if OpA == 'T' (or opB == 'T') its transpose. In this way, the transposed equation
can be solved without transposing the coefficient matrix before. The function
can also be called with (A, C) and (B, D)  decomposed into their generalized
Schur decompositions, (A,C) = (QA*AS*ZA',QA*CS*ZA') and (QB*BS*ZB',QB*DS*ZB')
If the function is called with nine return values,
the generalized Schur forms of (A, C) and (B,D) including their
unitary transformation matrices are returned.

The function uses the level-3 solver D/SLA_GGCSYLV_DUAL of MEPACK.

The involved generalized Schur decomposition is aware of matrix pairs in
Hessenberg-Triangular form. In this case the initial Hessenberg reduction
of the generalized Schur decomposition is skipped.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_csylv_dual

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
   - sign    -- swap the sign in the first equation
       = 1   -- solve the equation with '+'
       = -1  -- solve the equation with '-'
   - sign2   -- swap the sign in the second equation
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
    auto ParamDoubleC = Parameter<ArrayType<double>,Traits::SquareMatrix>("C", "double precision coefficient matrix");
    auto ParamDoubleD = Parameter<ArrayType<double>,Traits::SquareMatrix>("D", "double precision coefficient matrix");

    auto ParamDoubleAS = Parameter<ArrayType<double>,Traits::SquareMatrix>("AS", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleBS = Parameter<ArrayType<double>,Traits::SquareMatrix>("BS", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleCS = Parameter<ArrayType<double>,Traits::SquareMatrix>("CS", "double precision upper triangular coefficient matrix");
    auto ParamDoubleDS = Parameter<ArrayType<double>,Traits::SquareMatrix>("DS", "double precision upper triangular coefficient matrix");

    auto ParamDoubleQA = Parameter<ArrayType<double>,Traits::SquareMatrix>("QA", "double precision unitary matrix corresponding to A and C");
    auto ParamDoubleZA = Parameter<ArrayType<double>,Traits::SquareMatrix>("ZA", "double precision unitary matrix corresponding to A and C");
    auto ParamDoubleQB = Parameter<ArrayType<double>,Traits::SquareMatrix>("QB", "double precision unitary matrix corresponding to B and D");
    auto ParamDoubleZB = Parameter<ArrayType<double>,Traits::SquareMatrix>("ZB", "double precision unitary matrix corresponding to B and D");

    auto ParamDoubleE = Parameter<ArrayType<double>>("E", "first double precision right hand side");
    auto ParamDoubleF = Parameter<ArrayType<double>>("F", "second double precision right hand side");


    auto ParamFloatA = Parameter<ArrayType<float>,Traits::SquareMatrix>("A", "single precision coefficient matrix");
    auto ParamFloatB = Parameter<ArrayType<float>,Traits::SquareMatrix>("B", "single precision coefficient matrix");
    auto ParamFloatC = Parameter<ArrayType<float>,Traits::SquareMatrix>("C", "single precision coefficient matrix");
    auto ParamFloatD = Parameter<ArrayType<float>,Traits::SquareMatrix>("D", "single precision coefficient matrix");

    auto ParamFloatAS = Parameter<ArrayType<float>,Traits::SquareMatrix>("AS", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatBS = Parameter<ArrayType<float>,Traits::SquareMatrix>("BS", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatCS = Parameter<ArrayType<float>,Traits::SquareMatrix>("CS", "single precision upper triangular coefficient matrix");
    auto ParamFloatDS = Parameter<ArrayType<float>,Traits::SquareMatrix>("DS", "single precision upper triangular coefficient matrix");

    auto ParamFloatQA = Parameter<ArrayType<float>,Traits::SquareMatrix>("QA", "single precision unitary matrix corresponding to A and C");
    auto ParamFloatZA = Parameter<ArrayType<float>,Traits::SquareMatrix>("ZA", "single precision unitary matrix corresponding to A and C");
    auto ParamFloatQB = Parameter<ArrayType<float>,Traits::SquareMatrix>("QB", "single precision unitary matrix corresponding to B and D");
    auto ParamFloatZB = Parameter<ArrayType<float>,Traits::SquareMatrix>("ZB", "single precision unitary matrix corresponding to B and D");

    auto ParamFloatE = Parameter<ArrayType<float>>("E", "first single precision right hand side");
    auto ParamFloatF = Parameter<ArrayType<float>>("F", "second single precision right hand side");


    std::string return_str("[R, L, AS, BS, CS, DS, QA, ZA, QB, ZB]");


    auto cgsylv_dual_double                    = makeFunctionExt<2>(call_ggcsylv_dual<double>, "X",  ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF);
    auto cgsylv_dual_double_optset             = makeFunctionExt<2>(call_ggcsylv_dual_optset<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF, ParamOptSet);
    auto cgsylv_dual_double_op                 = makeFunctionExt<2>(call_ggcsylv_dual_op<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString);
    auto cgsylv_dual_double_op_optset          = makeFunctionExt<2>(call_ggcsylv_dual_op_optset<double>, "X", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString, ParamOptSet);
    auto cgsylv_dual_10ret_double               = makeFunctionExt<10>(call_10ret_ggcsylv_dual<double>, return_str, ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF);
    auto cgsylv_dual_10ret_double_optset        = makeFunctionExt<10>(call_10ret_ggcsylv_dual_optset<double>, return_str, ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF, ParamOptSet);
    auto cgsylv_dual_10ret_double_op            = makeFunctionExt<10>(call_10ret_ggcsylv_dual_op<double>, return_str , ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString);
    auto cgsylv_dual_10ret_double_op_optset     = makeFunctionExt<10>(call_10ret_ggcsylv_dual_op_optset<double>, return_str, ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD,
                                                                                                     ParamDoubleE, ParamDoubleF, ParamOpAString,
                                                                                                                   ParamOpBString, ParamOptSet);
    auto cgsylv_dual_prefact_double            = makeFunctionExt<2>(call_ggcsylv_dual_prefact<double>, "X",           ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS,
                                                                                                         ParamDoubleQA, ParamDoubleZA, ParamDoubleQB, ParamDoubleZB,
                                                                                                         ParamDoubleE, ParamDoubleF);
    auto cgsylv_dual_prefact_double_optset     = makeFunctionExt<2>(call_ggcsylv_dual_prefact_optset<double>, "X",    ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS,
                                                                                                         ParamDoubleQA, ParamDoubleZA, ParamDoubleQB, ParamDoubleZB,
                                                                                                         ParamDoubleE, ParamDoubleF, ParamOptSet);
    auto cgsylv_dual_prefact_double_op         = makeFunctionExt<2>(call_ggcsylv_dual_prefact_op<double>, "X",        ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS,
                                                                                                         ParamDoubleQA, ParamDoubleZA, ParamDoubleQB, ParamDoubleZB,
                                                                                                         ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString);
    auto cgsylv_dual_prefact_double_op_optset  = makeFunctionExt<2>(call_ggcsylv_dual_prefact_op_optset<double>, "X", ParamDoubleAS, ParamDoubleBS, ParamDoubleCS, ParamDoubleDS,
                                                                                                         ParamDoubleQA, ParamDoubleZA, ParamDoubleQB, ParamDoubleZB,
                                                                                                         ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString, ParamOptSet);

    auto cgsylv_dual_float                    = makeFunctionExt<2>(call_ggcsylv_dual<float>, "X",  ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF);
    auto cgsylv_dual_float_optset             = makeFunctionExt<2>(call_ggcsylv_dual_optset<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF, ParamOptSet);
    auto cgsylv_dual_float_op                 = makeFunctionExt<2>(call_ggcsylv_dual_op<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString);
    auto cgsylv_dual_float_op_optset          = makeFunctionExt<2>(call_ggcsylv_dual_op_optset<float>, "X", ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString, ParamOptSet);
    auto cgsylv_dual_10ret_float               = makeFunctionExt<10>(call_10ret_ggcsylv_dual<float>, return_str, ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF);
    auto cgsylv_dual_10ret_float_optset        = makeFunctionExt<10>(call_10ret_ggcsylv_dual_optset<float>, return_str, ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF, ParamOptSet);
    auto cgsylv_dual_10ret_float_op            = makeFunctionExt<10>(call_10ret_ggcsylv_dual_op<float>, return_str , ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString);
    auto cgsylv_dual_10ret_float_op_optset     = makeFunctionExt<10>(call_10ret_ggcsylv_dual_op_optset<float>, return_str, ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF, ParamOpAString,
                                                                                                                   ParamOpBString, ParamOptSet);
    auto cgsylv_dual_prefact_float            = makeFunctionExt<2>(call_ggcsylv_dual_prefact<float>, "X",           ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS,
                                                                                                         ParamFloatQA, ParamFloatZA, ParamFloatQB, ParamFloatZB,
                                                                                                         ParamFloatE, ParamFloatF);
    auto cgsylv_dual_prefact_float_optset     = makeFunctionExt<2>(call_ggcsylv_dual_prefact_optset<float>, "X",    ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS,
                                                                                                         ParamFloatQA, ParamFloatZA, ParamFloatQB, ParamFloatZB,
                                                                                                         ParamFloatE, ParamFloatF, ParamOptSet);
    auto cgsylv_dual_prefact_float_op         = makeFunctionExt<2>(call_ggcsylv_dual_prefact_op<float>, "X",        ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS,
                                                                                                         ParamFloatQA, ParamFloatZA, ParamFloatQB, ParamFloatZB,
                                                                                                         ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString);
    auto cgsylv_dual_prefact_float_op_optset  = makeFunctionExt<2>(call_ggcsylv_dual_prefact_op_optset<float>, "X", ParamFloatAS, ParamFloatBS, ParamFloatCS, ParamFloatDS,
                                                                                                         ParamFloatQA, ParamFloatZA, ParamFloatQB, ParamFloatZB,
                                                                                                         ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString, ParamOptSet);





	auto parser = makeArgParser("mepack_ggcsylv_dual", cgsylv_dual_double,  cgsylv_dual_double_optset, cgsylv_dual_double_op, cgsylv_dual_double_op_optset,
                                                 cgsylv_dual_10ret_double,  cgsylv_dual_10ret_double_optset, cgsylv_dual_10ret_double_op, cgsylv_dual_10ret_double_op_optset,
                                                 cgsylv_dual_prefact_double, cgsylv_dual_prefact_double_optset, cgsylv_dual_prefact_double_op, cgsylv_dual_prefact_double_op_optset,
                                                 cgsylv_dual_10ret_float,  cgsylv_dual_10ret_float_optset, cgsylv_dual_10ret_float_op, cgsylv_dual_10ret_float_op_optset,
                                                 cgsylv_dual_prefact_float, cgsylv_dual_prefact_float_optset, cgsylv_dual_prefact_float_op, cgsylv_dual_prefact_float_op_optset,
                                                 cgsylv_dual_float, cgsylv_dual_float_op, cgsylv_dual_float_optset, cgsylv_dual_float_op_optset);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


