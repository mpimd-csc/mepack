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

#include "FCMangle.h"
#include "mepack.h"
#include "mepack_internal.h"
#include <string.h>
#include <stdio.h>
/**
 * @addtogroup cggsylv
 * @{
 */
/**
 \brief Iterative Refinement for the Generalized Sylvester Equations

 \par Purpose:

 \verbatim
 mepack_double_ggsylv_refine solves a generalized Sylvester equation of the following forms

    op1(A) * X * op2(B) + SGN * op1(C) * X * op2(D) = Y                    (1)

 with iterative refinement, Thereby  (A,C) is a M-by-M matrix pencil and
 (B,D) is a N-by-N matrix pencil.
 The right hand side Y and the solution X are M-by-N matrices.
 The matrix pencils (A,C) and (B,D) need to be given in the original form as well
 as in their generalized Schur decomposition since both are required in the
 iterative refinement procedure.
 \endverbatim

 \remark This function is a wrapper for \ref dla_ggsylv_refine

 \see dla_ggsylv_refine

 \param[in] TRANSA
 \verbatim
          TRANSA is String
          Specifies the form of the system of equations with respect to A and C :
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is String
          Specifies the form of the system of equations with respect to B and D:
          == 'N':  op2(B) = B,
          == 'T':  op2(B) = B**T
 \endverbatim

 \param[in] GUESS
 \verbatim
          GUESS is String
          Specifies whether X  contains an initial guess on input or not.
          = 'I': X contains an initial guess for the solution
          == 'N': No initial guess is provided. X is set to zero.
 \endverbatim


 \param[in] SGN
 \verbatim
          SGN is DOUBLE PRECISION, allowed values: +/-1
          Specifies the sign between in the first equation.
 \endverbatim

 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrices A and C.  M >= 0.
 \endverbatim

 \param[in] N
 \verbatim
          N is INTEGER
          The order of the matrices B and D.  N >= 0.
 \endverbatim

 \param[in] A
 \verbatim
          A is DOUBLE PRECISION array, dimension (LDA,M)
          The array A contains the original matrix A defining the equation.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim

 \param[in] B
 \verbatim
          B is DOUBLE PRECISION array, dimension (LDB,N)
          The array B contains the original matrix B defining the equation.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
 \endverbatim

 \param[in] C
 \verbatim
          C is DOUBLE PRECISION array, dimension (LDC,M)
          The array C contains the original matrix C defining the equation.
 \endverbatim

 \param[in] LDC
 \verbatim
          LDC is INTEGER
          The leading dimension of the array C.  LDC >= max(1,M).
 \endverbatim

 \param[in] D
 \verbatim
          D is DOUBLE PRECISION array, dimension (LDD,N)
          The array D contains the original matrix D defining the equation.
 \endverbatim

 \param[in] LDD
 \verbatim
          LDD is INTEGER
          The leading dimension of the array D.  LDD >= max(1,N).
 \endverbatim


 \param[in,out] X
 \verbatim
          X is DOUBLE PRECISION array, dimension (LDX,N)
          On input, the array X contains the initial guess.
          On output, the array X contains the solution X.
 \endverbatim

 \param[in] LDX
 \verbatim
          LDX is INTEGER
          The leading dimension of the array X.  LDX >= max(1,M).
 \endverbatim

 \param[in] Y
 \verbatim
          Y is DOUBLE PRECISION array, dimension (LDY,N)
          On input, the array Y contains the right hand side.
 \endverbatim

 \param[in] LDY
 \verbatim
          LDY is INTEGER
          The leading dimension of the array Y.  LDY >= max(1,M).
 \endverbatim

 \param[in] AS
 \verbatim
          AS is DOUBLE PRECISION array, dimension (LDAS,M)
          The array AS contains the generalized Schur decomposition of the
          A.
 \endverbatim


 \param[in] LDAS
 \verbatim
          LDAS is INTEGER
          The leading dimension of the array AS.  LDAS >= max(1,M).
 \endverbatim


 \param[in] BS
 \verbatim
          BS is DOUBLE PRECISION array, dimension (LDBS,N)
          The array BS contains the generalized Schur decomposition of the
          B.
 \endverbatim


 \param[in] LDBS
 \verbatim
          LDBS is INTEGER
          The leading dimension of the array BS.  LDBS >= max(1,N).
 \endverbatim


 \param[in] CS
 \verbatim
          CS is DOUBLE PRECISION array, dimension (LDCS,M)
          The array CS contains the generalized Schur decomposition of the
          C.
 \endverbatim


 \param[in] LDCS
 \verbatim
          LDCS is INTEGER
          The leading dimension of the array CS.  LDCS >= max(1,M).
 \endverbatim


 \param[in] DS
 \verbatim
          DS is DOUBLE PRECISION array, dimension (LDDS,N)
          The array DS contains the generalized Schur decomposition of the
          D.
 \endverbatim


 \param[in] LDDS
 \verbatim
          LDDS is INTEGER
          The leading dimension of the array DS.  LDDS >= max(1,N).
 \endverbatim


 \param[in,out] Q
 \verbatim
          Q is DOUBLE PRECISION array, dimension (LDQ,M)
          The array Q contains the left generalized Schur vectors for (A,C) as returned by DGGES.
 \endverbatim

 \param[in] LDQ
 \verbatim
          LDQ is INTEGER
          The leading dimension of the array Q.  LDQ >= max(1,M).
 \endverbatim

 \param[in] Z
 \verbatim
          Z is DOUBLE PRECISION array, dimension (LDZ,M)
          The array Z contains the right generalized Schur vectors for (A,C) as returned by DGGES.
 \endverbatim

 \param[in] LDZ
 \verbatim
          LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= max(1,M).
 \endverbatim


 \param[in,out] U
 \verbatim
          U is DOUBLE PRECISION array, dimension (LDU,N)
          The array U contains the left generalized Schur vectors for (B,D) as returned by DGGES.
 \endverbatim

 \param[in] LDU
 \verbatim
          LDU is INTEGER
          The leading dimension of the array U.  LDU >= max(1,N).
 \endverbatim

 \param[in] V
 \verbatim
          V is DOUBLE PRECISION array, dimension (LDV,N)
          The array V contains the right generalized Schur vectors for (B,D) as returned by DGGES.
 \endverbatim

 \param[in] LDV
 \verbatim
          LDV is INTEGER
          The leading dimension of the array V.  LDV >= max(1,N).
 \endverbatim

 \param[in,out] MAXIT
 \verbatim
          MAXIT is INTEGER
          On input, MAXIT contains the maximum number of iteration that are performed, MAXIT <= 100
          On exit, MAXIT contains the number of iteration steps taken by the algorithm.
 \endverbatim

 \param[in,out] TAU
 \verbatim
          TAU is DOUBLE PRECISION
          On input, TAU contains the additional security factor for the stopping criterion, typical values are 0.1
          On exit, TAU contains the last relative residual when the stopping criterion got valid.
 \endverbatim

 \param[in,out] CONVLOG
 \verbatim
          CONVLOG is DOUBLE PRECISION array, dimension (MAXIT)
          The CONVLOG array contains the convergence history of the iterative refinement. CONVLOG(I) contains the maximum
          relative residual before it is solved for the I-Th time.
 \endverbatim


 \param[in,out] WORK
 \verbatim
          WORK is DOUBLE PRECISION array, dimension (MAX(1,LDWORK))
          Workspace for the algorithm.
 \endverbatim

 \param[in] LDWORK
 \verbatim
          LDWORK is INTEGER
          Size of the workspace for the algorithm. This can be determined by a call \ref mepack_memory_frontend.
 \endverbatim

 \param[out] INFO
 \verbatim
          INFO is INTEGER
          == 0:  Success
          > 0:  Iteration failed in step INFO
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date October 2023
*/
void mepack_double_ggsylv_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, double SGN, int M, int N , double * A, int LDA,
        double *B, int LDB, double *C, int LDC, double *D, int LDD, double *X, int LDX, double *Y, int LDY,
        double *AS, int LDAS, double *BS, int LDBS, double *CS, int LDCS, double *DS, int LDDS,
        double * Q, int LDQ, double *Z, int LDZ, double *U, int LDU, double *V, int LDV,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LDWORK, int *INFO)
{
    Int _M  = M;
    Int _N  = N;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDC = LDC;
    Int _LDD = LDD;
    Int _LDY = LDY;
    Int _LDX = LDX;
    Int _LDAS = LDAS;
    Int _LDBS = LDBS;
    Int _LDCS = LDCS;
    Int _LDDS = LDDS;
    Int _LDQ  = LDQ;
    Int _LDZ  = LDZ;
    Int _LDU  = LDU;
    Int _LDV  = LDV;
    Int _MAXIT = *MAXIT;
    Int _LDWORK = LDWORK;
    Int _INFO = * INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);
    FORTRAN_LEN_T c3len = strlen(GUESS);

    FC_GLOBAL_(dla_ggsylv_refine, DLA_GGSYLV_REFINE) (TRANSA, TRANSB, GUESS, &SGN, &_M, &_N, A, &_LDA, B, &_LDB, C, &_LDC, D, &_LDD,
            X, &_LDX, Y, &_LDY, AS, &_LDAS, BS, &_LDBS, CS, &_LDCS, DS, &_LDDS,
            Q, &_LDQ, Z, &_LDZ, U, &_LDU, V, &_LDV,
            &_MAXIT, TAU, CONVLOG, WORK, &_LDWORK, &_INFO, c1len, c2len, c3len);

    *INFO = _INFO;
    *MAXIT = _MAXIT;
    return;
}

/**
 \brief Iterative Refinement for the Generalized Sylvester Equations

 \par Purpose:

 \verbatim
 mepack_single_gesylv_refine solves a coupled generalized Sylvester equation of the following forms

    op1(A) * X * op2(B) + SGN * op1(C) * X * op2(D) = Y                    (1)

 with iterative refinement, Thereby  (A,C) is a M-by-M matrix pencil and
 (B,D) is a N-by-N matrix pencil.
 The right hand side Y and the solution X are M-by-N matrices.
 The matrix pencils (A,C) and (B,D) need to be given in the original form as well
 as in their generalized Schur decomposition since both are required in the
 iterative refinement procedure.
 \endverbatim

 \remark This function is a wrapper for \ref sla_ggsylv_refine

 \see sla_ggsylv_refine

 \param[in] TRANSA
 \verbatim
          TRANSA is String
          Specifies the form of the system of equations with respect to A and C :
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is String
          Specifies the form of the system of equations with respect to B and D:
          == 'N':  op2(B) = B,
          == 'T':  op2(B) = B**T
 \endverbatim

 \param[in] GUESS
 \verbatim
          GUESS is String
          Specifies whether X  contains an initial guess on input or not.
          = 'I': X contains an initial guess for the solution
          == 'N': No initial guess is provided. X is set to zero.
 \endverbatim


 \param[in] SGN
 \verbatim
          SGN is SINGLE PRECISION, allowed values: +/-1
          Specifies the sign between in the first equation.
 \endverbatim

 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrices A and C.  M >= 0.
 \endverbatim

 \param[in] N
 \verbatim
          N is INTEGER
          The order of the matrices B and D.  N >= 0.
 \endverbatim

 \param[in] A
 \verbatim
          A is SINGLE PRECISION array, dimension (LDA,M)
          The array A contains the original matrix A defining the equation.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim

 \param[in] B
 \verbatim
          B is SINGLE PRECISION array, dimension (LDB,N)
          The array B contains the original matrix B defining the equation.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
 \endverbatim

 \param[in] C
 \verbatim
          C is SINGLE PRECISION array, dimension (LDC,M)
          The array C contains the original matrix C defining the equation.
 \endverbatim

 \param[in] LDC
 \verbatim
          LDC is INTEGER
          The leading dimension of the array C.  LDC >= max(1,M).
 \endverbatim

 \param[in] D
 \verbatim
          D is SINGLE PRECISION array, dimension (LDA,N)
          The array D contains the original matrix D defining the equation.
 \endverbatim

 \param[in] LDD
 \verbatim
          LDD is INTEGER
          The leading dimension of the array D.  LDD >= max(1,N).
 \endverbatim


 \param[in,out] X
 \verbatim
          X is SINGLE PRECISION array, dimension (LDX,M)
          On input, the array X contains the initial guess.
          On output, the array X contains the solution X.
 \endverbatim

 \param[in] LDX
 \verbatim
          LDX is INTEGER
          The leading dimension of the array X.  LDX >= max(1,M).
 \endverbatim

 \param[in] Y
 \verbatim
          Y is SINGLE PRECISION array, dimension (LDY,M)
          On input, the array Y contains the right hand side.
 \endverbatim

 \param[in] LDY
 \verbatim
          LDY is INTEGER
          The leading dimension of the array Y.  LDY >= max(1,M).
 \endverbatim

 \param[in] AS
 \verbatim
          AS is SINGLE PRECISION array, dimension (LDAS,M)
          The array AS contains the generalized Schur decomposition of the
          A.
 \endverbatim


 \param[in] LDAS
 \verbatim
          LDAS is INTEGER
          The leading dimension of the array AS.  LDAS >= max(1,M).
 \endverbatim


 \param[in] BS
 \verbatim
          BS is SINGLE PRECISION array, dimension (LDBS,N)
          The array AS contains the generalized Schur decomposition of the
          B.
 \endverbatim


 \param[in] LDBS
 \verbatim
          LDBS is INTEGER
          The leading dimension of the array BS.  LDBS >= max(1,N).
 \endverbatim


 \param[in] CS
 \verbatim
          CS is SINGLE PRECISION array, dimension (LDCS,M)
          The array CS contains the generalized Schur decomposition of the
          C.
 \endverbatim


 \param[in] LDCS
 \verbatim
          LDCS is INTEGER
          The leading dimension of the array CS.  LDCS >= max(1,M).
 \endverbatim


 \param[in] DS
 \verbatim
          DS is SINGLE PRECISION array, dimension (LDDS,N)
          The array DS contains the generalized Schur decomposition of the
          D.
 \endverbatim


 \param[in] LDDS
 \verbatim
          LDDS is INTEGER
          The leading dimension of the array DS.  LDAS >= max(1,N).
 \endverbatim


 \param[in,out] Q
 \verbatim
          Q is SINGLE PRECISION array, dimension (LDQ,M)
          The array Q contains the left generalized Schur vectors for (A,C) as returned by DGGES.
 \endverbatim

 \param[in] LDQ
 \verbatim
          LDQ is INTEGER
          The leading dimension of the array Q.  LDQ >= max(1,M).
 \endverbatim

 \param[in] Z
 \verbatim
          Z is SINGLE PRECISION array, dimension (LDZ,M)
          The array Z contains the right generalized Schur vectors for (A,C) as returned by DGGES.
 \endverbatim

 \param[in] LDZ
 \verbatim
          LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= max(1,M).
 \endverbatim


 \param[in,out] U
 \verbatim
          U is SINGLE PRECISION array, dimension (LDU,N)
          The array U contains the left generalized Schur vectors for (B,D) as returned by DGGES.
 \endverbatim

 \param[in] LDU
 \verbatim
          LDU is INTEGER
          The leading dimension of the array U.  LDU >= max(1,N).
 \endverbatim

 \param[in] V
 \verbatim
          V is SINGLE PRECISION array, dimension (LDV,N)
          The array V contains the right generalized Schur vectors for (B,D) as returned by DGGES.
 \endverbatim

 \param[in] LDV
 \verbatim
          LDV is INTEGER
          The leading dimension of the array V.  LDV >= max(1,N).
 \endverbatim

 \param[in,out] MAXIT
 \verbatim
          MAXIT is INTEGER
          On input, MAXIT contains the maximum number of iteration that are performed, MAXIT <= 100
          On exit, MAXIT contains the number of iteration steps taken by the algorithm.
 \endverbatim

 \param[in,out] TAU
 \verbatim
          TAU is SINGLE PRECISION
          On input, TAU contains the additional security factor for the stopping criterion, typical values are 0.1
          On exit, TAU contains the last relative residual when the stopping criterion got valid.
 \endverbatim

 \param[in,out] CONVLOG
 \verbatim
          CONVLOG is SINGLE PRECISION array, dimension (MAXIT)
          The CONVLOG array contains the convergence history of the iterative refinement. CONVLOG(I) contains the maximum
          relative residual before it is solved for the I-Th time.
 \endverbatim


 \param[in,out] WORK
 \verbatim
          WORK is SINGLE PRECISION array, dimension (MAX(1,LDWORK))
          Workspace for the algorithm.
 \endverbatim

 \param[in] LDWORK
 \verbatim
          LDWORK is INTEGER
          Size of the workspace for the algorithm. This can be determined by a call \ref mepack_memory_frontend.
 \endverbatim

 \param[out] INFO
 \verbatim
          INFO is INTEGER
          == 0:  Success
          > 0:  Iteration failed in step INFO
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date October 2023
*/
void mepack_single_ggsylv_refine(const char * TRANSA, const char * TRANSB, const char *GUESS, float SGN, int M, int N , float * A, int LDA,
        float *B, int LDB, float *C, int LDC, float *D, int LDD, float *X, int LDX, float *Y, int LDY,
        float *AS, int LDAS, float *BS, int LDBS, float *CS, int LDCS, float *DS, int LDDS,
        float * Q, int LDQ, float *Z, int LDZ, float *U, int LDU, float *V, int LDV,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LDWORK, int *INFO)
{
    Int _M  = M;
    Int _N  = N;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDC = LDC;
    Int _LDD = LDD;
    Int _LDY = LDY;
    Int _LDX = LDX;
    Int _LDAS = LDAS;
    Int _LDBS = LDBS;
    Int _LDCS = LDCS;
    Int _LDDS = LDDS;
    Int _LDQ  = LDQ;
    Int _LDZ  = LDZ;
    Int _LDU  = LDU;
    Int _LDV  = LDV;
    Int _MAXIT = *MAXIT;
    Int _LDWORK = LDWORK;
    Int _INFO = * INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);
    FORTRAN_LEN_T c3len = strlen(GUESS);

    FC_GLOBAL_(sla_ggsylv_refine, SLA_GGSYLV_REFINE) (TRANSA, TRANSB, GUESS, &SGN, &_M, &_N, A, &_LDA, B, &_LDB, C, &_LDC, D, &_LDD,
            X, &_LDX, Y, &_LDY, AS, &_LDAS, BS, &_LDBS, CS, &_LDCS, DS, &_LDDS,
            Q, &_LDQ, Z, &_LDZ, U, &_LDU, V, &_LDV,
            &_MAXIT, TAU, CONVLOG, WORK, &_LDWORK, &_INFO, c1len, c2len, c3len);

    *INFO = _INFO;
    *MAXIT = _MAXIT;
    return;
}


/**
 * @}
 */

