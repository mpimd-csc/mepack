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

/**
 * @addtogroup cgglyap
 * @{
 */

/**
 \brief Iterative Refinement for the Generalized Lyapunov Equation

 \par Purpose:

 \verbatim
 mepack_double_gglyap_refine solves a generalized Lyapunov equation of the following forms

    A * X * B^T  + B * X * A^T = SCALE * Y                                              (1)

 or

    A^T * X * B  -  B^T * X * A =  SCALE * Y                                            (2)

 where A is a M-by-M matrix using iterative refinement.
 The right hand side Y and the solution X are M-by-M matrices.
 The matrix pencil (A,B) needs to be provided as the original data
 as well as in generalized Schur decomposition since both are required in the
 iterative refinement process.

 \endverbatim

 \remark This function is a wrapper for dla_gglyap_refine.

 \see dla_gglyap_refine

 \param[in] TRANS
 \verbatim
          TRANS is String
          Specifies the form of the system of equations with respect to A:
          == 'N':  Equation (1) is solved
          == 'T':  Equation (2) is solved
 \endverbatim

 \param[in] GUESS
 \verbatim
          GUESS is String
          Specifies whether X  contains an initial guess on input or not.
          = 'I': X contains an initial guess for the solution
          == 'N': No initial guess is provided. X is set to zero.
 \endverbatim


 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrix A.  M >= 0.
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
          B is DOUBLE PRECISION array, dimension (LDB,M)
          The array B contains the original matrix B defining the equation.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,M).
 \endverbatim


 \param[in,out] X
 \verbatim
          X is DOUBLE PRECISION array, dimension (LDX,M)
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
          Y is DOUBLE PRECISION array, dimension (LDY,M)
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
          matrix A.
 \endverbatim


 \param[in] LDAS
 \verbatim
          LDAS is INTEGER
          The leading dimension of the array AS.  LDAS >= max(1,M).
 \endverbatim

 \param[in] BS
 \verbatim
          BS is DOUBLE PRECISION array, dimension (LDBS,M)
          The array BS contains the generalized Schur decomposition of the
          matrix B.
 \endverbatim


 \param[in] LDBS
 \verbatim
          LDBS is INTEGER
          The leading dimension of the array BS.  LDBS >= max(1,M).
 \endverbatim


 \param[in] Q
 \verbatim
          Q is DOUBLE PRECISION array, dimension (LDQ,M)
          The array Q contains the left generalized Schur vectors of
          (A, B) as returned by DGGES3.
 \endverbatim

 \param[in] LDQ
 \verbatim
          LDQ is INTEGER
          The leading dimension of the array Q.  LDQ >= max(1,M).
 \endverbatim

 \param[in] Z
 \verbatim
          Z is DOUBLE PRECISION array, dimension (LDZ,M)
          The array Z contains the right generalized Schur vectors of
          (A,B) as returned by DGGES3.
 \endverbatim

 \param[in] LDZ
 \verbatim
          LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= max(1,M).
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

 \date June 2023
*/
void mepack_double_gglyap_refine(const char * TRANS, const char *GUESS, int M , double * A, int LDA,
        double *B, int LDB, double *X, int LDX, double *Y, int LDY, double *AS, int LDAS, double *BS, int LDBS,
        double * Q, int LDQ, double *Z, int LDZ,
        int *MAXIT, double *TAU, double *CONVLOG, double * WORK, size_t LDWORK, int *INFO)
{
    Int _M  = M;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDY = LDY;
    Int _LDX = LDX;
    Int _LDAS = LDAS;
    Int _LDBS = LDBS;
    Int _LDQ  = LDQ;
    Int _LDZ  = LDZ;
    Int _MAXIT = *MAXIT;
    Int _LDWORK = LDWORK;
    Int _INFO = * INFO;
    FORTRAN_LEN_T c1len = strlen(TRANS);
    FORTRAN_LEN_T c2len = strlen(GUESS);

    FC_GLOBAL_(dla_gglyap_refine, DLA_GGLYAP_REFINE) (TRANS, GUESS, &_M, A, &_LDA, B, &_LDB, X, &_LDX, Y, &_LDY, AS, &_LDAS, BS, &_LDBS,
            Q, &_LDQ, Z, &_LDZ,
            &_MAXIT, TAU, CONVLOG, WORK, &_LDWORK, &_INFO, c1len, c2len);

    *INFO = _INFO;
    *MAXIT = _MAXIT;
    return;
}

/**
 \brief Iterative Refinement for the Generalized Lyapunov Equation

 \par Purpose:

 \verbatim
 mepack_single_gglyap_refine solves a generalized Lyapunov equation of the following forms

    A * X * B^T  + B * X * A^T = SCALE * Y                                              (1)

 or

    A^T * X * B  +  B^T * X * A =  SCALE * Y                                            (2)

 where A is a M-by-M matrix using iterative refinement.
 The right hand side Y and the solution X are M-by-M matrices.
 The matrix pencil (A,B) needs to be provided as the original data
 as well as in generalized Schur decomposition since both are required in the
 iterative refinement process.

 \endverbatim

 \remark This function is a wrapper for dla_gglyap_refine.

 \see dla_gglyap_refine

 \param[in] TRANS
 \verbatim
          TRANS is String
          Specifies the form of the system of equations with respect to A:
          == 'N':  Equation (1) is solved
          == 'T':  Equation (2) is solved
 \endverbatim

 \param[in] GUESS
 \verbatim
          GUESS is String
          Specifies whether X  contains an initial guess on input or not.
          = 'I': X contains an initial guess for the solution
          == 'N': No initial guess is provided. X is set to zero.
 \endverbatim


 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrix A.  M >= 0.
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
          B is SINGLE PRECISION array, dimension (LDB,M)
          The array B contains the original matrix B defining the equation.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,M).
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
          matrix A.
 \endverbatim


 \param[in] LDAS
 \verbatim
          LDAS is INTEGER
          The leading dimension of the array AS.  LDAS >= max(1,M).
 \endverbatim

 \param[in] BS
 \verbatim
          BS is SINGLE PRECISION array, dimension (LDBS,M)
          The array BS contains the generalized Schur decomposition of the
          matrix B.
 \endverbatim


 \param[in] LDBS
 \verbatim
          LDBS is INTEGER
          The leading dimension of the array BS.  LDBS >= max(1,M).
 \endverbatim


 \param[in] Q
 \verbatim
          Q is SINGLE PRECISION array, dimension (LDQ,M)
          The array Q contains the left generalized Schur vectors of
          (A, B) as returned by SGGES3.
 \endverbatim

 \param[in] LDQ
 \verbatim
          LDQ is INTEGER
          The leading dimension of the array Q.  LDQ >= max(1,M).
 \endverbatim

 \param[in] Z
 \verbatim
          Z is SINGLE PRECISION array, dimension (LDZ,M)
          The array Z contains the right generalized Schur vectors of
          (A,B) as returned by SGGES3.
 \endverbatim

 \param[in] LDZ
 \verbatim
          LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= max(1,M).
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

 \date June 2023
*/
void mepack_single_gglyap_refine(const char * TRANS, const char *GUESS, int M , float * A, int LDA,
        float *B, int LDB, float *X, int LDX, float *Y, int LDY, float *AS, int LDAS, float *BS, int LDBS,
        float * Q, int LDQ, float *Z, int LDZ,
        int *MAXIT, float *TAU, float *CONVLOG, float * WORK, size_t LDWORK, int *INFO)
{
    Int _M  = M;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _LDY = LDY;
    Int _LDAS = LDAS;
    Int _LDBS = LDBS;
    Int _LDQ  = LDQ;
    Int _LDZ  = LDZ;
    Int _MAXIT = *MAXIT;
    Int _LDWORK = LDWORK;
    Int _INFO = * INFO;
    FORTRAN_LEN_T c1len = strlen(TRANS);
    FORTRAN_LEN_T c2len = strlen(GUESS);

    FC_GLOBAL_(sla_gglyap_refine, SLA_GGLYAP_REFINE) (TRANS, GUESS, &_M, A, &_LDA, B, &_LDB, X, &_LDX, Y, &_LDY, AS, &_LDAS, BS, &_LDBS,
            Q, &_LDQ, Z, &_LDZ,
            &_MAXIT, TAU, CONVLOG, WORK, &_LDWORK, &_INFO, c1len, c2len);

    *INFO = _INFO;
    *MAXIT = _MAXIT;
    return;
}




/**
 * @}
 */

