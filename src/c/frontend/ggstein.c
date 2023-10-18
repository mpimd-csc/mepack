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
 \brief Frontend for the solution of Generalized Stein Equations.

 \par Purpose:

 \verbatim
 mepack_double_ggstein solves a generalized Stein equation of the following forms

    A * X * A^T - B * X * B^T = SCALE * Y                                              (1)

 or

    A^T * X * A - B^T * X * B =  SCALE * Y                                             (2)

 where (A,B) is a M-by-M matrix pencil. The right hand side Y and the solution X are
 M-by-M matrices.  The matrix pencil (A,B) is either in general form, in generalized
 Hessenberg form, or in generalized Schur form where Q and Z also need to be provided.

 \endverbatim

 \remark This function is a wrapper for \ref dla_ggstein.

 \see dla_ggstein

 \param[in] FACT
 \verbatim
          FACT is String
          Specifies how the matrix A is given.
          == 'N':  The matrix pencil (A,B) is given as a general matrix pencil and its Schur decomposition
                  A = Q*S*Z**T, B = Q*R*Z**T will be computed.
          == 'F':  The matrix A is given as its Schur decomposition in terms of S and Q
                  form A = Q*S*Q**T
          == 'H': The matrix pencil (A,B) is given in generalized Hessenberg form and its Schur decomposition
                  A = Q*S*Z**T, B = Q*R*Z**T will be computed.
 \endverbatim

 \param[in] TRANS
 \verbatim
          TRANS is String
          Specifies the form of the system of equations with respect to A:
          == 'N':  Equation (1) is solved.
          == 'T':  Equation (2) is solved.
 \endverbatim

 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrices A, B, Y and X.  M >= 0.
 \endverbatim

 \param[in,out] A
 \verbatim
          A is DOUBLE PRECISION array, dimension (LDA,M)
          If FACT == "N", the matrix A is a general matrix and it is overwritten with the
          quasi upper triangular matrix S of the generalized schur decomposition.
          If FACT == "F", the matrix A contains its (quasi-) upper triangular matrix S of the
          generalized  Schur decomposition of (A,B).
          If FACT == "H", the matrix A is an upper Hessenberg matrix of the generalized
          Hessenberg form (A,B) and it is overwritten with the quasi upper triangular matrix S
          of the generalized Schur decomposition.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim

 \param[in,out] B
 \verbatim
          B is DOUBLE PRECISION array, dimension (LDB,M)
          If FACT == "N", the matrix B a general matrix and it is overwritten with the upper triangular
          matrix of the generalized Schur decomposition.
          If FACT == "F", the matrix B contains its upper triangular matrix R of the generalized schur
          Schur decomposition of (A,B).
          If FACT == "H", the matrix B is the upper triangular matrix of the generalized Hessenberg form
          (A,B) and it is overwritten with the upper triangular matrix of the generalized Schur decomposition.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,M).
 \endverbatim

 \param[in,out] Q
 \verbatim
          Q is DOUBLE PRECISION array, dimension (LDQ,M)
          If FACT == "N", the matrix Q is an empty M-by-M matrix on input and contains the
          left Schur vectors of (A,B) on output.
          If FACT == "F", the matrix Q contains the left Schur vectors of (A,B).
          If FACT == "H", the matrix Q is an empty M-by-M matrix on input and contains the
          left Schur vectors of (A,B) on output.
 \endverbatim

 \param[in] LDQ
 \verbatim
          LDQ is INTEGER
          The leading dimension of the array Q.  LDQ >= max(1,M).
 \endverbatim


 \param[in,out] Z
 \verbatim
          Z is DOUBLE PRECISION array, dimension (LDZ,M)
          If FACT == "N", the matrix Z is an empty M-by-M matrix on input and contains the
          right Schur vectors of (A,B) on output.
          If FACT == "F", the matrix Z contains the right Schur vectors of (A,B).
          If FACT == "H", the matrix Z is an empty M-by-M matrix on input and contains the
          right Schur vectors of (A,B) on output.
 \endverbatim

 \param[in] LDZ
 \verbatim
          LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= max(1,M).
 \endverbatim

 \param[in,out] X
 \verbatim
          X is DOUBLE PRECISION array, dimension (LDX,N)
          On input, the matrix X contains the right hand side Y.
          On output, the matrix X contains the solution of Equation (1) or (2)
          Right hand side Y and the solution X are symmetric M-by-M matrices.
 \endverbatim

 \param[in] LDX
 \verbatim
          LDX is INTEGER
          The leading dimension of the array X.  LDX >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
 \endverbatim

 \param[in,out] WORK
 \verbatim
          WORK is DOUBLE PRECISION array, dimension (MAX(1,LDWORK))
          Workspace for the algorithm. The optimal workspace is given either by \ref mepack_memory_frontend.
 \endverbatim

 \param[in] LDWORK
 \verbatim
          LDWORK is INTEGER
          Size of the workspace for the algorithm counted in floating point numbers of the actual precision.
          The C interface does not support the workspace query by setting LDWORK == -1 on input. In this case,
          the \ref mepack_memory_frontend function have to be used.
 \endverbatim


 \param[out] INFO
 \verbatim
          INFO is INTEGER
          == 0:  successful exit
          = 1:  DGGES failed
          = 2:  DLA_SORT_GEV failed
          = 3:  Inner solver failed
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date October 2023
*/

void mepack_double_ggstein(const char * FACT, const char *TRANS, int M, double * A,  int LDA, double *B, int LDB,
        double * Q, int LDQ, double *Z, int LDZ, double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int *INFO)
{
    Int _M = M;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDQ = LDQ;
    Int _LDZ = LDZ;
    Int _LDX = LDX;
    Int _INFO = * INFO;
    Int _LDWORK = LDWORK;
    FORTRAN_LEN_T c1len = strlen(FACT);
    FORTRAN_LEN_T c2len = strlen(TRANS);

    FC_GLOBAL_(dla_ggstein,DLA_GGSTEIN)(FACT, TRANS, &_M, A, &_LDA, B, &_LDB, Q, &_LDQ, Z, &_LDZ, X, &_LDX, SCALE, WORK, &_LDWORK, &_INFO, c1len, c2len);

    *INFO = _INFO;
    return;

}

/**
 \brief Frontend for the solution of Generalized Stein Equations.

 \par Purpose:

 \verbatim
 mepack_single_ggstein solves a generalized Stein equation of the following forms

    A * X * A^T - B * X * B^T = SCALE * Y                                              (1)

 or

    A^T * X * A - B^T * X * B =  SCALE * Y                                             (2)

 where (A,B) is a M-by-M matrix pencil. The right hand side Y and the solution X are
 M-by-M matrices.  The matrix pencil (A,B) is either in general form, in generalized
 Hessenberg form, or in generalized Schur form where Q and Z also need to be provided.

 \endverbatim

 \remark This function is a wrapper for \ref sla_ggstein.

 \see sla_ggstein

 \param[in] FACT
 \verbatim
          FACT is String
          Specifies how the matrix A is given.
          == 'N':  The matrix pencil (A,B) is given as a general matrix pencil and its Schur decomposition
                  A = Q*S*Z**T, B = Q*R*Z**T will be computed.
          == 'F':  The matrix A is given as its Schur decomposition in terms of S and Q
                  form A = Q*S*Q**T
          == 'H': The matrix pencil (A,B) is given in generalized Hessenberg form and its Schur decomposition
                  A = Q*S*Z**T, B = Q*R*Z**T will be computed.
 \endverbatim

 \param[in] TRANS
 \verbatim
          TRANS is String
          Specifies the form of the system of equations with respect to A:
          == 'N':  Equation (1) is solved.
          == 'T':  Equation (2) is solved.
 \endverbatim

 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrices A, B, Y and X.  M >= 0.
 \endverbatim

 \param[in,out] A
 \verbatim
          A is SINGLE PRECISION array, dimension (LDA,M)
          If FACT == "N", the matrix A is a general matrix and it is overwritten with the
          quasi upper triangular matrix S of the generalized schur decomposition.
          If FACT == "F", the matrix A contains its (quasi-) upper triangular matrix S of the
          generalized  Schur decomposition of (A,B).
          If FACT == "H", the matrix A is an upper Hessenberg matrix of the generalized
          Hessenberg form (A,B) and it is overwritten with the quasi upper triangular matrix S
          of the generalized Schur decomposition.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim

 \param[in,out] B
 \verbatim
          B is SINGLE PRECISION array, dimension (LDB,M)
          If FACT == "N", the matrix B a general matrix and it is overwritten with the upper triangular
          matrix of the generalized Schur decomposition.
          If FACT == "F", the matrix B contains its upper triangular matrix R of the generalized schur
          Schur decomposition of (A,B).
          If FACT == "H", the matrix B is the upper triangular matrix of the generalized Hessenberg form
          (A,B) and it is overwritten with the upper triangular matrix of the generalized Schur decomposition.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,M).
 \endverbatim

 \param[in,out] Q
 \verbatim
          Q is SINGLE PRECISION array, dimension (LDQ,M)
          If FACT == "N", the matrix Q is an empty M-by-M matrix on input and contains the
          left Schur vectors of (A,B) on output.
          If FACT == "F", the matrix Q contains the left Schur vectors of (A,B).
          If FACT == "H", the matrix Q is an empty M-by-M matrix on input and contains the
          left Schur vectors of (A,B) on output.
 \endverbatim

 \param[in] LDQ
 \verbatim
          LDQ is INTEGER
          The leading dimension of the array Q.  LDQ >= max(1,M).
 \endverbatim


 \param[in,out] Z
 \verbatim
          Z is SINGLE PRECISION array, dimension (LDZ,M)
          If FACT == "N", the matrix Z is an empty M-by-M matrix on input and contains the
          right Schur vectors of (A,B) on output.
          If FACT == "F", the matrix Z contains the right Schur vectors of (A,B).
          If FACT == "H", the matrix Z is an empty M-by-M matrix on input and contains the
          right Schur vectors of (A,B) on output.
 \endverbatim

 \param[in] LDZ
 \verbatim
          LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= max(1,M).
 \endverbatim

 \param[in,out] X
 \verbatim
          X is SINGLE PRECISION array, dimension (LDX,N)
          On input, the matrix X contains the right hand side Y.
          On output, the matrix X contains the solution of Equation (1) or (2)
          Right hand side Y and the solution X are symmetric M-by-M matrices.
 \endverbatim

 \param[in] LDX
 \verbatim
          LDX is INTEGER
          The leading dimension of the array X.  LDX >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is SINGLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
 \endverbatim

 \param[in,out] WORK
 \verbatim
          WORK is SINGLE PRECISION array, dimension (MAX(1,LDWORK))
          Workspace for the algorithm. The optimal workspace is given either by \ref mepack_memory_frontend.
 \endverbatim

 \param[in] LDWORK
 \verbatim
          LDWORK is INTEGER
          Size of the workspace for the algorithm counted in floating point numbers of the actual precision.
          The C interface does not support the workspace query by setting LDWORK == -1 on input. In this case,
          the \ref mepack_memory_frontend function have to be used.
 \endverbatim


 \param[out] INFO
 \verbatim
          INFO is INTEGER
          == 0:  successful exit
          = 1:  DGGES failed
          = 2:  DLA_SORT_GEV failed
          = 3:  Inner solver failed
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date October 2023
*/


void mepack_single_ggstein(const char * FACT, const char *TRANS, int M, float * A,  int LDA, float *B, int LDB,
        float * Q, int LDQ, float *Z, int LDZ, float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int *INFO)
{
    Int _M = M;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDQ = LDQ;
    Int _LDZ = LDZ;
    Int _LDX = LDX;
    Int _INFO = * INFO;
    Int _LDWORK = LDWORK;
    FORTRAN_LEN_T c1len = strlen(FACT);
    FORTRAN_LEN_T c2len = strlen(TRANS);

    FC_GLOBAL_(sla_ggstein,SLA_GGSTEIN)(FACT, TRANS, &_M, A, &_LDA, B, &_LDB, Q, &_LDQ, Z, &_LDZ, X, &_LDX, SCALE, WORK, &_LDWORK, &_INFO, c1len, c2len);

    *INFO = _INFO;
    return;

}








/**
 * @}
 */

