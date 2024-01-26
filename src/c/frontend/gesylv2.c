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
 * @addtogroup cgesylv
 * @{
 */

/**
 \brief Frontend for the solution of Standard Sylvester Equations.

 \par Purpose:

 \verbatim
 mepack_double_gesylv2 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) + X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) - X = SCALE * Y                               (2)

 where A is a M-by-M matrix and B is a N-by-N matrix. The right hand
 side Y and the solution X are M-by-N matrices. The matrices A and B can be
 either a general unreduced matrix or an upper Hessenberg form
 or a (quasi-) upper triangular factor. In the latter case QA and QB provide
 the Schur-vectors of the matrices A and B.
 \endverbatim

 \remark This function is a wrapper for \ref dla_gesylv2.

 \see dla_gesylv2

 \param[in] FACTA
 \verbatim
          FACTA is String
          Specifies how the matrix A is given.
          == 'N':  The matrix A is given as a general matrix and its Schur decomposition
                  A = QA*S*QA**T will be computed.
          == 'F':  The matrix A is given as its Schur decomposition in terms of S and QA
                  form A = QA*S*QA**T
          == 'H':  The matrix A is given in upper Hessenberg form and its Schur decomposition
                  A = QA*S*QA**T will be computed
 \endverbatim

 \param[in] FACTB
 \verbatim
          FACTB is String
          Specifies how the matrix B is given.
          == 'N':  The matrix B is given as a general matrix and its Schur decomposition
                  B = QB*R*QB**T will be computed.
          == 'F':  The matrix B is given as its Schur decomposition in terms of R and QB
                  form B = QB*R*QB**T
          == 'H':  The matrix B is given in upper Hessenberg form and its Schur decomposition
                  B = QB*R*QB**T will be computed
 \endverbatim

 \param[in] TRANSA
 \verbatim
          TRANSA is String
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is String
          Specifies the form of the system of equations with respect to B:
          == 'N':  op2(B) = B,
          == 'T':  op2(B) = B**T
 \endverbatim

 \param[in] SGN
 \verbatim
          SGN is DOUBLE PRECISION, allowed values: +/-1
          Specifies the sign between the two parts of the Sylvester equation.
          = 1 :  Solve Equation (1)
          == -1:  Solve Equation (2)
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

 \param[in,out] A
 \verbatim
          A is DOUBLE PRECISION array, dimension (LDA,M)
          If FACTA == "N", the matrix A is a general matrix and it is overwritten with its
          schur decomposition S.
          If FACTA == "F", the matrix A contains its (quasi-) upper triangular matrix S being the
          Schur decomposition of A.
          If FACTA == "H", the matrix A is an upper Hessenberg matrix and it is overwritten
          with its schur decomposition S.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim


 \param[in,out] B
 \verbatim
          B is DOUBLE PRECISION array, dimension (LDB,N)
          If FACTB == "N", the matrix B is a general matrix and it is overwritten with its
          schur decomposition R.
          If FACTB == "F", the matrix B contains its (quasi-) upper triangular matrix R beeping the
          Schur decomposition of B.
          If FACTB == "H", the matrix B is an upper Hessenberg matrix and it is overwritten with its
          schur decomposition R.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
 \endverbatim

 \param[in,out] QA
 \verbatim
          QA is DOUBLE PRECISION array, dimension (LDQA,M)
          If FACTA == "N", the matrix QA is an empty M-by-M matrix on input and contains the
          Schur vectors of A on output.
          If FACTA == "F", the matrix QA contains the Schur vectors of A.
          If FACTA == "H", the matrix QA is an empty M-by-M matrix on input and contains the
          Schur vectors of A on output.
 \endverbatim

 \param[in] LDQA
 \verbatim
          LDQA is INTEGER
          The leading dimension of the array QA.  LDQA >= max(1,M).
 \endverbatim


 \param[in,out] QB
 \verbatim
          QB is DOUBLE PRECISION array, dimension (LDQB,N)
          If FACTB == "N", the matrix QB is an empty N-by-N matrix on input and contains the
          Schur vectors of B on output.
          If FACTB == "F", the matrix QB contains the Schur vectors of B.
          QB is DOUBLE PRECISION array, dimension (LDQB,N)
          If FACTB == "H", the matrix QB is an empty N-by-N matrix on input and contains the
          Schur vectors of B on output.
 \endverbatim

 \param[in] LDQB
 \verbatim
          LDQB is INTEGER
          The leading dimension of the array QB.  LDQB >= max(1,N).
 \endverbatim

 \param[in,out] X
 \verbatim
          X is DOUBLE PRECISION array, dimension (LDX,N)
          On input, the matrix X contains the right hand side Y.
          On output, the matrix X contains the solution of Equation (1) or (2)
          Right hand side Y and the solution X are M-by-N matrices.
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
          = 1:  DHGEES failed
          = 2:  DLA_SORT_EV failed
          = 3:  DLA_TRLYAP_DAG failed
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date January 2024
*/
void mepack_double_gesylv2(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, double SGN,  int M, int N,
        double * A, int LDA,double * B, int LDB, double *QA, int LDQA, double *QB, int LDQB, double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int* INFO)
{
    Int _M = M;
    Int _N = N;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDQA = LDQA;
    Int _LDQB = LDQB;
    Int _LDX  = LDX;
    Int _LDWORK = LDWORK;
    Int _INFO = *INFO;
    FORTRAN_LEN_T c1len = strlen(FACTA);
    FORTRAN_LEN_T c2len = strlen(FACTB);
    FORTRAN_LEN_T c3len = strlen(TRANSA);
    FORTRAN_LEN_T c4len = strlen(TRANSB);

    FC_GLOBAL_(dla_gesylv2,DLA_GESYLV2)( FACTA, FACTB, TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, QA, &_LDQA, QB, &_LDQB, X, &_LDX,
            SCALE, WORK, &_LDWORK, &_INFO,
            c1len, c2len, c3len, c4len);
    *INFO = (int) _INFO;
    return;

}


/**
 \brief Frontend for the solution of Standard Sylvester Equations.

 \par Purpose:

 \verbatim
 mepack_single_gesylv2 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X = SCALE * Y                              (2)

 where A is a M-by-M matrix and B is a N-by-N matrix. The right hand
 side Y and the solution X are M-by-N matrices. The matrices A and B can be
 either a general unreduced matrix or an upper Hessenberg form
 or a (quasi-) upper triangular factor. In the latter case QA and QB provide
 the Schur-vectors of the matrices A and B.
 \endverbatim

 \remark This function is a wrapper for \ref sla_gesylv2.

 \see dla_gesylv2

 \param[in] FACTA
 \verbatim
          FACTA is String
          Specifies how the matrix A is given.
          == 'N':  The matrix A is given as a general matrix and its Schur decomposition
                  A = QA*S*QA**T will be computed.
          == 'F':  The matrix A is given as its Schur decomposition in terms of S and QA
                  form A = QA*S*QA**T
          == 'H':  The matrix A is given in upper Hessenberg form and its Schur decomposition
                  A = QA*S*QA**T will be computed
 \endverbatim

 \param[in] FACTB
 \verbatim
          FACTB is String
          Specifies how the matrix B is given.
          == 'N':  The matrix B is given as a general matrix and its Schur decomposition
                  B = QB*R*QB**T will be computed.
          == 'F':  The matrix B is given as its Schur decomposition in terms of R and QB
                  form B = QB*R*QB**T
          == 'H':  The matrix B is given in upper Hessenberg form and its Schur decomposition
                  B = QB*R*QB**T will be computed
 \endverbatim

 \param[in] TRANSA
 \verbatim
          TRANSA is String
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is String
          Specifies the form of the system of equations with respect to B:
          == 'N':  op2(B) = B,
          == 'T':  op2(B) = B**T
 \endverbatim

 \param[in] SGN
 \verbatim
          SGN is DOUBLE PRECISION, allowed values: +/-1
          Specifies the sign between the two parts of the Sylvester equation.
          = 1 :  Solve Equation (1)
          == -1:  Solve Equation (2)
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

 \param[in,out] A
 \verbatim
          A is DOUBLE PRECISION array, dimension (LDA,M)
          If FACTA == "N", the matrix A is a general matrix and it is overwritten with its
          schur decomposition S.
          If FACTA == "F", the matrix A contains its (quasi-) upper triangular matrix S being the
          Schur decomposition of A.
          If FACTA == "H", the matrix A is an upper Hessenberg matrix and it is overwritten
          with its schur decomposition S.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim


 \param[in,out] B
 \verbatim
          B is DOUBLE PRECISION array, dimension (LDB,N)
          If FACTB == "N", the matrix B is a general matrix and it is overwritten with its
          schur decomposition R.
          If FACTB == "F", the matrix B contains its (quasi-) upper triangular matrix R being the
          Schur decomposition of B.
          If FACTB == "H", the matrix B is an upper Hessenberg matrix and it is overwritten with its
          schur decomposition R.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
 \endverbatim

 \param[in,out] QA
 \verbatim
          QA is DOUBLE PRECISION array, dimension (LDA,M)
          If FACTA == "N", the matrix QA is an empty M-by-M matrix on input and contains the
          Schur vectors of A on output.
          If FACTA == "F", the matrix QA contains the Schur vectors of A.
          If FACTA == "H", the matrix QA is an empty M-by-M matrix on input and contains the
          Schur vectors of A on output.
 \endverbatim

 \param[in] LDQA
 \verbatim
          LDQA is INTEGER
          The leading dimension of the array QA.  LDQA >= max(1,M).
 \endverbatim


 \param[in,out] QB
 \verbatim
          QB is DOUBLE PRECISION array, dimension (LDA,M)
          If FACTB == "N", the matrix QB is an empty M-by-M matrix on input and contains the
          Schur vectors of B on output.
          If FACTB == "F", the matrix QB contains the Schur vectors of B.
          If FACTB == "H", the matrix QB is an empty M-by-M matrix on input and contains the
          Schur vectors of B on output.
 \endverbatim

 \param[in] LDQB
 \verbatim
          LDQB is INTEGER
          The leading dimension of the array QB.  LDQB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
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
          = 1:  DHGEES failed
          = 2:  DLA_SORT_EV failed
          = 3:  DLA_TRLYAP_DAG failed
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date January 2024
*/
void mepack_single_gesylv2(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, float SGN,  int M, int N,
        float * A, int LDA,float * B, int LDB, float *QA, int LDQA, float *QB, int LDQB, float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int* INFO)
{
    Int _M = M;
    Int _N = N;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDQA = LDQA;
    Int _LDQB = LDQB;
    Int _LDX  = LDX;
    Int _LDWORK = LDWORK;
    Int _INFO = *INFO;
    FORTRAN_LEN_T c1len = strlen(FACTA);
    FORTRAN_LEN_T c2len = strlen(FACTB);
    FORTRAN_LEN_T c3len = strlen(TRANSA);
    FORTRAN_LEN_T c4len = strlen(TRANSB);

    FC_GLOBAL_(sla_gesylv2,SLA_GESYLV2)( FACTA, FACTB, TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, QA, &_LDQA, QB, &_LDQB, X, &_LDX,
            SCALE, WORK, &_LDWORK, &_INFO,
            c1len, c2len, c3len, c4len);
    *INFO = (int) _INFO;
    return;

}






/**
 * @}
 */

