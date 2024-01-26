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
 * @addtogroup cggsylv
 * @{
 */
/**
 \brief Frontend for the solution of Generalized Sylvester Equations.

 \par Purpose:

 \verbatim
 mepack_double_ggsylv solves a generalized Sylvester equation of the following forms

    op1(A) * X * op2(B) + op1(C) * X * op2(D) = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) - op1(C) * X * op2(D) = SCALE * Y                              (2)

 where (A,C) is a M-by-M matrix pencil and (B,D) is a N-by-N matrix pencil.
 The right hand side Y and the solution X M-by-N matrices. The matrix pencils (A,C)
 and (B,D) can be either given as general unreduced matrices, as generalized
 Hessenberg form, or in terms of their generalized Schur decomposition. 
 If they are given as general matrices or as a generalized Hessenberg form
 their generalized Schur decomposition will be computed.

 \endverbatim

 \remark This function is a wrapper for \ref dla_ggsylv.

 \see dla_ggsylv

 \param[in] FACTA
 \verbatim
          FACTA is String
          Specifies how the matrix pencil (A,C) is given.
          == 'N':  The matrix pencil (A,C) is given as a general matrices and its Schur decomposition
                  A = QA*S*ZA**T, C = QA*R*ZA**T will be computed.
          == 'F':  The matrix pencil (A,C) is already in generalized Schur form and S, R, QA, and ZA
                  are given.
          == 'H': The matrix pencil (A,C) is given in generalized Hessenberg form and its Schur decomposition
                  A = QA*S*ZA**T, C = QA*R*ZA**T will be computed.
 \endverbatim

 \param[in] FACTB
 \verbatim
          FACTB is String
          Specifies how the matrix pencil (B,D) is given.
          == 'N':  The matrix pencil (B,D) is given as a general matrices and its Schur decomposition
                  B = QB*U*ZB**T, D = QB*V*ZB**T will be computed.
          == 'F':  The matrix pencil (B,D) is already in generalized Schur form and U, V, QB, and ZB
                  are given.
          == 'H': The matrix pencil (B,D) is given in generalized Hessenberg form and its Schur decomposition
                  B = QB*U*ZB**T, D = QB*V*ZB**T will be computed.
 \endverbatim

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
          If FACT == "N", the matrix A is a general matrix and it is overwritten with the
          (quasi-) upper triangular factor S of the Schur decomposition of (A,C).
          If FACT == "F", the matrix A contains its (quasi-) upper triangular matrix S of
          the Schur decomposition of (A,C).
          If FACT == "H", the matrix A is an upper Hessenberg matrix of the generalized
          Hessenberg form (A,C) and it is overwritten with the (quasi-) upper triangular
          factor S of the Schur decomposition of (A,C).
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim


 \param[in,out] B
 \verbatim
          B is DOUBLE PRECISION array, dimension (LDB,N)
          If FACT == "N",  the matrix B is a general matrix and it is overwritten with the
          (quasi-) upper triangular factor U of the Schur decomposition of (B,D).
          If FACT == "F", the matrix B contains its (quasi-) upper triangular matrix U of
          the Schur decomposition of (B,D).
          If FACT == "H", the matrix B is an upper Hessenberg matrix of the generalized
          Hessenberg form (B,D) and it is overwritten with the (quasi-) upper triangular
          factor U of the Schur decomposition of (B,D).
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
 \endverbatim

 \param[in,out] C
 \verbatim
          C is DOUBLE PRECISION array, dimension (LDC,M)
          If FACT == "N", the matrix C is a general matrix and it is overwritten with the
          upper triangular factor R of the Schur decomposition of (A,C).
          If FACT == "F", the matrix C contains its upper triangular matrix R of
          the Schur decomposition of (A,C).
          If FACT == "H", the matrix C is the upper triangular matrix of the generalized Hessenberg form
          (A,C) and it is overwritten with the upper triangular factor R of the Schur decomposition of (A,C).
 \endverbatim

 \param[in] LDC
 \verbatim
          LDC is INTEGER
          The leading dimension of the array C.  LDC >= max(1,M).
 \endverbatim


 \param[in,out] D
 \verbatim
          D is DOUBLE PRECISION array, dimension (LDD,N)
          If FACT == "N",  the matrix D is a general matrix and it is overwritten with the
          upper triangular factor V of the Schur decomposition of (B,D).
          If FACT == "F", the matrix D contains its upper triangular matrix V of
          the Schur decomposition of (B,D).
          If FACT == "H", the matrix D is the upper triangular matrix of the generalized Hessenberg form
          (B,D) and it is overwritten with the upper triangular factor V of the Schur decomposition of (B,D).
 \endverbatim

 \param[in] LDD
 \verbatim
          LDD is INTEGER
          The leading dimension of the array D.  LDD >= max(1,N).
 \endverbatim

 \param[in,out] QA
 \verbatim
          QA is DOUBLE PRECISION array, dimension (LDQA,M)
          If FACT == "N", the matrix QA is an empty M-by-M matrix on input and contains the
          left Schur vectors of (A,C) on output.
          If FACT == "F", the matrix QA contains the left Schur vectors of (A,C).
          If FACT == "H", the matrix QA is an empty M-by-M matrix on input and contains the
          left Schur vectors of (A,C) on output.
 \endverbatim

 \param[in] LDQA
 \verbatim
          LDQA is INTEGER
          The leading dimension of the array QA.  LDQA >= max(1,M).
 \endverbatim

 \param[in,out] ZA
 \verbatim
          ZA is DOUBLE PRECISION array, dimension (LDZA,M)
          If FACT == "N", the matrix ZA is an empty M-by-M matrix on input and contains the
          right Schur vectors of (A,C) on output.
          If FACT == "F", the matrix ZA contains the right Schur vectors of (A,C).
          If FACT == "H", the matrix ZA is an empty M-by-M matrix on input and contains the
          right Schur vectors of (A,C) on output.
 \endverbatim

 \param[in] LDZA
 \verbatim
          LDZA is INTEGER
          The leading dimension of the array ZA.  LDZA >= max(1,M).
 \endverbatim

 \param[in,out] QB
 \verbatim
          QB is DOUBLE PRECISION array, dimension (LDQB,N)
          If FACT == "N", the matrix QB is an empty N-by-N matrix on input and contains the
          left Schur vectors of (B,D) on output.
          If FACT == "F", the matrix QB contains the left Schur vectors of (B,D).
          If FACT == "H", the matrix QB is an empty M-by-M matrix on input and contains the
          left Schur vectors of (B,D) on output.
 \endverbatim

 \param[in] LDQB
 \verbatim
          LDQB is INTEGER
          The leading dimension of the array QB.  LDQB >= max(1,N).
 \endverbatim

 \param[in,out] ZB
 \verbatim
          ZB is DOUBLE PRECISION array, dimension (LDZB,N)
          If FACT == "N", the matrix ZB is an empty N-by-N matrix on input and contains the
          right Schur vectors of (B,D) on output.
          If FACT == "F", the matrix ZB contains the right Schur vectors of (B,D).
          If FACT == "H", the matrix ZB is an empty M-by-M matrix on input and contains the
          right Schur vectors of (B,D) on output.
 \endverbatim

 \param[in] LDZB
 \verbatim
          LDZB is INTEGER
          The leading dimension of the array ZB.  LDZB >= max(1,N).
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
          Workspace for the algorithm. The optimal workspace is given by \ref mepack_memory_frontend.
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
          = 1:  DHGGES failed
          = 2:  DLA_SORT_GEV failed
          = 3:  Inner solver failed
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date January 2024
*/
void mepack_double_ggsylv(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, double SGN,  int M, int N,
        double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
        double *QA, int LDQA, double *ZA, int LDZA, double *QB, int LDQB, double *ZB, int LDZB,
        double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int* INFO)
{
    Int _M = M;
    Int _N = N;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDC = LDC;
    Int _LDD = LDD;
    Int _LDQA = LDQA;
    Int _LDQB = LDQB;
    Int _LDZA = LDZA;
    Int _LDZB = LDZB;
    Int _LDX  = LDX;
    Int _LDWORK = LDWORK;
    Int _INFO = *INFO;
    FORTRAN_LEN_T c1len = strlen(FACTA);
    FORTRAN_LEN_T c2len = strlen(FACTB);
    FORTRAN_LEN_T c3len = strlen(TRANSA);
    FORTRAN_LEN_T c4len = strlen(TRANSB);

    FC_GLOBAL_(dla_ggsylv,DLA_GGSYLV)( FACTA, FACTB, TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, C, &_LDC, D, &_LDD,
            QA, &_LDQA, ZA, &_LDZA, QB, &_LDQB, ZB, &_LDZB, X, &_LDX,
            SCALE, WORK, &_LDWORK, &_INFO,
            c1len, c2len, c3len, c4len);
    *INFO = (int) _INFO;
    return;

}

/**
 \brief Frontend for the solution of Generalized Sylvester Equations.

 \par Purpose:

 \verbatim
 mepack_single_ggsylv solves a generalized Sylvester equation of the following forms

    op1(A) * X * op2(B) + op1(C) * X * op2(D) = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) - op1(C) * X * op2(D) = SCALE * Y                              (2)

 where (A,C) is a M-by-M matrix pencil and (B,D) is a N-by-N matrix pencil.
 The right hand side Y and the solution X M-by-N matrices. The matrix pencils (A,C)
 and (B,D) can be either given as general unreduced matrices, as generalized
 Hessenberg form, or in terms of their generalized Schur decomposition. 
 If they are given as general matrices or as a generalized Hessenberg form
 their generalized Schur decomposition will be computed.

 \endverbatim

 \remark This function is a wrapper for \ref sla_ggsylv.

 \see sla_ggsylv

 \param[in] FACTA
 \verbatim
          FACTA is String
          Specifies how the matrix pencil (A,C) is given.
          == 'N':  The matrix pencil (A,C) is given as a general matrices and its Schur decomposition
                  A = QA*S*ZA**T, C = QA*R*ZA**T will be computed.
          == 'F':  The matrix pencil (A,C) is already in generalized Schur form and S, R, QA, and ZA
                  are given.
          == 'H': The matrix pencil (A,C) is given in generalized Hessenberg form and its Schur decomposition
                  A = QA*S*ZA**T, C = QA*R*ZA**T will be computed.
 \endverbatim

 \param[in] FACTB
 \verbatim
          FACTB is String
          Specifies how the matrix pencil (B,D) is given.
          == 'N':  The matrix pencil (B,D) is given as a general matrices and its Schur decomposition
                  B = QB*U*ZB**T, D = QB*V*ZB**T will be computed.
          == 'F':  The matrix pencil (B,D) is already in generalized Schur form and U, V, QB, and ZB
                  are given.
          == 'H': The matrix pencil (B,D) is given in generalized Hessenberg form and its Schur decomposition
                  B = QB*U*ZB**T, D = QB*V*ZB**T will be computed.
 \endverbatim

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

 \param[in] SGN
 \verbatim
          SGN is SINGLE PRECISION, allowed values: +/-1
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
          A is SINGLE PRECISION array, dimension (LDA,M)
          If FACT == "N", the matrix A is a general matrix and it is overwritten with the
          (quasi-) upper triangular factor S of the Schur decomposition of (A,C).
          If FACT == "F", the matrix A contains its (quasi-) upper triangular matrix S of
          the Schur decomposition of (A,C).
          If FACT == "H", the matrix A is an upper Hessenberg matrix of the generalized
          Hessenberg form (A,C) and it is overwritten with the (quasi-) upper triangular
          factor S of the Schur decomposition of (A,C).
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim


 \param[in,out] B
 \verbatim
          B is SINGLE PRECISION array, dimension (LDB,N)
          If FACT == "N",  the matrix B is a general matrix and it is overwritten with the
          (quasi-) upper triangular factor U of the Schur decomposition of (B,D).
          If FACT == "F", the matrix B contains its (quasi-) upper triangular matrix U of
          the Schur decomposition of (B,D).
          If FACT == "H", the matrix B is an upper Hessenberg matrix of the generalized
          Hessenberg form (B,D) and it is overwritten with the (quasi-) upper triangular
          factor U of the Schur decomposition of (B,D).
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
 \endverbatim

 \param[in,out] C
 \verbatim
          C is SINGLE PRECISION array, dimension (LDC,M)
          If FACT == "N", the matrix C is a general matrix and it is overwritten with the
          upper triangular factor R of the Schur decomposition of (A,C).
          If FACT == "F", the matrix C contains its upper triangular matrix R of
          the Schur decomposition of (A,C).
          If FACT == "H", the matrix C is the upper triangular matrix of the generalized Hessenberg form
          (A,C) and it is overwritten with the upper triangular factor R of the Schur decomposition of (A,C).
 \endverbatim

 \param[in] LDC
 \verbatim
          LDC is INTEGER
          The leading dimension of the array C.  LDC >= max(1,M).
 \endverbatim


 \param[in,out] D
 \verbatim
          D is SINGLE PRECISION array, dimension (LDD,N)
          If FACT == "N",  the matrix D is a general matrix and it is overwritten with the
          upper triangular factor V of the Schur decomposition of (B,D).
          If FACT == "F", the matrix D contains its upper triangular matrix V of
          the Schur decomposition of (B,D).
          If FACT == "H", the matrix D is the upper triangular matrix of the generalized Hessenberg form
          (B,D) and it is overwritten with the upper triangular factor V of the Schur decomposition of (B,D).
 \endverbatim

 \param[in] LDD
 \verbatim
          LDD is INTEGER
          The leading dimension of the array D.  LDD >= max(1,N).
 \endverbatim

 \param[in,out] QA
 \verbatim
          QA is SINGLE PRECISION array, dimension (LDQA,M)
          If FACT == "N", the matrix QA is an empty M-by-M matrix on input and contains the
          left Schur vectors of (A,C) on output.
          If FACT == "F", the matrix QA contains the left Schur vectors of (A,C).
          If FACT == "H", the matrix QA is an empty M-by-M matrix on input and contains the
          left Schur vectors of (A,C) on output.
 \endverbatim

 \param[in] LDQA
 \verbatim
          LDQA is INTEGER
          The leading dimension of the array QA.  LDQA >= max(1,M).
 \endverbatim

 \param[in,out] ZA
 \verbatim
          ZA is SINGLE PRECISION array, dimension (LDZA,M)
          If FACT == "N", the matrix ZA is an empty M-by-M matrix on input and contains the
          right Schur vectors of (A,C) on output.
          If FACT == "F", the matrix ZA contains the right Schur vectors of (A,C).
          If FACT == "H", the matrix ZA is an empty M-by-M matrix on input and contains the
          right Schur vectors of (A,C) on output.
 \endverbatim

 \param[in] LDZA
 \verbatim
          LDZA is INTEGER
          The leading dimension of the array ZA.  LDZA >= max(1,M).
 \endverbatim

 \param[in,out] QB
 \verbatim
          QB is SINGLE PRECISION array, dimension (LDQB,M)
          If FACT == "N", the matrix QB is an empty M-by-M matrix on input and contains the
          left Schur vectors of (B,D) on output.
          If FACT == "F", the matrix QB contains the left Schur vectors of (B,D).
          If FACT == "H", the matrix QB is an empty M-by-M matrix on input and contains the
          left Schur vectors of (B,D) on output.
 \endverbatim

 \param[in] LDQB
 \verbatim
          LDQB is INTEGER
          The leading dimension of the array QB.  LDQB >= max(1,M).
 \endverbatim

 \param[in,out] ZB
 \verbatim
          ZB is SINGLE PRECISION array, dimension (LDZB,M)
          If FACT == "N", the matrix ZB is an empty M-by-M matrix on input and contains the
          right Schur vectors of (B,D) on output.
          If FACT == "F", the matrix ZB contains the right Schur vectors of (B,D).
          If FACT == "H", the matrix ZB is an empty M-by-M matrix on input and contains the
          right Schur vectors of (B,D) on output.
 \endverbatim

 \param[in] LDZB
 \verbatim
          LDZB is INTEGER
          The leading dimension of the array ZB.  LDZB >= max(1,M).
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
          The leading dimension of the array X.  LDB >= max(1,M).
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
          Workspace for the algorithm. The optimal workspace is given by \ref mepack_memory_frontend.
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
          = 1:  DHGGES failed
          = 2:  DLA_SORT_GEV failed
          = 3:  Inner solver failed
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date January 2024
*/
void mepack_single_ggsylv(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, float SGN,  int M, int N,
        float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
        float *QA, int LDQA, float *ZA, int LDZA, float *QB, int LDQB, float *ZB, int LDZB,
        float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int* INFO)
{
    Int _M = M;
    Int _N = N;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDC = LDC;
    Int _LDD = LDD;
    Int _LDQA = LDQA;
    Int _LDQB = LDQB;
    Int _LDZA = LDZA;
    Int _LDZB = LDZB;
    Int _LDX  = LDX;
    Int _LDWORK = LDWORK;
    Int _INFO = *INFO;
    FORTRAN_LEN_T c1len = strlen(FACTA);
    FORTRAN_LEN_T c2len = strlen(FACTB);
    FORTRAN_LEN_T c3len = strlen(TRANSA);
    FORTRAN_LEN_T c4len = strlen(TRANSB);

    FC_GLOBAL_(sla_ggsylv,SLA_GGSYLV)( FACTA, FACTB, TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, C, &_LDC, D, &_LDD,
            QA, &_LDQA, ZA, &_LDZA, QB, &_LDQB, ZB, &_LDZB, X, &_LDX,
            SCALE, WORK, &_LDWORK, &_INFO,
            c1len, c2len, c3len, c4len);
    *INFO = (int) _INFO;
    return;

}




/**
 * @}
 */

