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
 \brief Frontend for the solution of the dual Coupled Generalized Sylvester Equations.
 \par Purpose:

 \verbatim

 mepack_double_ggcsylv_dual solves a generalized coupled  Sylvester equation of the following form

    op1(A)**T * R  + op1(C)**T * L               = SCALE * E                             (1)
    SGN1 * R * op2(B)**T + SGN2 * L * op2(D)** T = SCALE * F

 where A and C are M-by-M matrices and B and D are N-by-N matrices.
 The right hand sides E, F  and the solutions R, L are  M-by-N matrices.
 The equation (1) is the dual to the generalized coupled Sylvester equation

    op1(A) * R + SGN1 * L * op2(B)  = SCALE * E                                          (2)
    op1(C) * R + SGN2 * L * op2(D)  = SCALE * F

 The equation (1) is the dual one to equation (2) with respect to the underlying linear system.
 Let Z be the matrix formed by rewriting (2) into its Kronecker form. This yields

         | kron(I, op1(A))   SGN1*kron(op2(B)**T, I) | | Vec R |   | Vec E |
   Z X = |                                           |*|       | = |       |
         | kron(I, op1(C))   SGN2*kron(op2(D)**T, I) | | Vec L |   | Vec F |

 Regarding Z**T one obtains

            | kron(I, op1(A)**T )    kron(I, op1(C)**T)   | | Vec R |   | Vec E |
   Z**T X = |                                             |*|       | = |       |
            | SGN1*kron(op2(B), I)   SGN2*kron(op2(D), I) | | Vec L |   | Vec F |

 which belongs to the Sylvester equation (1). For this reason the parameters TRANSA and TRANSB
 are expressed in terms of the Sylvester equation (2).

 \endverbatim
 \remark This function is a wrapper for \ref dla_ggcsylv_dual.

 \see dla_ggcsylv_dual

 \param[in] FACTA
 \verbatim
          FACTA is String
          Specifies how the matrix pencil (A,C) is given.
          == 'N':  The matrix pencil (A,C) is given as a general matrices and its Schur decomposition
                  A = QA*S*ZA**T, C = QA*R*ZA**T will be computed.
          == 'F':  The matrix pencil (A,C) is already in generalized Schur form and S, R, QA, and ZA
                  are given.
 \endverbatim

 \param[in] FACTB
 \verbatim
          FACTB is String
          Specifies how the matrix pencil (B,D) is given.
          == 'N':  The matrix pencil (B,D) is given as a general matrices and its Schur decomposition
                  B = QB*U*ZB**T, D = QB*V*ZB**T will be computed.
          == 'F':  The matrix pencil (B,D) is already in generalized Schur form and U, V, QB, and ZB
                  are given.
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

 \param[in] SGN1
 \verbatim
          SGN1 is DOUBLE PRECISION, allowed values: +/-1
          Specifies the sign between in the first equation.
 \endverbatim

 \param[in] SGN2
 \verbatim
          SGN2 is DOUBLE PRECISION, allowed values: +/-1
          Specifies the sign between in the second equation.
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
 \endverbatim

 \param[in] LDZB
 \verbatim
          LDZB is INTEGER
          The leading dimension of the array ZB.  LDZB >= max(1,N).
 \endverbatim

 \param[in,out] E
 \verbatim
          E is DOUBLE PRECISION array, dimension (LDE,N)
          On input, the matrix E contains the right hand side E.
          On output, the matrix E contains the solution R.
 \endverbatim

 \param[in] LDE
 \verbatim
          LDE is INTEGER
          The leading dimension of the array E.  LDE >= max(1,M).
 \endverbatim

 \param[in,out] F
 \verbatim
          F is DOUBLE PRECISION array, dimension (LDF,N)
          On input, the matrix F contains the right hand side F.
          On output, the matrix F contains the solution L.
 \endverbatim

 \param[in] LDF
 \verbatim
          LDF is INTEGER
          The leading dimension of the array F.  LDF >= max(1,M).
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
          = 1:  DGGES failed
          = 2:  DLA_SORT_GEV failed
          = 3:  Inner solver failed
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date June 2023
*/

void mepack_double_ggcsylv_dual(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
        double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
        double *QA, int LDQA, double *ZA, int LDZA, double *QB, int LDQB, double *ZB, int LDZB,
        double *E, int LDE, double *F, int LDF,
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
    Int _LDE  = LDE;
    Int _LDF  = LDF;
    Int _LDWORK = LDWORK;
    Int _INFO = *INFO;
    FORTRAN_LEN_T c1len = strlen(FACTA);
    FORTRAN_LEN_T c2len = strlen(FACTB);
    FORTRAN_LEN_T c3len = strlen(TRANSA);
    FORTRAN_LEN_T c4len = strlen(TRANSB);

    FC_GLOBAL_(dla_ggcsylv_dual,DLA_GGCSYLV_DUAL)( FACTA, FACTB, TRANSA, TRANSB, &SGN1, &SGN2, &_M, &_N,
            A, &_LDA, B, &_LDB, C, &_LDC, D, &_LDD,
            QA, &_LDQA, ZA, &_LDZA, QB, &_LDQB, ZB, &_LDZB, E, &_LDE, F, &_LDF,
            SCALE, WORK, &_LDWORK, &_INFO,
            c1len, c2len, c3len, c4len);
    *INFO = (int) _INFO;
    return;

}

/**
 \brief Frontend for the solution of the dual Coupled Generalized Sylvester Equations.
 \par Purpose:

 \verbatim

 mepack_single_ggcsylv_dual solves a generalized coupled  Sylvester equation of the following form

    op1(A)**T * R  + op1(C)**T * R               = SCALE * E                             (1)
    SGN1 * R * op2(B)**T + SGN2 * L * op2(D)** T = SCALE * F

 where A and C are M-by-M matrices and B and D are N-by-N matrices.
 The right hand sides E, F  and the solutions R, L are  M-by-N matrices.
 The equation (1) is the dual to the generalized coupled Sylvester equation

    op1(A) * R + SGN1 * L * op2(B)  = SCALE * E                                          (2)
    op1(C) * R + SGN2 * L * op2(D)  = SCALE * F

 The equation (1) is the dual one to equation (2) with respect to the underlying linear system.
 Let Z be the matrix formed by rewriting (2) into its Kronecker form. This yields

         | kron(I, op1(A))   SGN1*kron(op2(B)**T, I) | | Vec R |   | Vec E |
   Z X = |                                           |*|       | = |       |
         | kron(I, op1(C))   SGN2*kron(op2(D)**T, I) | | Vec L |   | Vec F |

 Regarding Z**T one obtains

            | kron(I, op1(A)**T )    kron(I, op1(C)**T)   | | Vec R |   | Vec E |
   Z**T X = |                                             |*|       | = |       |
            | SGN1*kron(op2(B), I)   SGN2*kron(op2(D), I) | | Vec L |   | Vec F |

 which belongs to the Sylvester equation (1). For this reason the parameters TRANSA and TRANSB
 are expressed in terms of the Sylvester equation (2).

 \endverbatim
 \remark This function is a wrapper for \ref sla_ggcsylv_dual.

 \see sla_ggcsylv_dual

 \param[in] FACTA
 \verbatim
          FACTA is String
          Specifies how the matrix pencil (A,C) is given.
          == 'N':  The matrix pencil (A,C) is given as a general matrices and its Schur decomposition
                  A = QA*S*ZA**T, C = QA*R*ZA**T will be computed.
          == 'F':  The matrix pencil (A,C) is already in generalized Schur form and S, R, QA, and ZA
                  are given.
 \endverbatim

 \param[in] FACTB
 \verbatim
          FACTB is String
          Specifies how the matrix pencil (B,D) is given.
          == 'N':  The matrix pencil (B,D) is given as a general matrices and its Schur decomposition
                  B = QB*U*ZB**T, D = QB*V*ZB**T will be computed.
          == 'F':  The matrix pencil (B,D) is already in generalized Schur form and U, V, QB, and ZB
                  are given.
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

 \param[in] SGN1
 \verbatim
          SGN1 is SINGLE PRECISION, allowed values: +/-1
          Specifies the sign between in the first equation.
 \endverbatim

 \param[in] SGN2
 \verbatim
          SGN2 is SINGLE PRECISION, allowed values: +/-1
          Specifies the sign between in the second equation.
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
 \endverbatim

 \param[in] LDZB
 \verbatim
          LDZB is INTEGER
          The leading dimension of the array ZB.  LDZB >= max(1,M).
 \endverbatim

 \param[in,out] E
 \verbatim
          E is SINGLE PRECISION array, dimension (LDE,N)
          On input, the matrix E contains the right hand side E.
          On output, the matrix E contains the solution R.
 \endverbatim

 \param[in] LDE
 \verbatim
          LDE is INTEGER
          The leading dimension of the array X.  LDE >= max(1,M).
 \endverbatim

 \param[in,out] F
 \verbatim
          F is SINGLE PRECISION array, dimension (LDF,N)
          On input, the matrix F contains the right hand side F.
          On output, the matrix F contains the solution L.
 \endverbatim

 \param[in] LDF
 \verbatim
          LDF is INTEGER
          The leading dimension of the array X.  LDF >= max(1,M).
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
          = 1:  DGGES failed
          = 2:  DLA_SORT_GEV failed
          = 3:  Inner solver failed
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date June 2023
*/


void mepack_single_ggcsylv_dual(const char *FACTA, const char *FACTB, const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
        float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
        float *QA, int LDQA, float *ZA, int LDZA, float *QB, int LDQB, float *ZB, int LDZB,
        float *E, int LDE, float *F, int LDF,
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
    Int _LDE  = LDE;
    Int _LDF  = LDF;
    Int _LDWORK = LDWORK;
    Int _INFO = *INFO;
    FORTRAN_LEN_T c1len = strlen(FACTA);
    FORTRAN_LEN_T c2len = strlen(FACTB);
    FORTRAN_LEN_T c3len = strlen(TRANSA);
    FORTRAN_LEN_T c4len = strlen(TRANSB);

    FC_GLOBAL_(sla_ggcsylv_dual,SLA_GGCSYLV_DUAL)( FACTA, FACTB, TRANSA, TRANSB, &SGN1, &SGN2, &_M, &_N,
            A, &_LDA, B, &_LDB, C, &_LDC, D, &_LDD,
            QA, &_LDQA, ZA, &_LDZA, QB, &_LDQB, ZB, &_LDZB, E, &_LDE, F, &_LDF,
            SCALE, WORK, &_LDWORK, &_INFO,
            c1len, c2len, c3len, c4len);
    *INFO = (int) _INFO;
    return;

}






/**
 * @}
 */

