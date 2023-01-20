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
 * @addtogroup ctgsylv
 * @{
 */
/**
 \brief Level-3 Bartels-Stewart Algorithm for the generalized Sylvester equation with DAG based parallelization.

 \par Purpose:

 \verbatim

 mepack_double_tgsylv_dag solves a generalized Sylvester equation of the following forms

    op1(A) * X * op2(B) + op1(C) * X * op2(D) = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) - op1(C) * X * op2(D) = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
 matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
 D can be quasi upper triangular as well. The right hand side Y and the solution X
 M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
 created by DGGES from LAPACK.
 \endverbatim

 \remark This function is a wrapper for \ref dla_tgsylv_dag.

 \see dla_tgsylv_dag

 \param[in] TRANSA
 \verbatim
          TRANSA is String
          Specifies the form of the system of equations with respect to A and C:
          == 'N':  op1(A) = A, op1(C) = C (No transpose for A and C)
          == 'T':  op1(A) = A**T, op1(C) = C **T (Transpose A and C)
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is String
          Specifies the form of the system of equations with respect to B and D:
          == 'N':  op2(B) = B, op2(D) = D (No transpose for B and D)
          == 'T':  op2(B) = B**T, op2(D) = D **T (Transpose B and D)
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

 \param[in] A
 \verbatim
          A is DOUBLE PRECISION array, dimension (LDA,M)
          The matrix A must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim

 \param[in] B
 \verbatim
          B is DOUBLE PRECISION array, dimension (LDB,N)
          The matrix B must be (quasi-) upper triangular. If the matrix D is already
          quasi-upper triangular the matrix B must be upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
 \endverbatim

 \param[in] C
 \verbatim
          C is DOUBLE PRECISION array, dimension (LDC,M)
          The matrix C must be upper triangular.
 \endverbatim

 \param[in] LDC
 \verbatim
          LDC is INTEGER
          The leading dimension of the array C.  LDC >= max(1,M).
 \endverbatim

 \param[in] D
 \verbatim
          D is DOUBLE PRECISION array, dimension (LDD,N)
          The matrix D must be (quasi-) upper triangular. If the matrix B is already
          quasi-upper triangular the matrix D must be upper triangular.
 \endverbatim

 \param[in] LDD
 \verbatim
          LDD is INTEGER
          The leading dimension of the array D.  LDD >= max(1,N).
 \endverbatim

 \param[in,out] X
 \verbatim
          X is DOUBLE PRECISION array, dimension (LDX,N)
          On input, the matrix X contains the right hand side Y.
          On output, the matrix X contains the solution of Equation (1) or (2)
          as selected by TRANSA, TRANSB, and SGN.
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
          If SCALE .NE. 1 the problem is no solved correctly in this case
          one have to use an other solver.
 \endverbatim

 \param[in] WORK
 \verbatim
          WORK is DOUBLE PRECISION array, dimension LDWORK
          Workspace for the algorithm.
          The workspace needs to queried before the running the computation.
          The query is performed by calling the subroutine with INFO == -1 on input.
          The required workspace is then returned in INFO.
 \endverbatim

 \param[inout] INFO
 \verbatim
          INFO is INTEGER

          On input:
            == -1 : Perform a workspace query
            <> -1 : normal operation

          On exit, workspace query:
            < 0 :  if INFO == -i, the i-Th argument had an illegal value
            >= 0:  The value of INFO is the required number of elements in the workspace.

          On exit, normal operation:
            == 0:  successful exit
            < 0:  if INFO == -i, the i-Th argument had an illegal value
            > 0:  The equation is not solved correctly. One of the arising inner
                  system got singular.
 \endverbatim

 \author Martin Koehler, MPI Magdeburg

 \date Januar 2023
*/
void mepack_double_tgsylv_dag(const char *TRANSA, const char*TRANSB, double SGN,
        int M, int N, double * A, int LDA, double * B, int LDB,
        double *C, int LDC, double *D, int LDD, double *X, int LDX,
        double * SCALE, double *WORK, int *INFO)
{
    Int _LDA  = LDA;
    Int _LDB  = LDB;
    Int _LDC  = LDC;
    Int _LDD  = LDD;
    Int _LDX  = LDX;
    Int _INFO = *INFO;
    Int _M    = M;
    Int _N    = N;

    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(dla_tgsylv_dag,DLA_TGSYLV_DAG)(TRANSA, TRANSB, &SGN, &_M, &_N, A, &_LDA, B, &_LDB,
            C, &_LDC, D, &_LDD, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);

    *INFO = (int) _INFO ;
    return;

}

/**
 \brief Level-3 Bartels-Stewart Algorithm for the generalized Sylvester equation with DAG based parallelization.

 \par Purpose:

 \verbatim

 mepack_single_tgsylv_dag solves a generalized Sylvester equation of the following forms

    op1(A) * X * op2(B) + op1(C) * X * op2(D) = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) - op1(C) * X * op2(D) = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
 matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
 D can be quasi upper triangular as well. The right hand side Y and the solution X
 M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
 created by SGGES from LAPACK.
 \endverbatim

 \remark This function is a wrapper for \ref sla_tgsylv_dag.

 \see sla_tgsylv_dag

 \param[in] TRANSA
 \verbatim
          TRANSA is String
          Specifies the form of the system of equations with respect to A and C:
          == 'N':  op1(A) = A, op1(C) = C (No transpose for A and C)
          == 'T':  op1(A) = A**T, op1(C) = C **T (Transpose A and C)
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is String
          Specifies the form of the system of equations with respect to B and D:
          == 'N':  op2(B) = B, op2(D) = D (No transpose for B and D)
          == 'T':  op2(B) = B**T, op2(D) = D **T (Transpose B and D)
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

 \param[in] A
 \verbatim
          A is SINGLE PRECISION array, dimension (LDA,M)
          The matrix A must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim

 \param[in] B
 \verbatim
          B is SINGLE PRECISION array, dimension (LDB,N)
          The matrix B must be (quasi-) upper triangular. If the matrix D is already
          quasi-upper triangular the matrix B must be upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
 \endverbatim

 \param[in] C
 \verbatim
          C is SINGLE PRECISION array, dimension (LDC,M)
          The matrix C must be upper triangular.
 \endverbatim

 \param[in] LDC
 \verbatim
          LDC is INTEGER
          The leading dimension of the array C.  LDC >= max(1,M).
 \endverbatim

 \param[in] D
 \verbatim
          D is SINGLE PRECISION array, dimension (LDD,N)
          The matrix D must be (quasi-) upper triangular. If the matrix B is already
          quasi-upper triangular the matrix D must be upper triangular.
 \endverbatim

 \param[in] LDD
 \verbatim
          LDD is INTEGER
          The leading dimension of the array D.  LDD >= max(1,N).
 \endverbatim

 \param[in,out] X
 \verbatim
          X is SINGLE PRECISION array, dimension (LDX,N)
          On input, the matrix X contains the right hand side Y.
          On output, the matrix X contains the solution of Equation (1) or (2)
          as selected by TRANSA, TRANSB, and SGN.
          Right hand side Y and the solution X are M-by-N matrices.
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
          If SCALE .NE. 1 the problem is no solved correctly in this case
          one have to use an other solver.
 \endverbatim

 \param[in] WORK
 \verbatim
          WORK is SINGLE PRECISION array, dimension LDWORK
          Workspace for the algorithm.
          The workspace needs to queried before the running the computation.
          The query is performed by calling the subroutine with INFO == -1 on input.
          The required workspace is then returned in INFO.
 \endverbatim

 \param[inout] INFO
 \verbatim
          INFO is INTEGER

          On input:
            == -1 : Perform a workspace query
            <> -1 : normal operation

          On exit, workspace query:
            < 0 :  if INFO == -i, the i-Th argument had an illegal value
            >= 0:  The value of INFO is the required number of elements in the workspace.

          On exit, normal operation:
            == 0:  successful exit
            < 0:  if INFO == -i, the i-Th argument had an illegal value
            > 0:  The equation is not solved correctly. One of the arising inner
                  system got singular.
 \endverbatim

 \author Martin Koehler, MPI Magdeburg

 \date Januar 2023
*/
void mepack_single_tgsylv_dag(const char *TRANSA, const char*TRANSB, float SGN,
        int M, int N, float * A, int LDA, float * B, int LDB,
        float *C, int LDC, float *D, int LDD, float *X, int LDX,
        float * SCALE, float *WORK, int *INFO)
{
    Int _LDA  = LDA;
    Int _LDB  = LDB;
    Int _LDC  = LDC;
    Int _LDD  = LDD;
    Int _LDX  = LDX;
    Int _INFO = *INFO;
    Int _M    = M;
    Int _N    = N;

    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(sla_tgsylv_dag,SLA_TGSYLV_DAG)(TRANSA, TRANSB, &SGN, &_M, &_N, A, &_LDA, B, &_LDB,
            C, &_LDC, D, &_LDD, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);

    *INFO = (int) _INFO ;
    return;

}

/**
 * @}
 */

