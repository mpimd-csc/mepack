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
 * @addtogroup ctrsylv
 * @{
 */
/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation.

 \par Purpose:

 \verbatim
 mepack_double_trsylv2_level2 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via DGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref dla_trsylv2_l2.

 \see dla_trsylv2_l2

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/
void mepack_double_trsylv2_level2(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(dla_trsylv2_l2,DLA_TRSYLV2_L2)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}


/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Reordered variant)

 \par Purpose:

 \verbatim
 mepack_double_trsylv2_level2_reorder solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via DGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref dla_trsylv2_l2_reorder.

 \see dla_trsylv2_l2_reorder

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/
void mepack_double_trsylv2_level2_reorder(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(dla_trsylv2_l2_reorder,DLA_TRSYLV2_L2_REORDER)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Unoptimized variant)

 \par Purpose:

 \verbatim
 mepack_double_trsylv2_level2_unopt solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via DGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref dla_trsylv2_l2_unopt.

 \see dla_trsylv2_l2_unopt

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/
void mepack_double_trsylv2_level2_unopt(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(dla_trsylv2_l2_unopt,DLA_TRSYLV2_L2_UNOPT)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Local copy variant)

 \par Purpose:

 \verbatim
 mepack_double_trsylv2_level2_local_copy solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via DGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref dla_trsylv2_l2_local_copy.

 \see dla_trsylv2_l2_local_copy

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/
void mepack_double_trsylv2_level2_local_copy(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(dla_trsylv2_l2_local_copy,DLA_TRSYLV2_L2_LOCAL_COPY)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Local copy variant, M,N <= 32)

 \par Purpose:

 \verbatim
 mepack_double_trsylv2_level2_local_copy_32 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via DGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref dla_trsylv2_l2_local_copy_32.

 \see dla_trsylv2_l2_local_copy_32

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/

void mepack_double_trsylv2_level2_local_copy_32(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(dla_trsylv2_l2_local_copy_32,DLA_TRSYLV2_L2_LOCAL_COPY_32)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Local copy variant, M,N <= 64)

 \par Purpose:

 \verbatim
 mepack_double_trsylv2_level2_local_copy_64 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via DGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref dla_trsylv2_l2_local_copy_64.

 \see dla_trsylv2_l2_local_copy_64

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/

void mepack_double_trsylv2_level2_local_copy_64(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(dla_trsylv2_l2_local_copy_64,DLA_TRSYLV2_L2_LOCAL_COPY_64)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Local copy variant, M,N<=96)

 \par Purpose:

 \verbatim
 mepack_double_trsylv2_level2_local_copy_96 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via DGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref dla_trsylv2_l2_local_copy_96.

 \see dla_trsylv2_l2_local_copy_96

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/

void mepack_double_trsylv2_level2_local_copy_96(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(dla_trsylv2_l2_local_copy_96,DLA_TRSYLV2_L2_LOCAL_COPY_96)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Local copy variant, M,N <= 128)

 \par Purpose:

 \verbatim
 mepack_double_trsylv2_level2_local_copy_128 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via DGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref dla_trsylv2_l2_local_copy_128.

 \see dla_trsylv2_l2_local_copy_128

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/

void mepack_double_trsylv2_level2_local_copy_128(const char *TRANSA, const char *TRANSB, double SGN,  int M, int N,
                                 double * A, int LDA, double * B, int LDB, double *X, int LDX,
                                 double * SCALE, double *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(dla_trsylv2_l2_local_copy_128,DLA_TRSYLV2_L2_LOCAL_COPY_128)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation.

 \par Purpose:

 \verbatim
 mepack_single_trsylv2_level2 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via SGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref sla_trsylv2_l2.

 \see sla_trsylv2_l2

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
          Specifies the form of the system of equations with respect to B:
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is SINGLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/


void mepack_single_trsylv2_level2(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(sla_trsylv2_l2,SLA_TRSYLV2_L2)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Reordered variant)

 \par Purpose:

 \verbatim
 mepack_single_trsylv2_level2_reorder solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via SGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref sla_trsylv2_l2_reorder.

 \see sla_trsylv2_l2_reorder

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
          Specifies the form of the system of equations with respect to B:
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is SINGLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/


void mepack_single_trsylv2_level2_reorder(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(sla_trsylv2_l2_reorder,SLA_TRSYLV2_L2_REORDER)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Unoptimized variant)

 \par Purpose:

 \verbatim
 mepack_single_trsylv2_level2_unopt solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via SGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref sla_trsylv2_l2_unopt.

 \see sla_trsylv2_l2_unopt

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
          Specifies the form of the system of equations with respect to B:
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is SINGLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/
void mepack_single_trsylv2_level2_unopt(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(sla_trsylv2_l2_unopt,SLA_TRSYLV2_L2_UNOPT)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Local copy variant)

 \par Purpose:

 \verbatim
 mepack_single_trsylv2_level2_local_copy solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via SGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref sla_trsylv2_l2_local_copy.

 \see sla_trsylv2_l2_local_copy

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
          Specifies the form of the system of equations with respect to B:
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is SINGLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/


void mepack_single_trsylv2_level2_local_copy(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(sla_trsylv2_l2_local_copy,SLA_TRSYLV2_L2_LOCAL_COPY)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Local copy variant, M,N <= 32 )

 \par Purpose:

 \verbatim
 mepack_single_trsylv2_level2_local_copy_32 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via SGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref sla_trsylv2_l2_local_copy_32.

 \see sla_trsylv2_l2_local_copy_32

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
          Specifies the form of the system of equations with respect to B:
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is SINGLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/
void mepack_single_trsylv2_level2_local_copy_32(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(sla_trsylv2_l2_local_copy_32,SLA_TRSYLV2_L2_LOCAL_COPY_32)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Local copy variant, M,N <= 64)

 \par Purpose:

 \verbatim
 mepack_double_trsylv2_level2_local_copy_64 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via DGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref dla_trsylv2_l2_local_copy_64.

 \see dla_trsylv2_l2_local_copy_64

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/
void mepack_single_trsylv2_level2_local_copy_64(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(sla_trsylv2_l2_local_copy_64,SLA_TRSYLV2_L2_LOCAL_COPY_64)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Local copy variant, M,N <= 96)

 \par Purpose:

 \verbatim
 mepack_single_trsylv2_level2_local_copy_96 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via SGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref sla_trsylv2_l2_local_copy_96.

 \see sla_trsylv2_l2_local_copy_96

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
          Specifies the form of the system of equations with respect to B:
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is SINGLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/


void mepack_single_trsylv2_level2_local_copy_96(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(sla_trsylv2_l2_local_copy_96,SLA_TRSYLV2_L2_LOCAL_COPY_96)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation. (Local copy variant, M,N <= 128)

 \par Purpose:

 \verbatim
 mepack_single_trsylv2_level2_local_copy_128 solves a Sylvester equation of the following forms

    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)

 or

    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
 matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
 A and B are generated via SGEES form LAPACK.
 \endverbatim

 \remark This function is a wrapper around \ref sla_trsylv2_l2_local_copy_128.

 \see sla_trsylv2_l2_local_copy_128

 \param[in] TRANSA
 \verbatim
          TRANSA is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is string
          Specifies the form of the system of equations with respect to B:
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
          The matrix B must be (quasi-) upper triangular.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
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
          The leading dimension of the array X.  LDB >= max(1,M).
 \endverbatim

 \param[out] SCALE
 \verbatim
          SCALE is SINGLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
          could not be solved correctly, 0 < SCALE <= 1 holds true.
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

 \date January 2024
*/


void mepack_single_trsylv2_level2_local_copy_128(const char *TRANSA, const char *TRANSB, float SGN,  int M, int N,
                                 float * A, int LDA, float * B, int LDB, float *X, int LDX,
                                 float * SCALE, float *WORK, int *INFO)
{
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDX = LDX;
    Int _M   = M;
    Int _N   = N;
    Int _INFO = (Int) *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(sla_trsylv2_l2_local_copy_128,SLA_TRSYLV2_L2_LOCAL_COPY_128)(TRANSA, TRANSB, &SGN, &_M, &_N,
            A, &_LDA, B, &_LDB, X, &_LDX, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO;
    return;
}

/**
 * @}
 */

