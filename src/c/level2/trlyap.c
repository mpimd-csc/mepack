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
 * Copyright (C) Martin Koehler, 2017-2022
 */

#include "FCMangle.h"
#include "mepack.h"
#include "mepack_internal.h"
#include <string.h>

/**
 * @addtogroup ctrlyap
 * @{
 */
/**

 \brief Level-2 Bartels-Stewart Algorithm for the Lyapunov Equation

 \par Purpose:

 \verbatim

 mepack_double_trlyap_level2 solves a Lyapunov equation of the following forms

    A * X  +  X * A**T = SCALE * Y                              (1)

 or

    A **T * X  -  X * A = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix.
 The right hand side Y and the solution X are M-by-N matrices.  Typically the matrix A
 is generated by DGEES from LAPACK.
 \endverbatim

 \remark The algorithm is implemented using BLAS level 2 operations.
 \remark This is a wrapper around \ref dla_trlyap_l2
 \see dla_trlyap_l2

 \param[in] TRANS
 \verbatim
          TRANS is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrices A and C.  M >= 0.
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

 \date Dezember 2022
*/
void mepack_double_trlyap_level2(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO)
{
    Int _M = M;
    Int _LDA = LDA;
    Int _LDX = LDX;
    Int _INFO = *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANS);
    FC_GLOBAL_(dla_trlyap_l2,DLA_TRLYAP_L2)(TRANS, &_M, A, &_LDA, X, &_LDX, SCALE, WORK, &_INFO, c1len);
    *INFO = (int) _INFO;
    return;
}

/**

 \brief Level-2 Bartels-Stewart Algorithm for the Lyapunov Equation (Optimized for small problems)

 \par Purpose:

 \verbatim

 mepack_double_trlyap_level2 solves a Lyapunov equation of the following forms

    A * X  +  X * A**T = SCALE * Y                              (1)

 or

    A **T * X  -  X * A = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix.
 The right hand side Y and the solution X are M-by-N matrices.  Typically the matrix A
 is generated by DGEES from LAPACK.
 \endverbatim

 \remark The algorithm is implemented using BLAS level 2 operations.
 \remark This is a wrapper around \ref dla_trlyap_l2_opt
 \see dla_trlyap_l2_opt

 \param[in] TRANS
 \verbatim
          TRANS is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrices A and C.  M >= 0.
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

 \date Dezember 2022
*/
void mepack_double_trlyap_level2_opt(const char *TRANS, int M, double * A, int LDA, double *X, int LDX, double * SCALE, double *WORK, int *INFO)
{
    Int _M = M;
    Int _LDA = LDA;
    Int _LDX = LDX;
    Int _INFO = *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANS);
    FC_GLOBAL_(dla_trlyap_l2_opt,DLA_TRLYAP_L2_OPT)(TRANS, &_M, A, &_LDA, X, &_LDX, SCALE, WORK, &_INFO, c1len);
    *INFO = (int) _INFO;
    return;
}

/**
 \brief Level-2 Bartels-Stewart Algorithm for the Lyapunov Equation

 \par Purpose:

 \verbatim

 mepack_single_trlyap_level2 solves a Lyapunov equation of the following forms

    A * X  +  X * A**T = SCALE * Y                              (1)

 or

    A **T * X  -  X * A = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix.
 The right hand side Y and the solution X are M-by-N matrices.  Typically the matrix A
 is generated by SGEES from LAPACK.
 \endverbatim

 \remark The algorithm is implemented using BLAS level 2 operations.
 \remark This is a wrapper around \ref sla_trlyap_l2
 \see sla_trlyap_l2

 \param[in] TRANS
 \verbatim
          TRANS is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrices A and C.  M >= 0.
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

 \date Dezember 2022
*/
void mepack_single_trlyap_level2(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO)
{
    Int _M = M;
    Int _LDA = LDA;
    Int _LDX = LDX;
    Int _INFO = *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANS);
    FC_GLOBAL_(sla_trlyap_l2,SLA_TRLYAP_L2)(TRANS, &_M, A, &_LDA, X, &_LDX, SCALE, WORK, &_INFO, c1len);
    *INFO = (int) _INFO;
    return;
}


/**
 \brief Level-2 Bartels-Stewart Algorithm for the Lyapunov Equation (Optimized for small problems)

 \par Purpose:

 \verbatim

 mepack_single_trlyap_level2_opt solves a Lyapunov equation of the following forms

    A * X  +  X * A**T = SCALE * Y                              (1)

 or

    A **T * X  -  X * A = SCALE * Y                              (2)

 where A is a M-by-M quasi upper triangular matrix.
 The right hand side Y and the solution X are M-by-N matrices.  Typically the matrix A
 is generated by SGEES from LAPACK.
 \endverbatim

 \remark The algorithm is implemented using BLAS level 2 operations.
 \remark This is a wrapper around \ref sla_trlyap_l2_opt
 \see sla_trlyap_l2_opt

 \param[in] TRANS
 \verbatim
          TRANS is string
          Specifies the form of the system of equations with respect to A:
          == 'N':  op1(A) = A
          == 'T':  op1(A) = A**T
 \endverbatim

 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrices A and C.  M >= 0.
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

 \date Dezember 2022
*/
void mepack_single_trlyap_level2_opt(const char *TRANS, int M, float * A, int LDA, float *X, int LDX, float * SCALE, float *WORK, int *INFO)
{
    Int _M = M;
    Int _LDA = LDA;
    Int _LDX = LDX;
    Int _INFO = *INFO;
    FORTRAN_LEN_T c1len = strlen(TRANS);
    FC_GLOBAL_(sla_trlyap_l2_opt,SLA_TRLYAP_L2_OPT)(TRANS, &_M, A, &_LDA, X, &_LDX, SCALE, WORK, &_INFO, c1len);
    *INFO = (int) _INFO;
    return;
}


/**
 * @}
 */

