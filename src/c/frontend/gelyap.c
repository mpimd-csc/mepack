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
 * @addtogroup cgelyap
 * @{
 */

/**
 \brief Frontend for the solution of Standard Lyapunov Equations.

 \par Purpose:

 \verbatim
 mepack_double_gelyap solves a Lyapunov equation of the following forms

    A  * X  +  X * A**T = SCALE * Y                              (1)

 or

    A ** T * X  +  X * A = SCALE * Y                              (2)

 where A is a M-by-M general matrix or a matrix in upper Hessenberg form.
 The right hand side Y and the solution X are M-by-M matrices.  
 The general matrix A can supplied factorized in terms of its
 Schur decomposition.

 \endverbatim

 \remark This function is a wrapper around \ref dla_gelyap.

 \see dla_gelyap

 \param[in] FACT
 \verbatim
          FACT is String
          Specifies how the matrix A is given.
          == 'N':  The matrix A is given as a general matrix and its Schur decomposition
                  A = Q*S*Q**T will be computed.
          == 'F':  The matrix A is given as its Schur decomposition in terms of S and Q
                  form A = Q*S*Q**T
          == 'H':  The matrix A is given in upper Hessenberg form and its Schur decomposition
                  A = Q*S*Q**T will be computed
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
          The order of the matrix A.  M >= 0.
 \endverbatim

 \param[in,out] A
 \verbatim
          A is DOUBLE PRECISION array, dimension (LDA,M)
          If FACT == "N", the matrix A is a general matrix and it is overwritten with its
          schur decomposition S.
          If FACT == "F", the matrix A contains its (quasi-) upper triangular matrix S being the
          Schur decomposition of A.
          If FACT == "H", the matrix A is an upper Hessenberg matrix and it is overwritten
          with its schur decomposition S.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim

 \param[in,out] Q
 \verbatim
          Q is DOUBLE PRECISION array, dimension (LDQ,M)
          If FACT == "N", the matrix Q is an empty M-by-M matrix on input and contains the
          Schur vectors of A on output.
          If FACT == "F", the matrix Q contains the Schur vectors of A.
          If FACT == "H", the matrix Q is an empty M-by-M matrix on input and contains the
          Schur vectors of A on output.
 \endverbatim

 \param[in] LDQ
 \verbatim
          LDQ is INTEGER
          The leading dimension of the array Q.  LDQ >= max(1,M).
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
          = 1:  DGEES failed
          = 2:  DLA_SORT_EV failed
          = 3:  Internal solver failed
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date October 2023
 */
void mepack_double_gelyap(const char * FACT, const char *TRANS, int M, double * A,  int LDA, double * Q, int LDQ, double *X, int LDX,
        double * SCALE, double *WORK, size_t LDWORK, int *INFO)
{
    Int _M = M;
    Int _LDA = LDA;
    Int _LDQ = LDQ;
    Int _LDX = LDX;
    Int _INFO = * INFO;
    Int _LDWORK = LDWORK;
    FORTRAN_LEN_T c1len = strlen(FACT);
    FORTRAN_LEN_T c2len = strlen(TRANS);

    FC_GLOBAL_(dla_gelyap,DLA_GELYAP)(FACT, TRANS, &_M, A, &_LDA, Q, &_LDQ, X, &_LDX, SCALE, WORK, &_LDWORK, &_INFO, c1len, c2len);

    *INFO = _INFO;
    return;

}


/**
 \brief Frontend for the solution of Standard Lyapunov Equations.

 \par Purpose:

 \verbatim
 mepack_single_gelyap solves a Lyapunov equation of the following forms

    A  * X  +  X * A**T = SCALE * Y                              (1)

 or

    A ** T * X  +  X * A = SCALE * Y                              (2)

 where A is a M-by-M general matrix or a matrix in upper Hessenberg form.
 The right hand side Y and the solution X are M-by-M matrices.
 The general matrix A can supplied factorized in terms of its
 Schur decomposition.

 \endverbatim

 \remark This function is a wrapper around \ref sla_gelyap.

 \see sla_gelyap

 \param[in] FACT
 \verbatim
          FACT is String
          Specifies how the matrix A is given.
          == 'N':  The matrix A is given as a general matrix and its Schur decomposition
                  A = Q*S*Q**T will be computed.
          == 'F':  The matrix A is given as its Schur decomposition in terms of S and Q
                  form A = Q*S*Q**T
          == 'H':  The matrix A is given in upper Hessenberg form and its Schur decomposition
                  A = Q*S*Q**T will be computed
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
          The order of the matrix A.  M >= 0.
 \endverbatim

 \param[in,out] A
 \verbatim
          A is SINGLE PRECISION array, dimension (LDA,M)
          If FACT == "N", the matrix A is a general matrix and it is overwritten with its
          schur decomposition S.
          If FACT == "F", the matrix A contains its (quasi-) upper triangular matrix S being the
          Schur decomposition of A.
          If FACT == "H", the matrix A is an upper Hessenberg matrix and it is overwritten
	  with its schur decomposition S.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim

 \param[in,out] Q
 \verbatim
          Q is SINGLE PRECISION array, dimension (LDQ,M)
          If FACT == "N", the matrix Q is an empty M-by-M matrix on input and contains the
          Schur vectors of A on output.
          If FACT == "F", the matrix Q contains the Schur vectors of A.
          If FACT == "H", the matrix Q is an empty M-by-M matrix on input and contains the
          Schur vectors of A on output.
 \endverbatim

 \param[in] LDQ
 \verbatim
          LDQ is INTEGER
          The leading dimension of the array Q.  LDQ >= max(1,M).
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
          Workspace for the algorithm. The optimal workspace is given either by \ref mepack_memory_frontend
          or a previous call to the this routine with LDWORK === -1.
          On exit, WORK(1) contains the size of the workspace.
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
          = 3:  Internal solver failed
          < 0:  if INFO == -i, the i-Th argument had an illegal value
 \endverbatim

 \attention The Fortran/LAPACK-like workspace query with setting LDWORK=-1 on input will not work in the C interface.
            One have to use the \ref mepack_memory_frontend function for this purpose.

 \author Martin Koehler, MPI Magdeburg

 \date October 2023
 */
void mepack_single_gelyap(const char * FACT, const char *TRANS, int M, float * A,  int LDA, float * Q, int LDQ, float *X, int LDX,
        float * SCALE, float *WORK, size_t LDWORK, int *INFO)
{
    Int _M = M;
    Int _LDA = LDA;
    Int _LDQ = LDQ;
    Int _LDX = LDX;
    Int _INFO = * INFO;
    Int _LDWORK = LDWORK;
    FORTRAN_LEN_T c1len = strlen(FACT);
    FORTRAN_LEN_T c2len = strlen(TRANS);

    FC_GLOBAL_(sla_gelyap,SLA_GELYAP)(FACT, TRANS, &_M, A, &_LDA, Q, &_LDQ, X, &_LDX, SCALE, WORK, &_LDWORK, &_INFO, c1len, c2len);

    *INFO = _INFO;
    return;

}




/**
 * @}
 */

