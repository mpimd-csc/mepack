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
 * @addtogroup cresidual
 * @{
 */
/**
 \brief Compute the relative residual for the generalized Sylvester equation.

 \par Purpose:

 \verbatim

 Compute the relative residual of the generalized Sylvester equation

     RelRes = || SCALE*Y - opA(A)*X *opB(B) - SGN * opA(C) * X * opB(D) || / ( (||A||*||B||+||C||*||D||)*||X|| + SCALE * ||Y|| )     (1)

 The norms are evaluated in terms of the Frobenius norm.
 \endverbatim
 \remark Auxiliary memory is managed by the function itself.

 \param[in] TRANSA
 \verbatim
          TRANSA is CHARACTER(1)
          Specifies the form of the system of equations with respect to A :
          == 'N':  opA(A) = A (No transpose for A)
          == 'T':  opA(A) = A**T (Transpose A)
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is CHARACTER(1)
          Specifies the form of the system of equations with respect to B :
          == 'N':  opB(B) = B (No transpose for B)
          == 'T':  opB(B) = B**T (Transpose B)
 \endverbatim

 \param[in] SGN
 \verbatim
          SGN is DOUBLE PRECISION, allowed values: +/-1
          Specifies the sign between the two parts of the Sylvester equation.
 \endverbatim

 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrix A.  M >= 0.
 \endverbatim

 \param[in] N
 \verbatim
          N is INTEGER
          The order of the matrix B.  N >= 0.
 \endverbatim

 \param[in] A
 \verbatim
          A is DOUBLE PRECISION array, dimension (LDA,M)
          The coefficient matrix A.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim

 \param[in] B
 \verbatim
          B is DOUBLE PRECISION array, dimension (LDB,N)
          The coefficient matrix B.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
 \endverbatim

 \param[in] C
 \verbatim
          C is DOUBLE PRECISION array, dimension (LDC,M)
          The coefficient matrix C.
 \endverbatim

 \param[in] LDC
 \verbatim
          LDC is INTEGER
          The leading dimension of the array C.  LDC >= max(1,M).
 \endverbatim

 \param[in] D
 \verbatim
          D is DOUBLE PRECISION array, dimension (LDD,N)
          The coefficient matrix D.
 \endverbatim

 \param[in] LDD
 \verbatim
          LDD is INTEGER
          The leading dimension of the array D.  LDD >= max(1,N).
 \endverbatim

 \param[in] X
 \verbatim
          X is DOUBLE PRECISION array, dimension (LDX,N)
          Solution of the Sylvester Equation.
 \endverbatim

 \param[in] LDX
 \verbatim
          LDX is INTEGER
          The leading dimension of the array X.  LDX >= max(1,M).
 \endverbatim

 \param[in] Y
 \verbatim
          Y is DOUBLE PRECISION array, dimension (LDY,N)
          The right hand side of the Sylvester Equation.
 \endverbatim

 \param[in] LDY
 \verbatim
          LDY is INTEGER
          The leading dimension of the array Y.  LDY >= max(1,M).
 \endverbatim

 \param[in] SCALE
 \verbatim
          SCALE is DOUBLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
 \endverbatim

 \return
 \verbatim
          The function returns the relative residual as defined in (1). If an error
          occurs a negative integer I is returned, signalizing that the parameter -I
          has a wrong/illegal value. If I = -1000, the auxiliary memory could not be allocated.
 \endverbatim

 \see dla_residual_gsylv

 \author Martin Koehler, MPI Magdeburg

 \date January 2024


 */

double mepack_double_residual_gsylv(const char * TRANSA, const char *TRANSB, double sgn, int M, int N, double *A, int LDA,
        double *B, int LDB, double *C, int LDC, double *D, int LDD,
        double *X, int LDX, double *Y, int LDY, double SCALE)
{
    Int _M   = M;
    Int _N   = N;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDC = LDC;
    Int _LDD = LDD;
    Int _LDX = LDX;
    Int _LDY = LDY;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    return FC_GLOBAL_(dla_residual_gsylv, DLA_RESIDUAL_GSYLV) (TRANSA, TRANSB, &sgn, &_M, &_N, A, &_LDA, B, &_LDB, C, &_LDC, D, &_LDD, X, &_LDX, Y, &_LDY, &SCALE, c1len, c2len);
}

/**
 \brief Compute the relative residual for the generalized Sylvester equation.

 \par Purpose:

 \verbatim

 Compute the relative residual of the generalized Sylvester equation

     RelRes = || SCALE*Y - opA(A)*X *opB(B) - SGN * opA(C) * X * opB(D) || / ( (||A||*||B||+||C||*||D||)*||X|| + SCALE * ||Y|| )     (1)

 The norms are evaluated in terms of the Frobenius norm.
 \endverbatim
 \remark Auxiliary memory is managed by the function itself.

 \param[in] TRANSA
 \verbatim
          TRANSA is CHARACTER(1)
          Specifies the form of the system of equations with respect to A :
          == 'N':  opA(A) = A (No transpose for A)
          == 'T':  opA(A) = A**T (Transpose A)
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is CHARACTER(1)
          Specifies the form of the system of equations with respect to B :
          == 'N':  opB(B) = B (No transpose for B)
          == 'T':  opB(B) = B**T (Transpose B)
 \endverbatim

 \param[in] SGN
 \verbatim
          SGN is SINGLE PRECISION, allowed values: +/-1
          Specifies the sign between the two parts of the Sylvester equation.
 \endverbatim

 \param[in] M
 \verbatim
          M is INTEGER
          The order of the matrix A.  M >= 0.
 \endverbatim

 \param[in] N
 \verbatim
          N is INTEGER
          The order of the matrix B.  N >= 0.
 \endverbatim

 \param[in] A
 \verbatim
          A is SINGLE PRECISION array, dimension (LDA,M)
          The coefficient matrix A.
 \endverbatim

 \param[in] LDA
 \verbatim
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
 \endverbatim

 \param[in] B
 \verbatim
          B is SINGLE PRECISION array, dimension (LDB,N)
          The coefficient matrix B.
 \endverbatim

 \param[in] LDB
 \verbatim
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
 \endverbatim

 \param[in] C
 \verbatim
          C is SINGLE PRECISION array, dimension (LDC,M)
          The coefficient matrix C.
 \endverbatim

 \param[in] LDC
 \verbatim
          LDC is INTEGER
          The leading dimension of the array C.  LDC >= max(1,M).
 \endverbatim

 \param[in] D
 \verbatim
          D is SINGLE PRECISION array, dimension (LDD,N)
          The coefficient matrix D.
 \endverbatim

 \param[in] LDD
 \verbatim
          LDD is INTEGER
          The leading dimension of the array D.  LDD >= max(1,N).
 \endverbatim

 \param[in] X
 \verbatim
          X is SINGLE PRECISION array, dimension (LDY,M)
          Solution of the Sylvester Equation.
 \endverbatim

 \param[in] LDX
 \verbatim
          LDX is INTEGER
          The leading dimension of the array X.  LDX >= max(1,M).
 \endverbatim

 \param[in] Y
 \verbatim
          Y is SINGLE PRECISION array, dimension (LDY,N)
          The right hand side of the Sylvester Equation.
 \endverbatim

 \param[in] LDY
 \verbatim
          LDY is INTEGER
          The leading dimension of the array Y.  LDY >= max(1,M).
 \endverbatim

 \param[in] SCALE
 \verbatim
          SCALE is SINGLE PRECISION
          SCALE is a scaling factor to prevent the overflow in the result.
 \endverbatim

 \return
 \verbatim
          The function returns the relative residual as defined in (1). If an error
          occurs a negative integer I is returned, signalizing that the parameter -I
          has a wrong/illegal value. If I = -1000, the auxiliary memory could not be allocated.
 \endverbatim

 \see sla_residual_gsylv

 \author Martin Koehler, MPI Magdeburg

 \date January 2024


 */

float mepack_single_residual_gsylv(const char * TRANSA, const char *TRANSB, float sgn, int M, int N, float *A, int LDA,
        float *B, int LDB, float *C, int LDC, float *D, int LDD,
        float *X, int LDX, float *Y, int LDY, float SCALE)
{
    Int _M   = M;
    Int _N   = N;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDC = LDC;
    Int _LDD = LDD;
    Int _LDX = LDX;
    Int _LDY = LDY;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    return FC_GLOBAL_(sla_residual_gsylv, SLA_RESIDUAL_GSYLV) (TRANSA, TRANSB, &sgn, &_M, &_N, A, &_LDA, B, &_LDB, C, &_LDC, D, &_LDD, X, &_LDX, Y, &_LDY, &SCALE, c1len, c2len);
}


/**
 * @}
 */

