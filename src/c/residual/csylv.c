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
 \brief Compute the relative residual for the coupled Sylvester equation.

 \par Purpose:

 \verbatim

 Compute the relative residual of the coupled Sylvester equation

 RelRes = max( || SCALE*E - opA(A)*R - SGN1 * L *opB(B) || / ( (||A||*||R||+||B||*||L||) + SCALE * ||E|| ),
               || SCALE*F - opA(C)*R - SGN2 * L *opB(D) || / ( (||C||*||R||+||D||*||L||) + SCALE * ||F|| ))    (1)

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

 \param[in] SGN1
 \verbatim
          SGN1 is DOUBLE PRECISION, allowed values: +/-1
          Specifies the sign in the first equation.
 \endverbatim

 \param[in] SGN2
 \verbatim
          SGN2 is DOUBLE PRECISION, allowed values: +/-1
          Specifies the sign in the second equation.
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

 \param[in] R
 \verbatim
          R is DOUBLE PRECISION array, dimension (LDR,N)
          The first of the Sylvester Equation.
 \endverbatim

 \param[in] LDR
 \verbatim
          LDR is INTEGER
          The leading dimension of the array R.  LDX >= max(1,M).
 \endverbatim

 \param[in] L
 \verbatim
          L is DOUBLE PRECISION array, dimension (LDL,N)
          The second solution the Sylvester Equation.
 \endverbatim

 \param[in] LDL
 \verbatim
          LDL is INTEGER
          The leading dimension of the array L.  LDL >= max(1,M).
 \endverbatim


 \param[in] E
 \verbatim
          E is DOUBLE PRECISION array, dimension (LDE,N)
          The first right hand side of the Sylvester Equation.
 \endverbatim

 \param[in] LDE
 \verbatim
          LDE is INTEGER
          The leading dimension of the array E.  LDE >= max(1,M).
 \endverbatim

 \param[in] F
 \verbatim
          F is DOUBLE PRECISION array, dimension (LDF,N)
          The second right hand side of the Sylvester Equation.
 \endverbatim

 \param[in] LDF
 \verbatim
          LDF is INTEGER
          The leading dimension of the array F.  LDF >= max(1,M).
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

 \see dla_residual_csylv

 \author Martin Koehler, MPI Magdeburg

 \date January 2024


 */

double mepack_double_residual_csylv(const char * TRANSA, const char *TRANSB, double sgn1, double sgn2, int M, int N, double *A, int LDA,
        double *B, int LDB, double *C, int LDC, double *D, int LDD,
        double *R, int LDR, double *L, int LDL,
        double *E, int LDE, double *F, int LDF, double SCALE)
{
    Int _M   = M;
    Int _N   = N;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDC = LDC;
    Int _LDD = LDD;
    Int _LDE = LDE;
    Int _LDF = LDF;
    Int _LDR = LDR;
    Int _LDL = LDL;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    return FC_GLOBAL_(dla_residual_csylv, DLA_RESIDUAL_CSYLV) (TRANSA, TRANSB, &sgn1, &sgn2, &_M, &_N,
            A, &_LDA, B, &_LDB, C, &_LDC, D, &_LDD,
            R, &_LDR, L, &_LDL, E, &_LDE, F, &_LDF, &SCALE, c1len, c2len);
}
/**
 \brief Compute the relative residual for the coupled Sylvester equation.

 \par Purpose:

 \verbatim

 Compute the relative residual of the coupled Sylvester equation

 RelRes = max( || SCALE*E - opA(A)*R - SGN1 * L *opB(B) || / ( (||A||*||R||+||B||*||L||) + SCALE * ||E|| ),
               || SCALE*F - opA(C)*R - SGN2 * L *opB(D) || / ( (||C||*||R||+||D||*||L||) + SCALE * ||F|| ))    (1)

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

 \param[in] SGN1
 \verbatim
          SGN1 is SINGLE PRECISION, allowed values: +/-1
          Specifies the sign in the first equation.
 \endverbatim

 \param[in] SGN2
 \verbatim
          SGN2 is SINGLE PRECISION, allowed values: +/-1
          Specifies the sign in the second equation.
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

 \param[in] R
 \verbatim
          R is SINGLE PRECISION array, dimension (LDR,N)
          The first of the Sylvester Equation.
 \endverbatim

 \param[in] LDR
 \verbatim
          LDR is INTEGER
          The leading dimension of the array R.  LDX >= max(1,M).
 \endverbatim

 \param[in] L
 \verbatim
          L is SINGLE PRECISION array, dimension (LDL,N)
          The second solution the Sylvester Equation.
 \endverbatim

 \param[in] LDL
 \verbatim
          LDL is INTEGER
          The leading dimension of the array L.  LDL >= max(1,M).
 \endverbatim


 \param[in] E
 \verbatim
          E is SINGLE PRECISION array, dimension (LDE,N)
          The first right hand side of the Sylvester Equation.
 \endverbatim

 \param[in] LDE
 \verbatim
          LDE is INTEGER
          The leading dimension of the array E.  LDE >= max(1,M).
 \endverbatim

 \param[in] F
 \verbatim
          F is SINGLE PRECISION array, dimension (LDF,N)
          The second right hand side of the Sylvester Equation.
 \endverbatim

 \param[in] LDF
 \verbatim
          LDF is INTEGER
          The leading dimension of the array F.  LDF >= max(1,M).
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

 \see sla_residual_csylv

 \author Martin Koehler, MPI Magdeburg

 \date January 2024


 */

float mepack_single_residual_csylv(const char * TRANSA, const char *TRANSB, float sgn1, float sgn2, int M, int N, float *A, int LDA,
        float *B, int LDB, float *C, int LDC, float *D, int LDD,
        float *R, int LDR, float *L, int LDL,
        float *E, int LDE, float *F, int LDF, float SCALE)
{
    Int _M   = M;
    Int _N   = N;
    Int _LDA = LDA;
    Int _LDB = LDB;
    Int _LDC = LDC;
    Int _LDD = LDD;
    Int _LDE = LDE;
    Int _LDF = LDF;
    Int _LDR = LDR;
    Int _LDL = LDL;
    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    return FC_GLOBAL_(sla_residual_csylv, SLA_RESIDUAL_CSYLV) (TRANSA, TRANSB, &sgn1, &sgn2, &_M, &_N,
            A, &_LDA, B, &_LDB, C, &_LDC, D, &_LDD,
            R, &_LDR, L, &_LDL, E, &_LDE, F, &_LDF, &SCALE, c1len, c2len);
}



/**
 * @}
 */

