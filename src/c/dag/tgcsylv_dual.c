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
 \brief OpenMP4 - DAG  Algorithm for the dual generalized coupled Sylvester equation

 \par Purpose:

 \verbatim

 mepack_double_tgcsylv_dual_dag solves a generalized coupled  Sylvester equation of the following form

    op1(A)**T * R  + op1(C)**T * L               = SCALE * E                             (1)
    SGN1 * R * op2(B)**T + SGN2 * L * op2(D)** T = SCALE * F

 where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
 matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
 D can be quasi upper triangular as well. The right hand side Y and the solution X
 M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
 created by DGGES from LAPACK.
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

 \remark This function is a wrapper for \ref dla_tgcsylv_dual_dag.

 \see dla_tgcsylv_dual_dag

 \param[in] TRANSA
 \verbatim
          TRANSA is CHARACTER(1)
          Specifies the form of the system of equations with respect to A and C:
          == 'N':  op1(A) = A, op1(C) = C (No transpose for A and C)
          == 'T':  op1(A) = A**T, op1(C) = C **T (Transpose A and C)
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is CHARACTER(1)
          Specifies the form of the system of equations with respect to B and D:
          == 'N':  op2(B) = B, op2(D) = D (No transpose for B and D)
          == 'T':  op2(B) = B**T, op2(D) = D **T (Transpose B and D)
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

 \param[in,out] E
 \verbatim
          E is DOUBLE PRECISION array, dimension (LDE,N)
          On input, the matrix E contains the right hand side E.
          On output, the matrix E contains the solution R of Equation (1)
          as selected by TRANSA, TRANSB, and SGN1/SGN2.
          Right hand side E and the solution R are M-by-N matrices.
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
          On output, the matrix F contains the solution L of Equation (1)
          as selected by TRANSA, TRANSB, and SGN1/SGN2.
          Right hand side F and the solution L are M-by-N matrices.
 \endverbatim

 \param[in] LDF
 \verbatim
          LDF is INTEGER
          The leading dimension of the array F.  LDE >= max(1,M).
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
          Workspace for the algorithm. The optimal workspace is given by \ref mepack_memory.
 \endverbatim

 \param[inout] INFO
 \verbatim
          INFO is INTEGER

          On input:
            == -1 : Perform a workspace query ( not recommended in the C interface)
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


void mepack_double_tgcsylv_dual_dag(const char *TRANSA, const char*TRANSB, double SGN1, double SGN2,  int M, int N,
        double * A, int LDA,double * B, int LDB, double *C, int LDC, double *D, int LDD,
        double *E, int LDE, double *F, int LDF ,  double * SCALE, double *WORK, int *INFO)
{
    Int _LDA  = LDA;
    Int _LDB  = LDB;
    Int _LDC  = LDC;
    Int _LDD  = LDD;
    Int _LDE  = LDE;
    Int _LDF  = LDF;
    Int _INFO = *INFO;
    Int _M    = M;
    Int _N    = N;

    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(dla_tgcsylv_dual_dag,DLA_TGCSYLV_DUAL_DAG)(TRANSA, TRANSB, &SGN1, &SGN2, &_M, &_N, A, &_LDA, B, &_LDB,
            C, &_LDC, D, &_LDD, E, &_LDE, F, &_LDF, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO ;
    return;


}

/**
 \brief OpenMP4 - DAG  Algorithm for the dual generalized coupled Sylvester equation

 \par Purpose:

 \verbatim

 mepack_single_tgcsylv_dual_dag solves a generalized coupled  Sylvester equation of the following form

    op1(A)**T * R  + op1(C)**T * R               = SCALE * E                             (1)
    SGN1 * R * op2(B)**T + SGN2 * L * op2(D)** T = SCALE * F

 where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
 matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
 D can be quasi upper triangular as well. The right hand side Y and the solution X
 M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
 created by DGGES from LAPACK.
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

 \remark This function is a wrapper for \ref sla_tgcsylv_dual_dag.

 \see sla_tgcsylv_dual_dag

 \param[in] TRANSA
 \verbatim
          TRANSA is CHARACTER(1)
          Specifies the form of the system of equations with respect to A and C:
          == 'N':  op1(A) = A, op1(C) = C (No transpose for A and C)
          == 'T':  op1(A) = A**T, op1(C) = C **T (Transpose A and C)
 \endverbatim

 \param[in] TRANSB
 \verbatim
          TRANSB is CHARACTER(1)
          Specifies the form of the system of equations with respect to B and D:
          == 'N':  op2(B) = B, op2(D) = D (No transpose for B and D)
          == 'T':  op2(B) = B**T, op2(D) = D **T (Transpose B and D)
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

 \param[in,out] E
 \verbatim
          E is SINGLE PRECISION array, dimension (LDE,N)
          On input, the matrix E contains the right hand side E.
          On output, the matrix E contains the solution R of Equation (1)
          as selected by TRANSA, TRANSB, and SGN1/SGN2.
          Right hand side E and the solution R are M-by-N matrices.
 \endverbatim

 \param[in] LDE
 \verbatim
          LDE is INTEGER
          The leading dimension of the array E.  LDE >= max(1,M).
 \endverbatim

 \param[in,out] F
 \verbatim
          F is SINGLE PRECISION array, dimension (LDF,N)
          On input, the matrix F contains the right hand side F.
          On output, the matrix F contains the solution L of Equation (1)
          as selected by TRANSA, TRANSB, and SGN1/SGN2.
          Right hand side F and the solution L are M-by-N matrices.
 \endverbatim

 \param[in] LDF
 \verbatim
          LDF is INTEGER
          The leading dimension of the array F.  LDE >= max(1,M).
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
          Workspace for the algorithm. The optimal workspace is given by \ref mepack_memory.
 \endverbatim

 \param[inout] INFO
 \verbatim
          INFO is INTEGER

          On input:
            == -1 : Perform a workspace query ( not recommended in the C interface)
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
void mepack_single_tgcsylv_dual_dag(const char *TRANSA, const char*TRANSB, float SGN1, float SGN2,  int M, int N,
        float * A, int LDA,float * B, int LDB, float *C, int LDC, float *D, int LDD,
        float *E, int LDE, float *F, int LDF ,  float * SCALE, float *WORK, int *INFO)
{
    Int _LDA  = LDA;
    Int _LDB  = LDB;
    Int _LDC  = LDC;
    Int _LDD  = LDD;
    Int _LDE  = LDE;
    Int _LDF  = LDF;
    Int _INFO = *INFO;
    Int _M    = M;
    Int _N    = N;

    FORTRAN_LEN_T c1len = strlen(TRANSA);
    FORTRAN_LEN_T c2len = strlen(TRANSB);

    FC_GLOBAL_(sla_tgcsylv_dual_dag,SLA_TGCSYLV_DUAL_DAG)(TRANSA, TRANSB, &SGN1, &SGN2, &_M, &_N, A, &_LDA, B, &_LDB,
            C, &_LDC, D, &_LDD, E, &_LDE, F, &_LDF, SCALE, WORK, &_INFO, c1len, c2len);
    *INFO = (int) _INFO ;
    return;


}

/**
 * @}
 */

