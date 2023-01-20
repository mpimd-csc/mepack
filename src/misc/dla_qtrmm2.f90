! 
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, see <http://www.gnu.org/licenses/>.
! 
! Copyright (C) Martin Koehler, 2017-2023
! 

!> \brief Multiply with the sub-diagonal entries of a quasi triangular matrix. 
!
!  Definition:
!  ===========
!          SUBROUTINE DLA_QTRMM2(SIDE, TRANS, M, N, ALPHA, A, LDA, B, LDB, C, LDC)
!              IMPLICIT NONE
!              INTEGER LDA,LDB,LDC,M,N
!              DOUBLE PRECISION A(LDA, *), B(LDB, *), C(LDC,*)
!              DOUBLE PRECISION ALPHA
!              CHARACTER(1) SIDE, TRANS
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>     Multiply a matrix with the sub-diagonal entries of a quasi-triangular matrix. 
!>     Consider a quasi-triangular matrix A where sub(A) is the matrix only containing 
!>     the sub-diagonal entries. Than the DLA_QTRMM2 routine computes 
!> 
!>        C = C + ALPHA * op(sub(A)) * B                                                (1) 
!> 
!>     or 
!> 
!>        C = C + ALPHA * B * op(sub(A))                                                (2) 
!> 
!> \endverbatim 
!>
!> \attention The routine does not check the input arguments. 
!
!  Arguments:
!  ==========
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER(1)
!>           On entry,  SIDE specifies whether  op( A ) multiplies B from
!>           the left or right as follows:
!>
!>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!>
!>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER(1)
!>           On entry, TRANS specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANS == 'N' or 'n'   op( A ) = A.
!>
!>              TRANS == 'T' or 't'   op( A ) = A**T.
!>
!>              TRANS = 'C' or 'c'   op( A ) = A**T.
!> \endverbatim
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of B. M must be at
!>           least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of B.  N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>           A is DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!>           The  leading  k by k 
!>           upper triangular part of the array  A must contain the upper
!>           quasi triangular matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!>           then LDA must be at least max( 1, n ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>           B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain the matrix  B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>           C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!>           Before entry,  the leading  m by n part of the array  C must
!>           contain the matrix  C,  and  on exit  is overwritten  by the
!>           update January 2021
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>           On entry, LDC specifies the first dimension of C as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg 
!
!> \date Januar 2023
!> \ingroup dblaux 
!
SUBROUTINE DLA_QTRMM2(SIDE, TRANS, M, N, ALPHA, A, LDA, B, LDB, C, LDC)
    ! Arguments
    IMPLICIT NONE
    INTEGER LDA,LDB,LDC,M,N
    DOUBLE PRECISION A(LDA, *), B(LDB, *), C(LDC,*)
    DOUBLE PRECISION ALPHA
    CHARACTER(1) SIDE, TRANS

    ! Local Variables 
    INTEGER I 
    DOUBLE PRECISION ZERO 
    PARAMETER (ZERO = 0.0D0) 

    ! External Routines 
    EXTERNAL DAXPY 

    IF (SIDE.EQ.'L') THEN
        IF (TRANS.EQ.'N') THEN
            DO I = 1,M-1
                IF (A(I+1,I).NE.ZERO) THEN
                    CALL DAXPY(N, ALPHA*A(I+1,I), B(I,1), LDB, C(I+1,1), LDC)
                END IF
            END DO
        ELSE
            DO I = 1,M-1
                IF (A(I+1,I).NE.ZERO) THEN
                    CALL DAXPY(N, ALPHA*A(I+1,I), B(I+1,1), LDB, C(I,1), LDC)
                END IF
            END DO
        END IF
    ELSE
        IF (TRANS.EQ.'N') THEN
            DO I = 1,N-1
                IF (A(I+1,I).NE. ZERO) THEN
                    CALL DAXPY(M, ALPHA*A(I+1,I), B(1,I+1), 1, C(1,I), 1)
                END IF
            END DO
        ELSE
            DO I = 1,N-1
                IF (A(I+1,I).NE.ZERO) THEN
                    CALL DAXPY(M, ALPHA*A(I+1,I), B(1,I), 1, C(1,I+1), 1)
                END IF
            END DO
        END IF
    END IF

END SUBROUTINE 

