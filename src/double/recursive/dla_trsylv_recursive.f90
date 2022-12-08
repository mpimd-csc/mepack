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
! Copyright (C) Martin Koehler, 2017-2022
!


!> \brief Level-3 Recursive Blocked Solver for the Sylvester equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TRSYLV_RECURSIVE ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANSA(1), TRANSB(1)
!           DOUBLE PRECISION SGN, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
!           DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
!           DOUBLE PRECISION WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> DLA_TRSYLV_RECURSIVE solves a Sylvester equation of the following forms
!>
!>    op1(A) * X  +  X * op2(B) = SCALE * Y                              (1)
!>
!> or
!>
!>    op1(A) * X  -  X * op2(B) = SCALE * Y                              (2)
!>
!> where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
!> matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
!> A and B are generated via DGEES form LAPACK.
!> \endverbatim
!> \attention The algorithm is implemented using BLAS level 3 operations and recursive blocking.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A:
!>          == 'N':  op1(A) = A
!>          == 'T':  op1(A) = A**T
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to B:
!>          == 'N':  op2(B) = B,
!>          == 'T':  op2(B) = B**T
!> \endverbatim
!>
!> \param[in] SGN
!> \verbatim
!>          SGN is DOUBLE PRECISION, allowed values: +/-1
!>          Specifies the sign between the two parts of the Sylvester equation.
!>          = 1 :  Solve Equation (1)
!>          == -1:  Solve Equation (2)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and C.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices B and D.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,M)
!>          The matrix A must be (quasi-) upper triangular.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          The matrix B must be (quasi-) upper triangular.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,N)
!>          On input, the matrix X contains the right hand side Y.
!>          On output, the matrix X contains the solution of Equation (1) or (2)
!>          as selected by TRANSA, TRANSB, and SGN.
!>          Right hand side Y and the solution X are M-by-N matrices.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          SCALE is a scaling factor to prevent the overflow in the result.
!>          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
!>          could not be solved correctly, 0 < SCALE <= 1 holds true.
!> \endverbatim
!>
!> \param[in] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension 1
!>          Workspace for the algorithm
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          == 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  The equation is not solved correctly. One of the arising inner
!>                system got singular.
!> \endverbatim
!>
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Dezember 2022
!> \ingroup dbltrsylv
!
RECURSIVE SUBROUTINE DLA_TRSYLV_RECURSIVE(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    DOUBLE PRECISION SGN, SCALE
    INTEGER M, N, LDA, LDB, LDX, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *) , X(LDX, *)
    DOUBLE PRECISION WORK(*)

    ! Local Variables
    INTEGER ININFO
    DOUBLE PRECISION SCAL

    LOGICAL BTRANSA, BTRANSB
    DOUBLE PRECISION ZERO , ONE
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
    INTEGER IZERO, IONE
    PARAMETER(IZERO = 0, IONE = 1 )

    INTEGER M1, M2, N1, N2, IINFO

    INTEGER BOUND
    PARAMETER(BOUND = 4)

    ! EXTERNAL FUNCTIONS
    EXTERNAL DGEMM
    EXTERNAL DLA_TRSYLV_KERNEL_44NN
    EXTERNAL DLA_TRSYLV_KERNEL_44NT
    EXTERNAL DLA_TRSYLV_KERNEL_44TN
    EXTERNAL DLA_TRSYLV_KERNEL_44TT
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME


    ! Check Input
    ININFO = INFO
    INFO = 0
    BTRANSA = LSAME(TRANSA, 'T')
    BTRANSB = LSAME(TRANSB, 'T')

    IF ( .NOT. BTRANSA .AND. .NOT. LSAME(TRANSA, 'N')) THEN
        INFO = -1
    ELSE IF ( .NOT. BTRANSB .AND. .NOT. LSAME(TRANSB, 'N')) THEN
        INFO = -2
    ELSE IF ( SGN .NE. -ONE .AND. SGN .NE. ONE ) THEN
        INFO = -3
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -4
    ELSE IF ( N .LT. 0 ) THEN
        INFO = -5
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -7
    ELSE IF ( LDB .LT. MAX(1, N)) THEN
        INFO = -9
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -11
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('DLA_TRSYLV_RECURSIVE', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1 ) THEN
        INFO = 1
        RETURN
    END IF


    ! Quick Return
    IF (M .LE. 0 .OR. N.LE.0 ) THEN
        RETURN
    ENDIF

    SCALE = ONE

    IF (M .LE. BOUND .AND. N .LE. BOUND) THEN
        IINFO = 0
        IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB) THEN
            !(N,N)
            CALL DLA_TRSYLV_KERNEL_44NN(SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, IINFO)
        ELSE IF ( .NOT. BTRANSA .AND. BTRANSB) THEN
            !(N,T)
            CALL DLA_TRSYLV_KERNEL_44NT(SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE,IINFO)
        ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
            !(T,N)
            CALL DLA_TRSYLV_KERNEL_44TN(SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE,IINFO)
        ELSE IF ( BTRANSA .AND. BTRANSB) THEN
            !(T,T)
            CALL DLA_TRSYLV_KERNEL_44TT(SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE,IINFO)
        ENDIF
        IF ( IINFO .NE. 0 ) THEN
            INFO = 100
        END IF
        RETURN
    ENDIF


    IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB ) THEN
        ! (N,N)
        IF ( N .LE. M/2 .OR. N .LE. BOUND ) THEN
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF

            CALL DLA_TRSYLV_RECURSIVE('N', 'N', SGN, M2, N, A(M1+1,M1+1), LDA, B, LDB, X(M1+1,1), LDX, SCALE, WORK, INFO )
            CALL DGEMM('N', 'N', M1, N, M2, -ONE, A(1,M1+1), LDA, X(M1+1, 1), LDX, SCALE, X(1,1), LDX)
            CALL DLA_TRSYLV_RECURSIVE('N', 'N', SGN, M1, N, A(1,1), LDA, B, LDB, X(1,1), LDX, SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(M1+1:M, 1:N) = SCAL*X(M1+1:M, 1:N)
                SCALE = SCAL * SCALE
            ENDIF
        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO ) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TRSYLV_RECURSIVE('N', 'N', SGN, M, N1, A, LDA, B(1,1), LDB, X(1,1), LDX, SCALE, WORK, INFO)
            CALL DGEMM('N', 'N', M, N2, N1, -SGN, X(1,1), LDX, B(1, N1+1), LDB, SCALE, X(1,N1+1), LDX)
            CALL DLA_TRSYLV_RECURSIVE('N', 'N', SGN, M, N2, A, LDA, B(N1+1,N1+1), LDB, X(1,N1+1), LDX, SCAL,  WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M, 1:N1) = SCAL*X(1:M, 1:N1)
                SCALE = SCAL * SCALE
            ENDIF
        ELSE
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO ) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TRSYLV_RECURSIVE('N','N', SGN, M2, N1, A(M1+1,M1+1), LDA, B(1,1), LDB, X(M1+1,1), LDX,SCALE,WORK,INFO)

            CALL DGEMM('N', 'N', M1, N1, M2, -ONE, A(1,M1+1), LDA, X(M1+1,1), LDX, ONE, X(1,1), LDX)

            CALL DLA_TRSYLV_RECURSIVE('N','N', SGN, M1, N1, A(1,1), LDA, B(1,1), LDB, X(1,1), LDX, SCAL, WORK, INFO)

            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(M1+1:M, 1:N1) = SCAL * X(M1+1:M, 1:N1)
            ENDIF

            CALL DGEMM('N', 'N', M2, N2, N1, -SGN, X(M1+1,1), LDX, B(1,N1+1), LDB, SCALE, X(M1+1,N1+1), LDX)

            CALL DLA_TRSYLV_RECURSIVE('N','N', SGN, M2, N2, A(M1+1, M1+1), LDA, B(N1+1,N1+1), LDB,&
                & X(M1+1, N1+1), LDX, SCAL, WORK, INFO)
            IF (SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M, 1:N1) = SCAL * X(1:M, 1:N1)
            ENDIF

            CALL DGEMM('N', 'N', M1, N2, M2, -ONE, A(1,M1+1), LDA, X(M1+1,N1+1), LDX, SCALE, X(1, N1+1), LDX)
            CALL DGEMM('N', 'N', M1, N2, N1, -SGN, X(1,1), LDX, B(1,N1+1), LDB, ONE, X(1,N1+1), LDX)

            CALL DLA_TRSYLV_RECURSIVE('N', 'N',SGN,  M1, N2, A(1,1), LDA, B(N1+1, N1+1), LDB, &
                & X(1,N1+1), LDX, SCAL, WORK, INFO)
            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M, 1:N1) = SCAL * X(1:M, 1:N1)
                X(M1+1:M, N1+1:N) = SCAL * X(M1+1:M, N1+1:N)
            ENDIF
        END IF


    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N, T)
        IF ( N .LE. M/2 .OR. N .LE. BOUND ) THEN
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF

            CALL DLA_TRSYLV_RECURSIVE('N', 'T', SGN, M2, N, A(M1+1,M1+1), LDA, B, LDB, X(M1+1,1), LDX, SCALE, WORK, INFO )
            CALL DGEMM('N', 'N', M1, N, M2, -ONE, A(1,M1+1), LDA, X(M1+1, 1), LDX, SCALE, X(1,1), LDX)

            CALL DLA_TRSYLV_RECURSIVE('N', 'T', SGN, M1, N, A(1,1), LDA, B, LDB, X(1,1), LDX, SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(M1+1:M, 1:N) = SCAL*X(M1+1:M, 1:N)
                SCALE = SCAL * SCALE
            ENDIF
        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO ) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TRSYLV_RECURSIVE('N', 'T', SGN, M, N2, A, LDA, B(N1+1, N1+1), LDB, X(1,N1+1), LDX, SCALE, WORK, INFO)
            CALL DGEMM('N', 'T', M, N1, N2, -SGN, X(1, N1+1), LDX, B(1, N1+1), LDB, SCALE, X(1,1), LDX)
            CALL DLA_TRSYLV_RECURSIVE('N', 'T',  SGN, M, N1, A, LDA, B(1,1), LDB, X(1,1), LDX, SCAL,  WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M, N1+1:N) = SCAL*X(1:M, N1+1:N)
                SCALE = SCAL * SCALE
            ENDIF
        ELSE
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO ) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TRSYLV_RECURSIVE('N','T', SGN, M2, N2, A(M1+1,M1+1), LDA, B(N1+1,N1+1), LDB, &
                & X(M1+1,N1+1), LDX, SCALE, WORK, INFO)

            CALL DGEMM('N', 'T', M2, N1, N2, -SGN, X(M1+1, N1+1), LDX, B(1, N1+1), LDB, SCALE, X(M1+1,1), LDX)
            CALL DLA_TRSYLV_RECURSIVE('N','T', SGN,  M2, N1, A(M1+1,M1+1), LDA, B(1,1), LDB, &
                & X(M1+1,1), LDX, SCAL, WORK, INFO)

            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(M1+1:M, N1+1:N) = SCAL * X(M1+1:M, N1+1:N)
            ENDIF

            CALL DGEMM('N', 'N', M1, N2, M2, -ONE, A(1,M1+1), LDA, X(M1+1, N1+1), LDX, SCALE, X(1,N1+1), LDX)
            CALL DLA_TRSYLV_RECURSIVE('N','T', SGN,  M1, N2, A(1,1), LDA, B(N1+1,N1+1), LDB, &
                & X(1, N1+1), LDX, SCAL, WORK, INFO)
            IF (SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(M1+1, 1:N) = SCAL * X(M1+1, 1:N)
            ENDIF

            CALL DGEMM('N', 'N', M1, N1, M2, -ONE, A(1,M1+1), LDA, X(M1+1,1), LDX, SCALE, X(1,1), LDX)
            CALL DGEMM('N', 'T', M1, N1, N2, -SGN, X(1,N1+1), LDX, B(1,N1+1), LDB, ONE, X(1,1), LDX)


            CALL DLA_TRSYLV_RECURSIVE('N', 'T',SGN,  M1, N1, A(1,1), LDA, B(1,1), LDB, X(1,1), LDX, SCAL, WORK, INFO)
            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M, 1:N1) = SCAL * X(1:M, 1:N1)
                X(M1+1:M, N1+1:N) = SCAL * X(M1+1:M, N1+1:N)
            ENDIF
        END IF

    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
        ! (T, N)
        IF ( N .LE. M/2 .OR. N .LE. BOUND ) THEN
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF

            CALL DLA_TRSYLV_RECURSIVE('T', 'N', SGN,  M1, N, A(1,1), LDA, B, LDB, X(1,1), LDX,  SCALE, WORK, INFO )
            CALL DGEMM('T', 'N', M2, N, M1, -ONE, A(1,M1+1), LDA, X(1, 1), LDX, SCALE, X(M1+1,1), LDX)
            CALL DLA_TRSYLV_RECURSIVE('T', 'N', SGN,  M2, N, A(M1+1,M1+1), LDA, B, LDB, X(M1+1,1), LDX, SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M1, 1:N) = SCAL*X(1:M1, 1:N)
                SCALE = SCAL * SCALE
            ENDIF

        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF ( B(N1+1,N1) .NE. ZERO ) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TRSYLV_RECURSIVE('T', 'N', SGN,  M, N1, A, LDA, B(1,1), LDB, X(1,1), LDX, SCALE, WORK, INFO)
            CALL DGEMM('N', 'N', M, N2, N1, -SGN, X(1,1), LDX, B(1, N1+1), LDB, SCALE, X(1,N1+1), LDX)

            CALL DLA_TRSYLV_RECURSIVE('T', 'N', SGN,  M, N2, A, LDA, B(N1+1,N1+1), LDB, X(1,N1+1), LDX, SCAL,  WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M, 1:N1) = SCAL*X(1:M, 1:N1)
                SCALE = SCAL * SCALE
            ENDIF
        ELSE
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF
            N1 = N/2
            N2 = N - N1
            IF ( B(N1+1,N1) .NE. ZERO ) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF


            CALL DLA_TRSYLV_RECURSIVE('T','N', SGN,  M1, N1, A(1,1), LDA, B(1,1), LDB , X(1,1), LDX, &
                & SCALE, WORK, INFO)

            CALL DGEMM("T", "N", M2, N1, M1, -ONE, A(1,M1+1), LDA, X(1,1), LDX, SCALE, X(M1+1,1), LDX)
            CALL DLA_TRSYLV_RECURSIVE('T','N', SGN,  M2, N1, A(M1+1, M1+1), LDA, B(1,1), LDB, &
                & X(M1+1, 1), LDX, SCAL, WORK, INFO)
            IF (SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, 1:N) = SCAL * X(1:M1, 1:N)
            ENDIF

            CALL DGEMM("N", "N", M1, N2, N1, -SGN, X(1,1), LDX, B(1,N1+1), LDB, SCALE, X(1,N1+1), LDX)
            CALL DLA_TRSYLV_RECURSIVE('T','N', SGN,  M1, N2, A(1,1), LDA, B(N1+1,N1+1), LDB, &
                & X(1,N1+1), LDX, SCAL, WORK, INFO)

            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, 1:N1) = SCAL * X(1:M1, 1:N1)
            ENDIF

            !X22
            CALL DGEMM('T', 'N', M2, N2, M1, -ONE, A(1,M1+1), LDA, X(1,N1+1), LDX, SCALE, X(M1+1, N1+1), LDX)
            CALL DGEMM('N', 'N', M2, N2, N1, -SGN, X(M1+1,1), LDX, B(1,N1+1), LDB, ONE, X(M1+1, N1+1), LDX)

            CALL DLA_TRSYLV_RECURSIVE('T', 'N', SGN,  M2, N2, A(M1+1,M1+1), LDA, B(N1+1, N1+1), LDB,  &
                & X(M1+1,N1+1), LDX, SCAL, WORK, INFO)
            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, 1:N) = SCAL * X(1:M1, 1:N)
                X(M1+1:M, 1:N1) = SCAL * X(M1+1:M, 1:N1)
            ENDIF
        END IF
    ELSE
        ! (T,T)
        IF ( N .LE. M/2 .OR. N .LE. BOUND ) THEN
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF

            !X11
            CALL DLA_TRSYLV_RECURSIVE('T', 'T', SGN,  M1, N, A(1,1), LDA, B, LDB, X(1,1), LDX, SCALE, WORK, INFO )
            CALL DGEMM('T', 'N', M2, N, M1, -ONE, A(1,M1+1), LDA, X(1, 1), LDX, SCALE, X(M1+1,1), LDX)

            ! X22
            CALL DLA_TRSYLV_RECURSIVE('T', 'T', SGN,  M2, N, A(M1+1,M1+1), LDA, B, LDB, X(M1+1,1), LDX, SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M1, 1:N) = SCAL*X(1:M1, 1:N)
                SCALE = SCAL * SCALE
            ENDIF

        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO ) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TRSYLV_RECURSIVE('T', 'T', SGN,  M, N2, A, LDA, B(N1+1, N1+1), LDB, X(1,N1+1), LDX, SCALE, WORK, INFO)
            CALL DGEMM('N', 'T', M, N1, N2, -SGN, X(1, N1+1), LDX, B(1, N1+1), LDB, SCALE, X(1,1), LDX)

            CALL DLA_TRSYLV_RECURSIVE('T', 'T', SGN,  M, N1, A, LDA, B(1,1), LDB, X(1,1), LDX, SCAL,  WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M, N1+1:N) = SCAL*X(1:M, N1+1:N)
                SCALE = SCAL * SCALE
            ENDIF

        ELSE
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO ) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            !X12
            CALL DLA_TRSYLV_RECURSIVE('T','T', SGN,  M1, N2, A(1,1), LDA, B(N1+1,N1+1), LDB , &
                & X(1,N1+1), LDX, SCALE, WORK, INFO)

            !X11
            CALL DGEMM('N', 'T', M1, N1, N2, -SGN, X(1, N1+1), LDX, B(1, N1+1), LDB, SCALE, X(1,1), LDX)
            CALL DLA_TRSYLV_RECURSIVE('T','T', SGN,  M1, N1, A(1,1), LDA, B(1,1), LDB,   &
                & X(1,1), LDX, SCAL, WORK, INFO)

            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, N1+1) = SCAL * X(1:M1, N1+1)
            ENDIF

            !X22
            CALL DGEMM("T", "N", M2, N2, M1, -ONE, A(1, M1+1), LDA, X(1, N1+1), LDX, SCALE, X(M1+1, N1+1), LDX)

            CALL DLA_TRSYLV_RECURSIVE('T','T', SGN,  M2, N2, A(M1+1, M1+1), LDA, B(N1+1,N1+1), LDB, &
                & X(M1+1, N1+1), LDX, SCAL, WORK, INFO)
            IF (SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, 1:N) = SCAL * X(1:M1, 1:N)
            ENDIF

            ! X21
            CALL DGEMM('T', 'N', M2, N1, M1, -ONE, A(1,M1+1), LDA, X(1,1), LDX, SCALE, X(M1+1, 1), LDX)
            CALL DGEMM('N', 'T', M2, N1, N2, -SGN, X(M1+1,N1+1), LDX, B(1,N1+1), LDB, ONE, X(M1+1,1), LDX)



            CALL DLA_TRSYLV_RECURSIVE('T', 'T',SGN,  M2, N1, A(M1+1,M1+1), LDA, B(1,1), LDB,&
                & X(M1+1,1), LDX, SCAL, WORK, INFO)
            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, 1:N) = SCAL * X(1:M1, 1:N)
                X(M1+1:M, N1+1:N) = SCAL * X(M1+1:M, N1+1:N)
            ENDIF

        ENDIF

    ENDIF
END SUBROUTINE

