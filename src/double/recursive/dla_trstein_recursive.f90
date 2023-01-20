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


!> \brief Recursive Blocked Algorithm for the Stein equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TRSTEIN_RECURSIVE ( TRANS, M,  A, LDA, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANS(1)
!           DOUBLE PRECISION SCALE
!           INTEGER M, N, LDA, LDX, INFO
!           DOUBLE PRECISION A(LDA, *), X(LDX, *)
!           DOUBLE PRECISION WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLA_TRSTEIN_RECURSIVE solves a generalized Lyapunov  equation of the following forms
!>
!>    A * X * A^T - X  = SCALE * Y                                              (2)
!>
!> or
!>
!>    A^T * X * A - X  =  SCALE * Y                                             (1)
!>
!> where A is a M-by-M quasi upper triangular matrix.
!> and X and Y are symmetric  M-by-M matrices.
!> Typically the matrix A is created by DGEES from LAPACK.
!> \endverbatim
!> \attention The algorithm is implemented using recursive blocking.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A and C:
!>          == 'N':  Equation (1) is solved.
!>          == 'T':  Equation (2) is solved.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and C.  M >= 0.
!> \endverbatim
!>
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
!>          The leading dimension of the array X.  LDX >= max(1,M).
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
!>          WORK is DOUBLE PRECISION array, dimension M*N
!>          Workspace for the algorithm.
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
!> \date Januar 2023
!> \ingroup dbltrlyap
!
RECURSIVE SUBROUTINE DLA_TRSTEIN_RECURSIVE ( TRANS, M, A, LDA, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANS(1)
    DOUBLE PRECISION SCALE
    INTEGER M, LDA, LDX, INFO
    DOUBLE PRECISION A(LDA, *), X(LDX, *)
    DOUBLE PRECISION WORK(*)


    ! Local Variables
    INTEGER ININFO
    DOUBLE PRECISION SCAL
    DOUBLE PRECISION MAT(4,4)
    DOUBLE PRECISION RHS(4)
    INTEGER IPIV(4), JPIV(4)

    LOGICAL BTRANS
    CHARACTER(1) TRANSA, TRANSB
    DOUBLE PRECISION ZERO , ONE, SGN, HALF, TWO
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, SGN = -1.0D0, HALF = 0.5D0, TWO = 2.0D0 )
    INTEGER IZERO, IONE, INFO1, K, KH, KB
    PARAMETER(IZERO = 0, IONE = 1 )

    INTEGER M1, M2, IINFO, I


    INTEGER BOUND
    PARAMETER(BOUND = 2)

    ! EXTERNAL FUNCTIONS
    EXTERNAL DGEMM
    EXTERNAL DGESC2
    EXTERNAL DGETC2X
    EXTERNAL DLA_DTRMM_FIX
    EXTERNAL DLA_ITRANSPOSE
    EXTERNAL DLA_SMALL_SOLVE4
    EXTERNAL DLA_TRSYLV2_RECURSIVE
    EXTERNAL DSYR2K
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME

    ! Check Input
    ININFO = INFO
    INFO = 0
    BTRANS = LSAME(TRANS, 'N')

    IF ( .NOT. BTRANS .AND. .NOT. LSAME(TRANS, 'T')) THEN
        INFO = -1
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -2
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -4
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -6
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('DLA_TRSTEIN_RECURSIVE', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1 ) THEN
        INFO = (M/2+2)**2
        RETURN
    END IF


    ! Quick Return
    SCALE = ONE
    SCAL  = ONE
    INFO = 0
    INFO1 = 0
    IF ( M .EQ. 0 ) THEN
        RETURN
    END IF

    IF ( BTRANS ) THEN
        TRANSA = 'N'
        TRANSB = 'T'
    ELSE
        TRANSA = 'T'
        TRANSB = 'N'
    END IF

    IF ( M .LE. BOUND ) THEN
        IF ( BTRANS) THEN

            SELECT CASE(M)
            CASE(1)
                X(1,1) = X(1,1) / (A(1,1)**TWO-ONE)
            CASE(2)

                IF ( A(2,1) .EQ. ZERO ) THEN
                    X(2,2) = X(2,2) / (A(2,2)**TWO-ONE)
                    X(1,2) = (X(1,2) - A(1,2) * X(2,2) * A(2,2)) / (A(1,1) * A(2,2) - ONE)
                    X(2,1) = X(1,2)
                    X(1,1) = (X(1,1) - TWO*A(1,2)*X(1,2)*A(1,1) - A(1,2)*A(1,2)*X(2,2)) / (A(1,1)**TWO-ONE)

                ELSE
                    I = 0
                    DO

                        MAT(1,1) = A(1,1)*A(1,1) - ONE
                        MAT(1,2) = A(1,1)*A(1,2)
                        MAT(1,3) = A(1,2)*A(1,1)
                        MAT(1,4) = A(1,2)*A(1,2)

                        MAT(2,1) = A(1,1)*A(2,1)
                        MAT(2,2) = A(1,1)*A(2,2) - ONE
                        MAT(2,3) = A(1,2)*A(2,1)
                        MAT(2,4) = A(1,2)*A(2,2)

                        MAT(3,1) = A(2,1)*A(1,1)
                        MAT(3,2) = A(2,1)*A(1,2)
                        MAT(3,3) = A(2,2)*A(1,1) - ONE
                        MAT(3,4) = A(2,2)*A(1,2)

                        MAT(4,1) = A(2,1)*A(2,1)
                        MAT(4,2) = A(2,1)*A(2,2)
                        MAT(4,3) = A(2,2)*A(2,1)
                        MAT(4,4) = A(2,2)*A(2,2) - ONE

                        RHS(1) = X(1,1)
                        RHS(2) = X(1,2)
                        RHS(3) = X(1,2)
                        RHS(4) = X(2,2)

                        IF ( I .EQ. 0 ) THEN
                            CALL DLA_SMALL_SOLVE4( 4, MAT, RHS, IINFO, MACHINE_PRECISION, SMALL_NUMBER)
                            IF ( IINFO .EQ. 0 ) EXIT
                            I = 1
                        END IF
                        CALL DGETC2X(4, MAT, 4, IPIV, JPIV, IINFO)
                        INFO = IINFO
                        CALL DGESC2(4, MAT, 4, RHS, IPIV, JPIV, SCAL)
                        IF ( SCAL.NE.ONE) THEN
                            INFO = 1
                        END IF
                    END DO
                    X(1,1) = RHS(1)
                    X(2,1) = RHS(2)
                    X(1,2) = RHS(2)
                    X(2,2) = RHS(4)

                END IF

            END SELECT
        ELSE
           SELECT CASE(M)
            CASE(1)
                X(1,1) = X(1,1) / (A(1,1)**TWO-ONE)
            CASE(2)
                IF ( A(2,1) .EQ. ZERO ) THEN
                     X(1,1) = X(1,1) / (A(1,1)**TWO-ONE)
                     X(2,1) = (X(2,1) - A(1,2)*X(1,1)*A(1,1))/(A(1,1)*A(2,2) - ONE)
                     X(1,2) = X(2,1)
                     X(2,2) = (X(2,2) - A(1,2)**TWO*X(1,1) - TWO*A(2,2)*A(1,2)*X(2,1))/(A(2,2)**TWO-ONE)
                ELSE
                    I = 0
                    DO
                        MAT(1,1) = A(1,1)*A(1,1) - ONE
                        MAT(1,2) = A(1,1)*A(2,1)
                        MAT(1,3) = A(2,1)*A(1,1)
                        MAT(1,4) = A(2,1)*A(2,1)

                        MAT(2,1) = A(1,1)*A(1,2)
                        MAT(2,2) = A(1,1)*A(2,2) - ONE
                        MAT(2,3) = A(2,1)*A(1,2)
                        MAT(2,4) = A(2,1)*A(2,2)

                        MAT(3,1) = A(1,2)*A(1,1)
                        MAT(3,2) = A(1,2)*A(2,1)
                        MAT(3,3) = A(2,2)*A(1,1) - ONE
                        MAT(3,4) = A(2,2)*A(2,1)

                        MAT(4,1) = A(1,2)*A(1,2)
                        MAT(4,2) = A(1,2)*A(2,2)
                        MAT(4,3) = A(2,2)*A(1,2)
                        MAT(4,4) = A(2,2)*A(2,2) - ONE

                        RHS(1) = X(1,1)
                        RHS(2) = X(2,1)
                        RHS(3) = X(2,1)
                        RHS(4) = X(2,2)

                        IF ( I .EQ. 0 ) THEN
                            CALL DLA_SMALL_SOLVE4( 4, MAT, RHS, IINFO, MACHINE_PRECISION, SMALL_NUMBER)
                            IF ( IINFO .EQ. 0 ) EXIT
                            I = 1
                        END IF
                        CALL DGETC2X(4, MAT, 4, IPIV, JPIV, IINFO)
                        if ( IINFO .NE. 0 ) INFO = IINFO
                        CALL DGESC2(4, MAT, 4, RHS, IPIV, JPIV, SCAL)
                        IF ( SCAL.NE.ONE) THEN
                            INFO = 1
                        END IF
                    END DO

                    X(1,1) = RHS(1)
                    X(2,1) = RHS(3)
                    X(1,2) = RHS(3)
                    X(2,2) = RHS(4)
                END IF
            END SELECT
        END IF

        IF ( SCAL .NE. ONE) THEN
            SCALE = SCAL
        END IF
        RETURN
    END IF

    ! Prepare
    INFO = 0

    M1 = M/2
    M2 = M - M1
    IF ( A(M1+1, M1) .NE. ZERO) THEN
        M1 = M1 +1
        M2 = M2 -1
    ENDIF

    IF ( BTRANS ) THEN
        ! (N,T)
        CALL DLA_TRSTEIN_RECURSIVE(TRANS, M2, A(M1+1,M1+1), LDA, X(M1+1,M1+1), LDX, SCALE, WORK, INFO1)
        CALL DGEMM("N", "T", M2, M2, M2, ONE, X(M1+1,M1+1), LDX, A(M1+1,M1+1), LDA, ZERO, WORK, M2)
        CALL DGEMM("N", "N", M1, M2, M2, -ONE, A(1,M1+1), LDA, WORK, M2,  SCALE, X(1,M1+1), LDX)

        CALL DLA_TRSYLV2_RECURSIVE(TRANSA, TRANSB, -ONE, M1, M2, A(1,1), LDA, A(M1+1,M1+1), LDA, X(1,M1+1), LDX, &
            & SCAL, WORK, INFO1)
        IF (SCAL .NE. ONE) THEN
            SCALE = SCAL *SCALE
            X(M1+1:M,M1+1:M) = SCAL * X(M1+1:M,M1+1:M)
        END IF

        CALL DLA_ITRANSPOSE(1,M1,M1+1,M,X, LDX)


        CALL DGEMM("N", "N", M1, M2, M2, ONE, A(1,M1+1), LDA, X(M1+1,M1+1), LDX, ZERO, WORK, M1)
        CALL DGEMM("N", "T", M1, M1, M2, -ONE, WORK, M1, A(1,M1+1), LDA, SCALE, X(1,1), LDX)

        CALL DLA_DTRMM_FIX("N", M1, M2, ONE, A(1,1), LDA, X(1,M1+1), LDX, ZERO, WORK, M1)

        CALL DSYR2K("U", "N", M1, M2, -ONE, WORK, M1, A(1,M1+1), LDA, ONE, X(1,1), LDX)


        CALL DLA_TRSTEIN_RECURSIVE(TRANS, M1, A(1,1), LDA, X(1,1), LDX, SCAL, WORK, INFO1)

        IF (SCAL .NE. ONE) THEN
            SCALE = SCAL *SCALE
            X(1:M,M1+1:M) = SCAL * X(1:M,M1+1:M)
            X(M1+1:M,1:M1) = SCAL * X(M1+1:M,1:M1)
        END IF

    ELSE
        ! (T, N)
        CALL DLA_TRSTEIN_RECURSIVE(TRANS, M1, A(1,1), LDA, X(1,1), LDX, SCALE, WORK, INFO1)

        CALL DGEMM("T", "N", M2, M1, M1, ONE, A(1,M1+1), LDA, X(1,1), LDX, ZERO, WORK, M2)
        ! CALL DGEMM("N", "N", M2, M1, M1, -ONE, WORK, M2, A(1,1), LDA, SCALE, X(M1+1,1), LDX)
        DO K = 1, M1, 64
            KH = MIN(M1, K + 64  - 1)
            KB = KH - K + 1

            IF ( K .LT. M1 ) THEN
                IF ( A(KH+1,KH) .NE. ZERO ) THEN
                    CALL DGEMM("N", "N", M2, KB, KH+1, -ONE, WORK, M2, A(1,K), LDA, SCALE, X(M1+1,K), LDX)
                ELSE
                    CALL DGEMM("N", "N", M2, KB, KH, -ONE, WORK, M2, A(1,K), LDA, SCALE, X(M1+1,K), LDX)

                END IF
            ELSE
                CALL DGEMM("N", "N", M2, KB, KH, -ONE, WORK, M2, A(1,K), LDA, SCALE, X(M1+1,K), LDX)
            END IF
        END DO

        CALL DLA_TRSYLV2_RECURSIVE(TRANSA, TRANSB, -ONE, M2, M1, A(M1+1,M1+1), LDA, A(1,1), LDA, X(M1+1,1), LDX, &
            & SCAL, WORK, INFO1 )
        IF (SCAL .NE. ONE) THEN
            SCALE = SCAL *SCALE
            X(1:M1,1:M1) = SCAL * X(1:M1,1:M1)
        END IF
        CALL DLA_ITRANSPOSE(M1+1,M, 1, M1, X, LDX)

        CALL DGEMM("T", "N", M2, M1, M1, ONE, A(1,M1+1), LDA, X(1,1), LDX, ZERO, WORK, M2)
        CALL DGEMM("N", "N", M2, M2, M1, -ONE, WORK, M2, A(1,M1+1), LDA, SCALE, X(M1+1,M1+1), LDX)

        CALL DGEMM("N", "N", M1, M2, M2, ONE, X(1,M1+1), LDX, A(M1+1,M1+1), LDA, ZERO, WORK, M1)
        ! DO K = 1, M2, 64
        !     KH = MIN(M2, K + 64  - 1)
        !     KB = KH - K + 1
        !     IF ( K .LT. M2 ) THEN
        !         IF ( A(M1+KH+1,M1+KH) .NE. ZERO ) THEN
        !             CALL DGEMM("N", "N", M1, KB, KH+1, ONE, X(1,M1+1), LDX, A(M1+1,M1+K), LDA, ZERO, WORK((K-1)*M1+1), M1)
        !         ELSE
        !             CALL DGEMM("N", "N", M1, KB, KH, ONE, X(1,M1+1), LDX, A(M1+1,M1+K), LDA, ZERO, WORK((K-1)*M1+1), M1)
        !         END IF
        !
        !     ELSE
        !         CALL DGEMM("N", "N", M1, KB, KH, ONE, X(1,M1+1), LDX, A(M1+1,M1+K), LDA, ZERO, WORK((K-1)*M1+1), M1)
        !     END IF
        ! END DO
        CALL DSYR2K("L","T", M2, M1, -ONE, A(1,M1+1), LDA, WORK, M1, ONE, X(M1+1,M1+1), LDX)


        CALL DLA_TRSTEIN_RECURSIVE(TRANS, M2, A(M1+1,M1+1), LDA, X(M1+1,M1+1), LDX, SCALE, WORK, INFO1)

        IF (SCAL .NE. ONE) THEN
            SCALE = SCAL *SCALE
            X(1:M,1:M1) = SCAL * X(1:M,1:M1)
            X(1:M1,M1+1:M) = SCAL * X(1:M1,M1+1:M)
        END IF

    END IF
END SUBROUTINE



