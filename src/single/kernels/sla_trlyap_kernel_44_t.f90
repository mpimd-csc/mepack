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

!> \brief Solver for a 4x4 standard Lyapunov equation (TRANS = T)
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TRLYAP_KERNEL_44T (M, A, LDA, X, LDX, SCALE, INFO)
!           REAL SCALE
!           INTEGER M, LDA, LDX, INFO
!           REAL A(LDA, *), X(LDX, *)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TRLYAP_KERNEL_44T solves a Lyapunov  equation of the following form
!>
!>    A **T * X  + X * A = SCALE * Y                              (1)
!>
!> where A is a M-by-M quasi upper triangular matrix.
!> The right hand side Y and the solution X M-by-M matrices.
!> Typically the matrix  A is created by SGEES from LAPACK.
!> The algorithm is implemented without BLAS level 2
!> operations. Thereby the order of M and N is at most 4. Furthermore, for fast execution
!> the function does not check the input arguments.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and C.  4 >= M >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,M)
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
!>          X is REAL array, dimension (LDX,N)
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
!>          SCALE is REAL
!>          SCALE is a scaling factor to prevent the overflow in the result.
!>          If INFO == 0 then SCALE is 1.0 otherwise if one of the inner systems
!>          could not be solved correctly, 0 < SCAL <= 1 holds true.
!> \endverbatim
!>
!> \param[in,out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          On output:
!>          == 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  The equation is not solved correctly. One of the arising inner
!>                system got singular.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date October 2023
!> \ingroup sgltrlyap
!
SUBROUTINE SLA_TRLYAP_KERNEL_44T (M, A, LDA, X, LDX, SCALE, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
    ! Arguments
    REAL SCALE
    INTEGER M, LDA, LDX, INFO
    REAL A(LDA, *), X(LDX, *)


    ! Local Variables
    INTEGER K, L, KB, LB, KH, LH, IT, II, JJ
    INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
    REAL MAT(4,4), RHS(4)
    REAL A11, A12, A21, A22
    REAL B11, B12, B21, B22
    REAL SCAL
    REAL AL(4,4)
    REAL TMP(4)

    REAL ZERO , ONE, HALF
    PARAMETER (ZERO = 0.0, ONE = 1.0, HALF=0.5)
    INTEGER IONE
    PARAMETER (IONE = 1)

    !dir$ attributes align: 64:: AL

    ! EXTERNAL FUNCTIONS
    EXTERNAL SLA_SMALL_SOLVE4
    EXTERNAL SGETC2X
    EXTERNAL SGESC2


    ! Quick Return
    SCALE = ONE
    IF ( M .EQ. 0  ) THEN
        RETURN
    END IF
    IF ( M .GT. 4 ) THEN
        INFO = -1
        RETURN
    END IF

    ! Prepare
    INFO = 0

    AL = TRANSPOSE(A(1:M,1:M))
    ! (N,T)
    K = 1
    DO WHILE ( K .LE. M )
        KB = 1
        IF ( K .LT. M ) THEN
            IF (AL(K,K+1) .NE. ZERO ) THEN
                KB = 2
            END IF
        END IF
        KH = K + KB - 1


        L = K
        DO WHILE ( L .LE. M )
            LB = 1
            IF ( L .LT. M ) THEN
                IF ( ( AL(L,L+1) .NE. ZERO ) ) THEN
                    LB = 2
                END IF
            END IF
            LH = L + LB - 1


            IF (K.EQ.L .AND. KB.EQ.2) THEN
                X(KH,K) = X(K,KH)
            END IF

            IT = 0
            DO

                ! Solve Inner System
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 1

                    MAT(1,1) = AL(K,K) + AL(L,L)
                    RHS(1)   = X(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 2

                    A11 = AL(K,K)


                    B11 = AL(L,L)
                    B21 = AL(L,LH)
                    B12 = AL(LH,L)
                    B22 = AL(LH,LH)


                    MAT(1,1) = A11 + B11
                    MAT(1,2) =       B21
                    MAT(2,1) =       B12
                    MAT(2,2) = A11 + B22

                    RHS(1) = X(K,L)
                    RHS(2) = X(K,LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2

                    A11 = AL(K,K)
                    A12 = AL(KH,K)
                    A21 = AL(K,KH)
                    A22 = AL(KH,KH)

                    B11 = AL(L,L)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A21
                    MAT(2,1) = A12
                    MAT(2,2) = A22 + B11

                    RHS(1) = X(K,L)
                    RHS(2) = X(KH,L)

                ELSE
                    DIMMAT = 4

                    A11 = AL(K,K)
                    A12 = AL(KH,K)
                    A21 = AL(K,KH)
                    A22 = AL(KH,KH)


                    B11 = AL(L,L)
                    B12 = AL(LH,L)
                    B21 = A(LH, L)
                    B22 = AL(LH,LH)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A21
                    MAT(1,3) =       B21
                    MAT(1,4) = ZERO

                    MAT(2,1) = A12
                    MAT(2,2) = A22 + B11
                    MAT(2,3) = ZERO
                    MAT(2,4) =       B21

                    MAT(3,1) =        B12
                    MAT(3,2) = ZERO
                    MAT(3,3) = A11 +  B22
                    MAT(3,4) = A21

                    MAT(4,1) = ZERO
                    MAT(4,2) =        B12
                    MAT(4,3) = A12
                    MAT(4,4) = A22 + B22

                    RHS(1) = X(K,L)
                    RHS(2) = X(KH,L)
                    RHS(3) = X(K, LH)
                    RHS(4) = X(KH,LH)

                END IF
                IF ( IT .EQ. 0 ) THEN
                    CALL SLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL SGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL SGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
                    IF ( SCAL .NE. ONE ) THEN
                        X(1:M, 1:M) = SCAL*X(1:M, 1:M)
                        SCALE = SCALE * SCAL
                        INFO = 1
                    END IF
                    EXIT
                END IF

            END DO

            IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                X(K,L) = RHS(1)
            ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                X(K,L) = RHS(1)
                X(K,LH) = RHS(2)
            ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                X(K,L) = RHS(1)
                X(KH,L) = RHS(2)
            ELSE
                X(K,L) = RHS(1)
                X(KH,L) = RHS(2)
                X(K, LH) = RHS(3)
                X(KH,LH) = RHS(4)
            END IF

            IF ( K.NE.L) THEN
                DO II = K, KH
                    DO JJ = II+1,LH
                        TMP(1) = X(JJ,II)
                        X(JJ,II) = X(II, JJ)
                        X(II,JJ) = TMP(1)
                    END DO
                END DO
                ! The transpose operation produces an temporary copy

                ! X(L:LH,K:KH) = TRANSPOSE(X(K:KH,L:LH))
            END IF

            ! in current row
            IF ( LH .LT. M ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    X(K,LH+1:M) = X(K,LH+1:M) - X(K,L) * AL(LH+1:M,L)
                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    X(K,LH+1:M) = X(K,LH+1:M) - X(K,L) * AL(LH+1:M,L) - X(K,LH) * AL(LH+1:M,LH)
                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    X(K,LH+1:M) = X(K,LH+1:M) - X(K,L) * AL(LH+1:M,L)
                    X(KH,LH+1:M) = X(KH,LH+1:M) - X(KH,L) * AL(LH+1:M,L)
                ELSE
                    X(K,LH+1:M) = X(K,LH+1:M) - X(K,L) * AL(LH+1:M,L) - X(K,LH) * AL(LH+1:M,LH)
                    X(KH,LH+1:M) = X(KH,LH+1:M) - X(KH,L) * AL(LH+1:M,L) - X(KH,LH) * AL(LH+1:M,LH)
                END IF
            END IF


            L = LH + 1
        END DO

        IF ( KH .LT. M ) THEN
            IF ( KB.EQ.IONE) THEN
                SELECT CASE(M-KH)
                CASE(1)
                    X(KH+1,KH+1) = X(KH+1,KH+1) - 2* AL(KH+1,K) * X(K,KH+1)
                CASE(2)
                    X(KH+1,KH+1:M) = X(KH+1,KH+1:M) - AL(KH+1,K) * X(K,KH+1:M) - X(K,KH+1) * AL(KH+1:M,K)
                    X(KH+2,KH+2:M) = X(KH+2,KH+2:M) - AL(KH+2,K) * X(K,KH+2) - X(K,KH+2) * AL(KH+2,K)
                CASE(3)
                    X(KH+1,KH+1:M) = X(KH+1,KH+1:M) - AL(KH+1,K) * X(K,KH+1:M) - X(K,KH+1) * AL(KH+1:M,K)
                    X(KH+2,KH+2:M) = X(KH+2,KH+2:M) - AL(KH+2,K) * X(K,KH+2:M) - X(K,KH+2) * AL(KH+2:M,K)
                    X(KH+3,KH+3:M) = X(KH+3,KH+3:M) - AL(KH+3,K) * X(K,KH+3) - X(K,KH+3) * AL(KH+3,K)
                END SELECT

            ELSE
                SELECT CASE(M-KH)
                CASE(1)
                    X(KH+1,KH+1) = X(KH+1,KH+1) - 2* AL(KH+1,K) * X(K,KH+1) - 2 * AL(KH+1,KH) * X(KH,KH+1)
                CASE(2)
                    X(KH+1,KH+1:M) = X(KH+1,KH+1:M) - AL(KH+1,K) * X(K,KH+1:M) - X(K,KH+1) * AL(KH+1:M,K) &
                        & - AL(KH+1,KH) * X(KH,KH+1:M) - X(KH,KH+1) * AL(KH+1:M,KH)
                    X(KH+2,KH+2:M) = X(KH+2,KH+2:M) - AL(KH+2,K) * X(K,KH+2) - X(K,KH+2) * AL(KH+2,K) &
                        & - AL(KH+2,KH) * X(KH,KH+2) - X(KH,KH+2) * AL(KH+2,KH)
                CASE(3)
                    X(KH+1,KH+1:M) = X(KH+1,KH+1:M) - AL(KH+1,K) * X(K,KH+1:M) - X(K,KH+1) * AL(KH+1:M,K) &
                        & - AL(KH+1,KH) * X(KH,KH+1:M) - X(KH,KH+1) * AL(KH+1:M,KH)
                    X(KH+2,KH+2:M) = X(KH+2,KH+2:M) - AL(KH+2,K) * X(K,KH+2:M) - X(K,KH+2) * AL(KH+2:M,K) &
                        & - AL(KH+2,KH) * X(KH,KH+2:M) - X(KH,KH+2) * AL(KH+2:M,KH)
                    X(KH+3,KH+3:M) = X(KH+3,KH+3:M) - AL(KH+3,K) * X(K,KH+3) - X(K,KH+3) * AL(KH+3,K) &
                        & - AL(KH+3,KH) * X(KH,KH+3) - X(KH,KH+3) * AL(KH+3,KH)
                END SELECT
            ENDIF
        END IF

        K = KH + 1
    END DO

END SUBROUTINE

