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

!> \brief Solver for a 4x4 standard Lyapunov equation (TRANS = N)
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TRLYAP_KERNEL_44N (M, A, LDA, X, LDX, SCALE, INFO)
!           REAL SCALE
!           INTEGER M, LDA, LDX, INFO
!           REAL A(LDA, *), X(LDX, *)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TRLYAP_KERNEL_44N solves a Lyapunov  equation of the following form
!>
!>    A * X  + X * A**T = SCALE * Y                              (1)
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
!> \date Januar 2023
!> \ingroup sgltrlyap
!
SUBROUTINE SLA_TRLYAP_KERNEL_44N (M, A, LDA, X, LDX, SCALE, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
    ! Arguments
    REAL SCALE
    INTEGER M, LDA, LDX, INFO
    REAL A(LDA, *), X(LDX, *)


    ! Local Variables
    INTEGER K, L, KB, LB, KH, LH, KE, IT, LE
    INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
    REAL MAT(4,4), RHS(4)
    REAL TMP(4)
    REAL A11, A12, A21, A22
    REAL B11, B12, B21, B22
    REAL SCAL

    REAL ZERO , ONE, HALF
    PARAMETER (ZERO = 0.0, ONE = 1.0, HALF=0.5)
    INTEGER IONE
    PARAMETER (IONE = 1)

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

    ! (N,T)
    L = M
    DO WHILE ( L .GT. 0 )
        LB = 1
        LH = L
        IF ( L .GT. 1 ) THEN
            IF (A(L,L-1) .NE. ZERO) THEN
                LB = 2
                L = L - 1
            END IF
        END IF
        LE = L - 1

        K = LH
        DO WHILE (K .GT. 0)
            KB = 1
            KH = K
            IF ( K .GT. 1 ) THEN
                IF (A(K,K-1) .NE. ZERO ) THEN
                    KB = 2
                    K = K - 1
                END IF
            END IF
            KE = K - 1

            IT = 0

            IF (K .EQ. L .AND. KB .EQ. 2 ) THEN
                X(KH, K ) = X(K,KH)
            END IF

            DO

                ! Solve Inner System
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 1

                    MAT(1,1) = A(K,K) + A(L,L)
                    RHS(1)   = X(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 2

                    A11 = A(K,K)

                    B11 = A(L,L)
                    B21 = A(LH,L)
                    B12 = A(L,LH)
                    B22 = A(LH,LH)


                    MAT(1,1) = A11 + B11
                    MAT(1,2) =       B12
                    MAT(2,1) =       B21
                    MAT(2,2) = A11 + B22

                    RHS(1) = X(K,L)
                    RHS(2) = X(K,LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2

                    A11 = A(K,K)
                    A12 = A(K,KH)
                    A21 = A(KH,K)
                    A22 = A(KH,KH)

                    B11 = A(L,L)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A12
                    MAT(2,1) = A21
                    MAT(2,2) = A22 + B11

                    RHS(1) = X(K,L)
                    RHS(2) = X(KH,L)

                ELSE
                    DIMMAT = 4

                    A11 = A(K,K)
                    A12 = A(K,KH)
                    A21 = A(KH,K)
                    A22 = A(KH,KH)

                    B11 = A(L,L)
                    B12 = A(L,LH)
                    B21 = A(LH, L)
                    B22 = A(LH,LH)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A12
                    MAT(1,3) =       B12
                    MAT(1,4) = ZERO

                    MAT(2,1) = A21
                    MAT(2,2) = A22 + B11
                    MAT(2,3) = ZERO
                    MAT(2,4) =     + B12

                    MAT(3,1) = B21
                    MAT(3,2) = ZERO
                    MAT(3,3) = A11 + B22
                    MAT(3,4) = A12

                    MAT(4,1) = ZERO
                    MAT(4,2) = B21
                    MAT(4,3) = A21
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

            IF ( K .EQ. L ) THEN
                ! Diagonal Block
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    X(K,L) = RHS(1)
                ELSE
                    X(K,L) = RHS(1)
                    X(KH,L) = HALF*(RHS(2)+RHS(3))
                    X(K, LH) = X(KH,L)
                    X(KH,LH) = RHS(4)
                END IF
            ELSE
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    X(K,L) = RHS(1)
                    X(L,K) = RHS(1)
                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    X(K,L) = RHS(1)
                    X(K,LH) = RHS(2)
                    X(L,K) = RHS(1)
                    X(LH,K) = RHS(2)
                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    X(K,L) = RHS(1)
                    X(KH,L) = RHS(2)
                    X(L, K) = RHS(1)
                    X(L, KH) = RHS(2)
                ELSE
                    X(K,L) = RHS(1)
                    X(KH,L) = RHS(2)
                    X(K, LH) = RHS(3)
                    X(KH,LH) = RHS(4)

                    X(L,K) =   RHS(1)
                    X(L, KH) = RHS(2)
                    X(LH, K) = RHS(3)
                    X(LH, KH) = RHS(4)
                END IF

            END IF

            IF ( K .GT. 1 ) THEN
                IF ( KB .EQ. IONE .AND. LB .EQ. IONE ) THEN
                    X(1:KE, L ) = X(1:KE, L) - X(K,L)* A(1:KE, K)
                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    X(1:KE,L) = X(1:KE,L) - X(K,L) * A(1:KE,K)
                    X(1:KE,LH) = X(1:KE,LH) - X(K,LH) *A(1:KE,K)
                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    X(1:KE, L) = X(1:KE,L) - A(1:KE,K) *X(K,L) - A(1:KE,KH) *X(KH,L)
                ELSE
                    X(1:KE,L) = X(1:KE,L) -   X(K,L)  * A(1:KE,K) - X(KH,L) * A(1:KE,KH)
                    X(1:KE,LH) = X(1:KE,LH) - X(K,LH) * A(1:KE,K) - X(KH,LH) * A(1:KE,KH)
                END IF

            END IF

            K = K - 1
        END DO

        ! Update January 2021
        IF ( L .GT. 1 ) THEN
            IF ( LB .EQ. IONE ) THEN
                ! CALL SSYR2("U", L-1, -ONE, X(1,L), IONE, A(1,L), IONE, X(1,1), LDX)

                SELECT CASE (L-1)
                CASE (1)
                    X(1,1) = X(1,1) - 2 * X(1,L) * A(1,L)
                case (2)
                    ! X(1,1:2) = X(1,1:2) - X(1,L) * A(1:2,L) - A(1,L) * X(1:2,L)
                    TMP(1:2) = X(1,1:2) - X(1,L) * A(1:2,L) - A(1,L) * X(1:2,L)
                    X(1,1:2) = TMP(1:2)
                    X(2,2)   = X(2,2)   - X(2,L) * A(2,L)   - A(2,L) * X(2,L)
                case (3)
                    TMP(1:3) = X(1,1:3) - X(1,L) * A(1:3,L) - A(1,L) * X(1:3,L)
                    X(1,1:3) = TMP(1:3)
                    TMP(2:3) = X(2,2:3) - X(2,L) * A(2:3,L) - A(2,L) * X(2:3,L)
                    X(2,2:3) = TMP(2:3)
                    X(3,3)   = X(3,3)   - X(3,L) * A(3,L)   - A(3,L) * X(3,L)
                END SELECT
            ELSE
                ! CALL SSYR2("U", L-1, -ONE, X(1,L), IONE, A(1,L), IONE, X(1,1), LDX)
                ! CALL SSYR2("U", L-1, -ONE, X(1,LH), IONE, A(1,LH), IONE, X(1,1), LDX)

                SELECT CASE (L-1)
                CASE (1)
                    X(1,1) = X(1,1) - 2 * X(1,L) * A(1,L) - 2 * X(1,LH) * A(1,LH)
                case (2)
                    TMP(1:2) = X(1,1:2) - X(1,L) * A(1:2,L) - A(1,L) * X(1:2,L) - X(1,LH) * A(1:2,LH) - A(1,LH) * X(1:2,LH)
                    X(1,1:2) = TMP(1:2)
                    X(2,2)   = X(2,2)   - X(2,L) * A(2,L)   - A(2,L) * X(2,L)   - X(2,LH) * A(2,LH)   - A(2,LH) * X(2,LH)
                case (3)
                    TMP(1:3) = X(1,1:3) - X(1,L) * A(1:3,L) - A(1,L) * X(1:3,L) - X(1,LH) * A(1:3,LH) - A(1,LH) * X(1:3,LH)
                    X(2,1:3) = TMP(1:3)
                    TMP(2:3) = X(2,2:3) - X(2,L) * A(2:3,L) - A(2,L) * X(2:3,L) - X(2,LH) * A(2:3,LH) - A(2,LH) * X(2:3,LH)
                    X(2,2:3) = TMP(2:3)

                    X(3,3)   = X(3,3)   - X(3,L) * A(3,L)   - A(3,L) * X(3,L)   - X(3,LH) * A(3,LH)   - A(3,LH) * X(3,LH)
                END SELECT

            ENDIF
        END IF
        L = L - 1
    END DO

END SUBROUTINE

