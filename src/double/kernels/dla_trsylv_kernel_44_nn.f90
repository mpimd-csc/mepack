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


!> \brief Solver for a 4x4 Sylvester equation (TRANSA = N, TRANSB = N)
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TRSYLV_KERNEL_44NN (SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, INFO)
!           DOUBLE PRECISION SGN, SCALE
!           INTEGER M, N, LDA, LDB, LDX, INFO
!           DOUBLE PRECISION A(LDA, *), B(LDB, *),  X(LDX, *)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLA_TRSYLV_KERNEL_44NN solves a Sylvester equation of the following form
!>
!>    A * X  + SGN * X * B = SCALE * Y                              (1)
!>
!> where A is a M-by-M quasi upper triangular matrix and B is N-by-N quasi
!> upper triangular matrices. The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrices A and B are create by DGEES from LAPACK.
!> The algorithm is implemented without using BLAS level 2
!> operations. Thereby the order of M and N is at most 4. Furthermore, for fast execution
!> the function does not check the input arguments.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SGN
!> \verbatim
!>          SGN is DOUBLE PRECISION, allowed values: +/-1
!>          Specifies the sign between the two parts of the Sylvester equation.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and C.  4 >= M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices B and D.  4 >= N >= 0.
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
!>          The matrix B must be (quasi-) upper triangular. If the matrix D is already
!>          quasi-upper triangular the matrix B must be upper triangular.
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
!> \ingroup dbltrsylv
!
SUBROUTINE DLA_TRSYLV_KERNEL_44NN (SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Arguments
    DOUBLE PRECISION SGN, SCALE
    INTEGER M, N, LDA, LDB,LDX, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *),  X(LDX, *)


    ! Local Variables
    INTEGER K, L, KB, LB, KH, LH, KE, IT
    INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
    DOUBLE PRECISION MAT(4,4), RHS(4)
    DOUBLE PRECISION A11, A12, A21, A22
    DOUBLE PRECISION B11, B12, B21, B22
    DOUBLE PRECISION SCAL

    DOUBLE PRECISION ZERO , ONE
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
    INTEGER IONE
    PARAMETER (IONE = 1)

    ! EXTERNAL FUNCTIONS
    EXTERNAL DLA_SMALL_SOLVE4
    EXTERNAL DGETC2X
    EXTERNAL DGESC2

    ! Quick Return
    SCALE = ONE
    IF ( M .EQ. 0 .OR. N .EQ. 0 ) THEN
        RETURN
    END IF

    ! Prepare
    INFO = 0

    ! ( N, N )
    K = M
    DO WHILE (K .GT. 0)
        KB = 1
        KH = K
        IF ( K .GT. 1 ) THEN
            IF ( A(K,K-1) .NE. ZERO ) THEN
                KB = 2
                K = K - 1
            END IF
        END IF
        KE = K - 1


        L = 1
        DO WHILE ( L .LE. N )
            LB = 1
            IF ( L .LT. N ) THEN
                IF ( B(L+1,L) .NE. ZERO ) THEN
                    LB = 2
                END IF
            END IF
            LH = L + LB - 1

            IT = 0
            DO
                ! Solve inner system
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 1

                    MAT(1,1) = A(K,K) + SGN * B(L,L)
                    RHS(1)   = X(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 2

                    A11 = A(K,K)

                    B11 = SGN *B(L,L)
                    B12 = SGN *B(L,LH)
                    B21 = SGN *B(LH,L)
                    B22 = SGN *B(LH,LH)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) =       B21
                    MAT(2,1) =       B12
                    MAT(2,2) = A11 + B22

                    RHS(1) = X(K,L)
                    RHS(2) = X(K,LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2

                    A11 = A(K,K)
                    A12 = A(K,KH)
                    A21 = A(KH,K)
                    A22 = A(KH,KH)

                    B11 = SGN*B(L,L)

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

                    B11 = SGN* B(L,L)
                    B12 = SGN* B(L,LH)
                    B21 = SGN* B(LH, L)
                    B22 = SGN* B(LH,LH)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A12
                    MAT(1,3) =       B21
                    MAT(1,4) = ZERO

                    MAT(2,1) = A21
                    MAT(2,2) = A22 + B11
                    MAT(2,3) = ZERO
                    MAT(2,4) =       B21

                    MAT(3,1) =       B12
                    MAT(3,2) = ZERO
                    MAT(3,3) = A11 + B22
                    MAT(3,4) = A12

                    MAT(4,1) = ZERO
                    MAT(4,2) =       B12
                    MAT(4,3) = A21
                    MAT(4,4) = A22 + B22

                    RHS(1) = X(K,L)
                    RHS(2) = X(KH,L)
                    RHS(3) = X(K, LH)
                    RHS(4) = X(KH,LH)

                END IF

                IF ( IT .EQ. 0 ) THEN
                    CALL DLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL DGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL DGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
                    IF ( SCAL .NE. ONE ) THEN
                        X(1:M, 1:N) = SCAL*X(1:M, 1:N)
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

            ! Update January 2021
            ! in current row
            IF ( LH .LT. N ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    X(K,LH+1:N) = X(K,LH+1:N) - SGN * X(K,L) * B(L,LH+1:N)
                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    X(K,LH+1:N) = X(K,LH+1:N) - SGN*X(K,L) * B(L,LH+1:N) - SGN* X(K,LH) * B(LH,LH+1:N)
                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    X(K,LH+1:N)  = X(K,LH+1:N)  - SGN * X(K,L)  * B(L,LH+1:N)
                    X(KH,LH+1:N) = X(KH,LH+1:N) - SGN * X(KH,L) * B(L,LH+1:N)
                ELSE
                    X(K,LH+1:N) = X(K,LH+1:N) - SGN * X(K,L) * B(L,LH+1:N) - SGN * X(K,LH) * B(LH, LH+1:N)
                    X(KH,LH+1:N) = X(KH,LH+1:N) - SGN * X(KH,L) * B(L,LH+1:N) - SGN * X(KH,LH) * B(LH, LH+1:N)
                END IF
            END IF

            L = LH + 1
        END DO

        ! Update January 2021
        IF ( K .GT. 1 ) THEN
            IF ( KB .EQ. 1 ) THEN
                SELECT CASE(N)
                CASE(1)
                    X(1:KE,1) = X(1:KE, 1) -       A(1:KE,K) * X(K,1)
                CASE(2)
                    X(1:KE,1) = X(1:KE, 1) -       A(1:KE,K) * X(K,1)
                    X(1:KE,2) = X(1:KE, 2) -       A(1:KE,K) * X(K,2)

                CASE(3)
                    X(1:KE,1) = X(1:KE, 1) -       A(1:KE,K) * X(K,1)
                    X(1:KE,2) = X(1:KE, 2) -       A(1:KE,K) * X(K,2)
                    X(1:KE,3) = X(1:KE, 3) -       A(1:KE,K) * X(K,3)

                CASE(4)
                    X(1:KE,1) = X(1:KE, 1) -       A(1:KE,K) * X(K,1)
                    X(1:KE,2) = X(1:KE, 2) -       A(1:KE,K) * X(K,2)
                    X(1:KE,3) = X(1:KE, 3) -       A(1:KE,K) * X(K,3)
                    X(1:KE,4) = X(1:KE, 4) -       A(1:KE,K) * X(K,4)


                END SELECT
            ELSE
                SELECT CASE(N)
                CASE(1)
                    X(1:KE,1) = X(1:KE, 1) -       (A(1:KE,K) * X(K,1) + A(1:KE,KH) * X(KH,1))

                CASE(2)
                    X(1:KE,1) = X(1:KE, 1) -       (A(1:KE,K) * X(K,1) + A(1:KE,KH) * X(KH,1))
                    X(1:KE,2) = X(1:KE, 2) -       (A(1:KE,K) * X(K,2) + A(1:KE,KH) * X(KH,2))

                CASE(3)
                    X(1:KE,1) = X(1:KE, 1) -       (A(1:KE,K) * X(K,1) + A(1:KE,KH) * X(KH,1))
                    X(1:KE,2) = X(1:KE, 2) -       (A(1:KE,K) * X(K,2) + A(1:KE,KH) * X(KH,2))
                    X(1:KE,3) = X(1:KE, 3) -       (A(1:KE,K) * X(K,3) + A(1:KE,KH) * X(KH,3))


                CASE(4)
                    X(1:KE,1) = X(1:KE, 1) -       (A(1:KE,K) * X(K,1) + A(1:KE,KH) * X(KH,1))
                    X(1:KE,2) = X(1:KE, 2) -       (A(1:KE,K) * X(K,2) + A(1:KE,KH) * X(KH,2))
                    X(1:KE,3) = X(1:KE, 3) -       (A(1:KE,K) * X(K,3) + A(1:KE,KH) * X(KH,3))
                    X(1:KE,4) = X(1:KE, 4) -       (A(1:KE,K) * X(K,4) + A(1:KE,KH) * X(KH,4))
                END SELECT
            END IF
        END IF

        K = K - 1
    END DO

END SUBROUTINE



