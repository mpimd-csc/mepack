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

!> \brief Solver for a 4x4 discrete time Sylvester equation (TRANSA = N, TRANSB = N)
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TRSYLV2_KERNEL_44NN (SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, INFO)
!           REAL SGN, SCALE
!           INTEGER M, N, LDA, LDB, LDX, INFO
!           REAL A(LDA, *), B(LDB, *) , X(LDX, *)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TRSYLV2_KERNEL_44NN solves a discrete time Sylvester equation of the following form
!>
!>    A * X * B + SGN * X  = SCALE * Y                              (1)
!>
!> where A is a M-by-M quasi upper triangular matrix and B is a N-by-N quasi upper
!> triangular matrix.  The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrices A and B are created by SGEES from LAPACK.
!> The algorithm is implemented using BLAS level 2
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
!>          SGN is REAL, allowed values: +/-1
!>          Specifies the sign between the two parts of the Sylvester equation.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrix A.  4 >= M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B.  4 >= N >= 0.
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
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
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
!>          X is REAL array, dimension (LDX,N)
!>          On input, the matrix X contains the right hand side Y.
!>          On output, the matrix X contains the solution of Equation (1)
!>          as selected by SGN.
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
!>          On input:
!>          = 1 :  Skip the input data checks
!>          <> 1:  Check input data like normal LAPACK like routines.
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
!> \ingroup sgltrsylv
!
SUBROUTINE SLA_TRSYLV2_KERNEL_44NN (SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
    ! Arguments
    REAL SGN, SCALE
    INTEGER M, N, LDA, LDB, LDX, INFO
    REAL A(LDA, *), B(LDB, *) , X(LDX, *)


    ! Local Variables
    INTEGER K, L, KB, LB, KH, LH, KE, IT
    INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
    REAL MAT(4,4), RHS(4)
    REAL A11, A12, A21, A22
    REAL B11, B12, B21, B22
    REAL SCAL

    REAL WORK(4,4)

    REAL ZERO , ONE
    PARAMETER (ZERO = 0.0, ONE = 1.0)

    ! EXTERNAL FUNCTIONS
    EXTERNAL SLA_SMALL_SOLVE4
    EXTERNAL SGETC2X
    EXTERNAL SGESC2


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
                IF (( B(L+1,L) .NE. ZERO ) ) THEN
                    LB = 2
                END IF
            END IF
            LH = L + LB - 1

            IT = 0
            DO
                ! Solve inner system
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 1

                    MAT(1,1) = B(L,L) * A(K,K) + SGN
                    RHS(1)   = X(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 2

                    A11 = A(K,K)

                    B11 = B(L,L)
                    B21 = B(LH,L)
                    B12 = B(L,LH)
                    B22 = B(LH,LH)


                    MAT(1,1) = B11*A11 + SGN
                    MAT(1,2) = B21*A11
                    MAT(2,1) = B12*A11
                    MAT(2,2) = B22*A11 + SGN

                    RHS(1) = X(K,L)
                    RHS(2) = X(K,LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2

                    A11 = A(K,K)
                    A12 = A(K,KH)
                    A21 = A(KH,K)
                    A22 = A(KH,KH)
                    B11 = B(L,L)

                    MAT(1,1) = B11*A11 + SGN
                    MAT(1,2) = B11*A12
                    MAT(2,1) = B11*A21
                    MAT(2,2) = B11*A22 + SGN

                    RHS(1) = X(K,L)
                    RHS(2) = X(KH,L)

                ELSE
                    DIMMAT = 4

                    A11 = A(K,K)
                    A12 = A(K,KH)
                    A21 = A(KH,K)
                    A22 = A(KH,KH)

                    B11 = B(L,L)
                    B12 = B(L,LH)
                    B21 = B(LH, L)
                    B22 = B(LH,LH)

                    MAT(1,1) = B11*A11 + SGN
                    MAT(1,2) = B11*A12
                    MAT(1,3) = B21*A11
                    MAT(1,4) = B21*A12

                    MAT(2,1) = B11*A21
                    MAT(2,2) = B11*A22 + SGN
                    MAT(2,3) = B21*A21
                    MAT(2,4) = B21*A22

                    MAT(3,1) = B12*A11
                    MAT(3,2) = B12*A12
                    MAT(3,3) = B22*A11 + SGN
                    MAT(3,4) = B22*A12

                    MAT(4,1) = B12*A21
                    MAT(4,2) = B12*A22
                    MAT(4,3) = B22*A21
                    MAT(4,4) = B22*A22 + SGN

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
                IF ( KB .EQ. 1 .AND. LB .EQ. 1) THEN
                    X(K, LH+1:N) = X(K, LH+1:N) - A(K,K) * X(K,L) * B(L,LH+1:N)
                ELSE IF (KB .EQ. 1 .AND. LB .EQ. 2) THEN
                    X(K, LH+1:N) = X(K, LH+1:N) - A(K,K) * (X(K,L) * B(L,LH+1:N) + X(K,LH) * B(LH,LH+1:N))

                ELSE IF (KB .EQ. 2 .AND. LB .EQ. 1) THEN
                    X(K, LH+1:N) = X(K, LH+1:N) - (A(K,K) * X(K,L) + A(K,KH) * X(KH,L)) * B(L,LH+1:N)
                    X(KH, LH+1:N) = X(KH, LH+1:N) - (A(KH,K) * X(K,L) + A(KH,KH) * X(KH,L)) * B(L,LH+1:N)
                ELSE
                    WORK(1,1) = A(K,K) * X(K,L) + A(K,KH) * X(KH,L)
                    WORK(1,2) = A(K,K) * X(K,LH) + A(K,KH) * X(KH,LH)
                    WORK(2,1) = A(KH,K) * X(K,L) + A(KH,KH) * X(KH,L)
                    WORK(2,2) = A(KH,K) * X(K,LH) + A(KH,KH) * X(KH,LH)

                    X(K, LH+1:N) = X(K, LH+1:N) - (WORK(1,1) *  B(L,LH+1:N) + WORK(1,2) * B(LH,LH+1:N))
                    X(KH, LH+1:N) = X(KH, LH+1:N) - (WORK(2,1) *  B(L,LH+1:N) + WORK(2,2) * B(LH,LH+1:N))
                ENDIF

            END IF

            L = LH + 1
        END DO

        ! Update January 2021
        IF ( K .GT. 1 ) THEN
            IF ( KB .EQ. 1 ) THEN
                SELECT CASE(N)
                CASE(1)
                    WORK(1,1) = X(K,1) * B(1,1)
                    X(1:KE,1) = X(1:KE, 1) -       A(1:KE,K) * WORK(1,1)

                CASE(2)
                    WORK(1,1) = X(K,1) * B(1,1) + X(K,2) * B(2,1)
                    WORK(2,1) = X(K,1) * B(1,2) + X(K,2) * B(2,2)

                    X(1:KE,1) = X(1:KE, 1) -       A(1:KE,K) * WORK(1,1)
                    X(1:KE,2) = X(1:KE, 2) -       A(1:KE,K) * WORK(2,1)

                CASE(3)
                    WORK(1,1) = X(K,1) * B(1,1) + X(K,2) * B(2,1)
                    WORK(2,1) = X(K,1) * B(1,2) + X(K,2) * B(2,2) + X(K,3) * B(3,2)
                    WORK(3,1) = X(K,1) * B(1,3) + X(K,2) * B(2,3) + X(K,3) * B(3,3)

                    X(1:KE,1) = X(1:KE, 1) -       A(1:KE,K) * WORK(1,1)
                    X(1:KE,2) = X(1:KE, 2) -       A(1:KE,K) * WORK(2,1)
                    X(1:KE,3) = X(1:KE, 3) -       A(1:KE,K) * WORK(3,1)

                CASE(4)
                    WORK(1,1) = X(K,1) * B(1,1) + X(K,2) * B(2,1)
                    WORK(2,1) = X(K,1) * B(1,2) + X(K,2) * B(2,2) + X(K,3) * B(3,2)
                    WORK(3,1) = X(K,1) * B(1,3) + X(K,2) * B(2,3) + X(K,3) * B(3,3) + X(K,4) * B(4,3)
                    WORK(4,1) = X(K,1) * B(1,4) + X(K,2) * B(2,4) + X(K,3) * B(3,4) + X(K,4) * B(4,4)

                    X(1:KE,1) = X(1:KE, 1) -       A(1:KE,K) * WORK(1,1)
                    X(1:KE,2) = X(1:KE, 2) -       A(1:KE,K) * WORK(2,1)
                    X(1:KE,3) = X(1:KE, 3) -       A(1:KE,K) * WORK(3,1)
                    X(1:KE,4) = X(1:KE, 4) -       A(1:KE,K) * WORK(4,1)


                END SELECT
            ELSE
                SELECT CASE(N)
                CASE(1)

                    WORK(1,1) = X(K,1)  * B(1,1)
                    WORK(1,2) = X(KH,1) * B(1,1)

                    X(1:KE,1) = X(1:KE, 1) -       (A(1:KE,K) * WORK(1,1) + A(1:KE,KH) * WORK(1,2))

                CASE(2)
                    WORK(1,1) = X(K,1)  * B(1,1) + X(K,2)  * B(2,1)
                    WORK(1,2) = X(KH,1) * B(1,1) + X(KH,2) * B(2,1)
                    WORK(2,1) = X(K,1)  * B(1,2) + X(K,2)  * B(2,2)
                    WORK(2,2) = X(KH,1) * B(1,2) + X(KH,2) * B(2,2)

                    X(1:KE,1) = X(1:KE, 1) -       (A(1:KE,K) * WORK(1,1) + A(1:KE,KH) * WORK(1,2))
                    X(1:KE,2) = X(1:KE, 2) -       (A(1:KE,K) * WORK(2,1) + A(1:KE,KH) * WORK(2,2))

                CASE(3)
                    WORK(1,1) = X(K,1)  * B(1,1) + X(K,2)  * B(2,1)
                    WORK(1,2) = X(KH,1) * B(1,1) + X(KH,2) * B(2,1)
                    WORK(2,1) = X(K,1)  * B(1,2) + X(K,2)  * B(2,2) + X(K,3)  * B(3,2)
                    WORK(2,2) = X(KH,1) * B(1,2) + X(KH,2) * B(2,2) + X(KH,3) * B(3,2)
                    WORK(3,1) = X(K,1)  * B(1,3) + X(K,2)  * B(2,3) + X(K,3)  * B(3,3)
                    WORK(3,2) = X(KH,1) * B(1,3) + X(KH,2) * B(2,3) + X(KH,3) * B(3,3)

                    X(1:KE,1) = X(1:KE, 1) -       (A(1:KE,K) * WORK(1,1) + A(1:KE,KH) * WORK(1,2))
                    X(1:KE,2) = X(1:KE, 2) -       (A(1:KE,K) * WORK(2,1) + A(1:KE,KH) * WORK(2,2))
                    X(1:KE,3) = X(1:KE, 3) -       (A(1:KE,K) * WORK(3,1) + A(1:KE,KH) * WORK(3,2))

                CASE(4)
                    WORK(1,1) = X(K,1)  * B(1,1) + X(K,2)  * B(2,1)
                    WORK(1,2) = X(KH,1) * B(1,1) + X(KH,2) * B(2,1)
                    WORK(2,1) = X(K,1)  * B(1,2) + X(K,2)  * B(2,2) + X(K,3)  * B(3,2)
                    WORK(2,2) = X(KH,1) * B(1,2) + X(KH,2) * B(2,2) + X(KH,3) * B(3,2)
                    WORK(3,1) = X(K,1)  * B(1,3) + X(K,2)  * B(2,3) + X(K,3)  * B(3,3) + X(K,4)  * B(4,3)
                    WORK(3,2) = X(KH,1) * B(1,3) + X(KH,2) * B(2,3) + X(KH,3) * B(3,3) + X(KH,4) * B(4,3)
                    WORK(4,1) = X(K,1)  * B(1,4) + X(K,2)  * B(2,4) + X(K,3)  * B(3,4) + X(K,4)  * B(4,4)
                    WORK(4,2) = X(KH,1) * B(1,4) + X(KH,2) * B(2,4) + X(KH,3) * B(3,4) + X(KH,4) * B(4,4)


                    X(1:KE,1) = X(1:KE, 1) -       (A(1:KE,K) * WORK(1,1) + A(1:KE,KH) * WORK(1,2))
                    X(1:KE,2) = X(1:KE, 2) -       (A(1:KE,K) * WORK(2,1) + A(1:KE,KH) * WORK(2,2))
                    X(1:KE,3) = X(1:KE, 3) -       (A(1:KE,K) * WORK(3,1) + A(1:KE,KH) * WORK(3,2))
                    X(1:KE,4) = X(1:KE, 4) -       (A(1:KE,K) * WORK(4,1) + A(1:KE,KH) * WORK(4,2))


                END SELECT
            END IF
        END IF

        K = K - 1
    END DO

END SUBROUTINE



