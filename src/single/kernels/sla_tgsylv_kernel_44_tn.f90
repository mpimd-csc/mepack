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

!> \brief Solver for a 4x4 generalized Sylvester equation (TRANSA = T, TRANSB = N)
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TGSYLV_KERNEL_44TN (SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, INFO)
!           REAL SGN, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
!           REAL A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TGSYLV_KERNEL_44TN solves a generalized Sylvester equation of the following form
!>
!>    A**T * X * B + SGN * C**T * X * D = SCALE * Y                              (1)
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by SGGES from LAPACK. The algorithm is implemented using BLAS level 2
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
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (LDC,M)
!>          The matrix C must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,M).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (LDD,N)
!>          The matrix D must be (quasi-) upper triangular. If the matrix B is already
!>          quasi-upper triangular the matrix D must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array D.  LDD >= max(1,N).
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
!> \date October 2023
!> \ingroup sgltghelp
!

SUBROUTINE SLA_TGSYLV_KERNEL_44TN (SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
    ! Arguments
    REAL SGN, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
    REAL A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)


    ! Local Variables
    INTEGER K, L, KB, LB, KH, LH, KE, IT
    INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
    REAL MAT(4,4), RHS(4)
    REAL A11, A12, A21, A22, C11, C12, C22
    REAL B11, B12, B21, B22, D11, D12, D21, D22
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

    ! Thresholds for internal solve

    ! (T, N)
    K = 1
    DO WHILE ( K .LE. M )
        KB = 1
        IF ( K .LT. M ) THEN
            IF ( ( A(K+1,K) .NE. ZERO )) THEN
                KB = 2
            END IF
        END IF
        KH = K + KB - 1
        KE = KH + 1

        L = 1
        DO WHILE ( L .LE. N )
            LB = 1
            IF ( L .LT. N ) THEN
                IF ( ( B(L+1,L) .NE. ZERO .OR. D(L+1,L) .NE. ZERO ) ) THEN
                    LB = 2
                END IF
            END IF
            LH = L + LB - 1

            IT = 0
            DO

                ! Solve Inner System
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 1

                    MAT(1,1) = B(L,L) * A(K,K) + SGN * D(L,L) * C(K,K)
                    RHS(1)   = X(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 2

                    A11 = A(K,K)
                    C11 = SGN * C(K,K)

                    D11 = D(L,L)
                    D12 = D(L,LH)
                    D21 = D(LH,L)
                    D22 = D(LH,LH)

                    B11 = B(L,L)
                    B21 = B(LH,L)
                    B12 = B(L,LH)
                    B22 = B(LH,LH)


                    MAT(1,1) = B11*A11 + D11*C11
                    MAT(1,2) = B21*A11 + D21*C11
                    MAT(2,1) = B12*A11 + D12*C11
                    MAT(2,2) = B22*A11 + D22*C11

                    RHS(1) = X(K,L)
                    RHS(2) = X(K,LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2

                    A11 = A(K,K)
                    A12 = A(K,KH)
                    A21 = A(KH,K)
                    A22 = A(KH,KH)
                    C11 = C(K,K)
                    C12 = C(K,KH)
                    C22 = C(KH,KH)

                    D11 = SGN * D(L,L)
                    B11 = B(L,L)

                    MAT(1,1) = B11*A11 + D11*C11
                    MAT(1,2) = B11*A21
                    MAT(2,1) = B11*A12 + D11*C12
                    MAT(2,2) = B11*A22 + D11*C22

                    RHS(1) = X(K,L)
                    RHS(2) = X(KH,L)

                ELSE
                    DIMMAT = 4

                    A11 = A(K,K)
                    A12 = A(K,KH)
                    A21 = A(KH,K)
                    A22 = A(KH,KH)
                    C11 = SGN * C(K,K)
                    C12 = SGN * C(K,KH)
                    C22 = SGN * C(KH,KH)

                    D11 = D(L,L)
                    D12 = D(L,LH)
                    D21 = D(LH,L)
                    D22 = D(LH,LH)

                    B11 = B(L,L)
                    B12 = B(L,LH)
                    B21 = B(LH, L)
                    B22 = B(LH,LH)

                    MAT(1,1) = B11*A11 + D11*C11
                    MAT(1,2) = B11*A21 ! + D11*C21
                    MAT(1,3) = B21*A11 + D21*C11
                    MAT(1,4) = B21*A21 ! + D21*C21

                    MAT(2,1) = B11*A12 + D11*C12
                    MAT(2,2) = B11*A22 + D11*C22
                    MAT(2,3) = B21*A12 + D21*C12
                    MAT(2,4) = B21*A22 + D21*C22

                    MAT(3,1) = B12*A11 + D12*C11
                    MAT(3,2) = B12*A21
                    MAT(3,3) = B22*A11 + D22*C11
                    MAT(3,4) = B22*A21

                    MAT(4,1) = B12*A12 + D12*C12
                    MAT(4,2) = B12*A22 + D12*C22
                    MAT(4,3) = B22*A12 + D22*C12
                    MAT(4,4) = B22*A22 + D22*C22

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

            ! in current row
            IF ( LH .LT. N ) THEN
                IF ( KB .EQ. 1 .AND. LB .EQ. 1) THEN
                    X(K, LH+1:N) = X(K, LH+1:N) - A(K,K) * X(K,L) * B(L,LH+1:N)
                    X(K, LH+1:N) = X(K, LH+1:N) - SGN * C(K,K) * X(K,L) * D(L,LH+1:N)
                ELSE IF (KB .EQ. 1 .AND. LB .EQ. 2) THEN
                    X(K, LH+1:N) = X(K, LH+1:N) - A(K,K) * (X(K,L) * B(L,LH+1:N) + X(K,LH) * B(LH,LH+1:N))
                    X(K, LH+1:N) = X(K, LH+1:N) - SGN * C(K,K) * (X(K,L) * D(L,LH+1:N) + X(K,LH) * D(LH,LH+1:N))

                ELSE IF (KB .EQ. 2 .AND. LB .EQ. 1) THEN
                    X(K, LH+1:N) = X(K, LH+1:N) - (A(K,K) * X(K,L) + A(KH,K) * X(KH,L)) * B(L,LH+1:N)
                    X(K, LH+1:N) = X(K, LH+1:N) - SGN * (C(K,K) * X(K,L) + C(KH,K) * X(KH,L)) * D(L,LH+1:N)

                    X(KH, LH+1:N) = X(KH, LH+1:N) - (A(K,KH) * X(K,L) + A(KH,KH) * X(KH,L)) * B(L,LH+1:N)
                    X(KH, LH+1:N) = X(KH, LH+1:N) - SGN * (C(K,KH) * X(K,L) + C(KH,KH) * X(KH,L)) * D(L,LH+1:N)
                ELSE
                    WORK(1,1) = A(K,K) * X(K,L) + A(KH,K) * X(KH,L)
                    WORK(1,2) = A(K,K) * X(K,LH) + A(KH,K) * X(KH,LH)
                    WORK(2,1) = A(K,KH) * X(K,L) + A(KH,KH) * X(KH,L)
                    WORK(2,2) = A(K,KH) * X(K,LH) + A(KH,KH) * X(KH,LH)

                    X(K, LH+1:N) = X(K, LH+1:N) - (WORK(1,1) *  B(L,LH+1:N) + WORK(1,2) * B(LH,LH+1:N))
                    X(KH, LH+1:N) = X(KH, LH+1:N) - (WORK(2,1) *  B(L,LH+1:N) + WORK(2,2) * B(LH,LH+1:N))

                    WORK(1,1) = C(K,K) * X(K,L) + C(KH,K) * X(KH,L)
                    WORK(1,2) = C(K,K) * X(K,LH) + C(KH,K) * X(KH,LH)
                    WORK(2,1) = C(K,KH) * X(K,L) + C(KH,KH) * X(KH,L)
                    WORK(2,2) = C(K,KH) * X(K,LH) + C(KH,KH) * X(KH,LH)

                    X(K, LH+1:N) = X(K, LH+1:N) - SGN * (WORK(1,1) *  D(L,LH+1:N) + WORK(1,2) * D(LH,LH+1:N))
                    X(KH, LH+1:N) = X(KH, LH+1:N) - SGN * (WORK(2,1) *  D(L,LH+1:N) + WORK(2,2) * D(LH,LH+1:N))
                ENDIF
            END IF


            L = LH + 1
        END DO

        IF ( KH .LT. M ) THEN
            IF ( KB .EQ. 1 ) THEN
                SELECT CASE(N)
                CASE(1)
                    WORK(1,1) = X(K,1) * B(1,1)
                    WORK(1,2) = SGN* X(K,1) * D(1,1)

                    X(KE:M, 1) = X(KE:M, 1) - A(K,KE:M) * WORK(1,1) &
                        &   - C(K,KE:M) * WORK(1,2)
                CASE(2)
                    WORK(1,1) = X(K,1) * B(1,1) + X(K,2) * B(2,1)
                    WORK(2,1) = X(K,1) * B(1,2) + X(K,2) * B(2,2)

                    WORK(1,2) = X(K,1) * D(1,1) + X(K,2) * D(2,1)
                    WORK(2,2) = X(K,1) * D(1,2) + X(K,2) * D(2,2)

                    X(KE:M, 1) = X(KE:M, 1) - A(K,KE:M) * WORK(1,1) &
                        &   - SGN * C(K,KE:M) * WORK(1,2)
                    X(KE:M, 2) = X(KE:M, 2) - A(K,KE:M) * WORK(2,1) &
                        &   - SGN * C(K,KE:M) * WORK(2,2)


                CASE(3)
                    WORK(1,1) = X(K,1) * B(1,1) + X(K,2) * B(2,1)
                    WORK(2,1) = X(K,1) * B(1,2) + X(K,2) * B(2,2) + X(K,3) * B(3,2)
                    WORK(3,1) = X(K,1) * B(1,3) + X(K,2) * B(2,3) + X(K,3) * B(3,3)

                    WORK(1,2) = X(K,1) * D(1,1) + X(K,2) * D(2,1)
                    WORK(2,2) = X(K,1) * D(1,2) + X(K,2) * D(2,2) + X(K,3) * D(3,2)
                    WORK(3,2) = X(K,1) * D(1,3) + X(K,2) * D(2,3) + X(K,3) * D(3,3)

                    X(KE:M, 1) = X(KE:M, 1) - A(K,KE:M) * WORK(1,1) &
                        &   - SGN * C(K,KE:M) * WORK(1,2)
                    X(KE:M, 2) = X(KE:M, 2) - A(K,KE:M) * WORK(2,1) &
                        &   - SGN * C(K,KE:M) * WORK(2,2)
                    X(KE:M, 3) = X(KE:M, 3) - A(K,KE:M) * WORK(3,1) &
                        &   - SGN * C(K,KE:M) * WORK(3,2)

                CASE(4)
                    WORK(1,1) = X(K,1) * B(1,1) + X(K,2) * B(2,1)
                    WORK(2,1) = X(K,1) * B(1,2) + X(K,2) * B(2,2) + X(K,3) * B(3,2)
                    WORK(3,1) = X(K,1) * B(1,3) + X(K,2) * B(2,3) + X(K,3) * B(3,3) + X(K,4) * B(4,3)
                    WORK(4,1) = X(K,1) * B(1,4) + X(K,2) * B(2,4) + X(K,3) * B(3,4) + X(K,4) * B(4,4)

                    WORK(1,2) = X(K,1) * D(1,1) + X(K,2) * D(2,1)
                    WORK(2,2) = X(K,1) * D(1,2) + X(K,2) * D(2,2) + X(K,3) * D(3,2)
                    WORK(3,2) = X(K,1) * D(1,3) + X(K,2) * D(2,3) + X(K,3) * D(3,3) + X(K,4) * D(4,3)
                    WORK(4,2) = X(K,1) * D(1,4) + X(K,2) * D(2,4) + X(K,3) * D(3,4) + X(K,4) * D(4,4)

                    X(KE:M, 1) = X(KE:M, 1) - A(K,KE:M) * WORK(1,1) &
                        &   - SGN * C(K,KE:M) * WORK(1,2)
                    X(KE:M, 2) = X(KE:M, 2) - A(K,KE:M) * WORK(2,1) &
                        &   - SGN * C(K,KE:M) * WORK(2,2)
                    X(KE:M, 3) = X(KE:M, 3) - A(K,KE:M) * WORK(3,1) &
                        &   - SGN * C(K,KE:M) * WORK(3,2)
                    X(KE:M, 4) = X(KE:M, 4) - A(K,KE:M) * WORK(4,1) &
                        &   - SGN * C(K,KE:M) * WORK(4,2)
                END SELECT
            ELSE
                SELECT CASE(N)
                CASE(1)
                    WORK(1,1) = X(K,1)  * B(1,1)
                    WORK(1,2) = X(KH,1) * B(1,1)

                    WORK(1,3) = X(K,1) * D(1,1)
                    WORK(1,4) = X(KH,1) * D(1,1)

                    X(KE:M, 1) = X(KE:M, 1) - (A(K,KE:M) * WORK(1,1) + A(KH, KE:M) * WORK(1,2)) &
                        &   - SGN * (C(K,KE:M) * WORK(1,3)  + C(KH, KE:M) * WORK(1,4))



                CASE(2)
                    WORK(1,1) = X(K,1)  * B(1,1) + X(K,2)  * B(2,1)
                    WORK(1,2) = X(KH,1) * B(1,1) + X(KH,2) * B(2,1)
                    WORK(2,1) = X(K,1)  * B(1,2) + X(K,2)  * B(2,2)
                    WORK(2,2) = X(KH,1) * B(1,2) + X(KH,2) * B(2,2)

                    WORK(1,3) = X(K,1)  * D(1,1) + X(K,2)  * D(2,1)
                    WORK(1,4) = X(KH,1) * D(1,1) + X(KH,2) * D(2,1)
                    WORK(2,3) = X(K,1)  * D(1,2) + X(K,2)  * D(2,2)
                    WORK(2,4) = X(KH,1) * D(1,2) + X(KH,2) * D(2,2)

                    X(KE:M, 1) = X(KE:M, 1) - (A(K,KE:M) * WORK(1,1) + A(KH, KE:M) * WORK(1,2)) &
                        &   - SGN * (C(K,KE:M) * WORK(1,3)  + C(KH, KE:M) * WORK(1,4))
                    X(KE:M, 2) = X(KE:M, 2) - (A(K,KE:M) * WORK(2,1) + A(KH, KE:M) * WORK(2,2)) &
                        &   - SGN * (C(K,KE:M) * WORK(2,3)  + C(KH, KE:M) * WORK(2,4))

                CASE(3)
                    WORK(1,1) = X(K,1)  * B(1,1) + X(K,2)  * B(2,1)
                    WORK(1,2) = X(KH,1) * B(1,1) + X(KH,2) * B(2,1)
                    WORK(2,1) = X(K,1)  * B(1,2) + X(K,2)  * B(2,2) + X(K,3)  * B(3,2)
                    WORK(2,2) = X(KH,1) * B(1,2) + X(KH,2) * B(2,2) + X(KH,3) * B(3,2)
                    WORK(3,1) = X(K,1)  * B(1,3) + X(K,2)  * B(2,3) + X(K,3)  * B(3,3)
                    WORK(3,2) = X(KH,1) * B(1,3) + X(KH,2) * B(2,3) + X(KH,3) * B(3,3)

                    WORK(1,3) = X(K,1)  * D(1,1) + X(K,2)  * D(2,1)
                    WORK(1,4) = X(KH,1) * D(1,1) + X(KH,2) * D(2,1)
                    WORK(2,3) = X(K,1)  * D(1,2) + X(K,2)  * D(2,2) + X(K,3)  * D(3,2)
                    WORK(2,4) = X(KH,1) * D(1,2) + X(KH,2) * D(2,2) + X(KH,3) * D(3,2)
                    WORK(3,3) = X(K,1)  * D(1,3) + X(K,2)  * D(2,3) + X(K,3)  * D(3,3)
                    WORK(3,4) = X(KH,1) * D(1,3) + X(KH,2) * D(2,3) + X(KH,3) * D(3,3)

                    X(KE:M, 1) = X(KE:M, 1) - (A(K,KE:M) * WORK(1,1) + A(KH, KE:M) * WORK(1,2)) &
                        &   - SGN * (C(K,KE:M) * WORK(1,3)  + C(KH, KE:M) * WORK(1,4))
                    X(KE:M, 2) = X(KE:M, 2) - (A(K,KE:M) * WORK(2,1) + A(KH, KE:M) * WORK(2,2)) &
                        &   - SGN * (C(K,KE:M) * WORK(2,3)  + C(KH, KE:M) * WORK(2,4))
                    X(KE:M, 3) = X(KE:M, 3) - (A(K,KE:M) * WORK(3,1) + A(KH, KE:M) * WORK(3,2)) &
                        &   - SGN * (C(K,KE:M) * WORK(3,3)  + C(KH, KE:M) * WORK(3,4))

                CASE(4)
                    WORK(1,1) = X(K,1)  * B(1,1) + X(K,2)  * B(2,1)
                    WORK(1,2) = X(KH,1) * B(1,1) + X(KH,2) * B(2,1)
                    WORK(2,1) = X(K,1)  * B(1,2) + X(K,2)  * B(2,2) + X(K,3)  * B(3,2)
                    WORK(2,2) = X(KH,1) * B(1,2) + X(KH,2) * B(2,2) + X(KH,3) * B(3,2)
                    WORK(3,1) = X(K,1)  * B(1,3) + X(K,2)  * B(2,3) + X(K,3)  * B(3,3) + X(K,4)  * B(4,3)
                    WORK(3,2) = X(KH,1) * B(1,3) + X(KH,2) * B(2,3) + X(KH,3) * B(3,3) + X(KH,4) * B(4,3)
                    WORK(4,1) = X(K,1)  * B(1,4) + X(K,2)  * B(2,4) + X(K,3)  * B(3,4) + X(K,4)  * B(4,4)
                    WORK(4,2) = X(KH,1) * B(1,4) + X(KH,2) * B(2,4) + X(KH,3) * B(3,4) + X(KH,4) * B(4,4)

                    WORK(1,3) = X(K,1)  * D(1,1) + X(K,2)  * D(2,1)
                    WORK(1,4) = X(KH,1) * D(1,1) + X(KH,2) * D(2,1)
                    WORK(2,3) = X(K,1)  * D(1,2) + X(K,2)  * D(2,2) + X(K,3)  * D(3,2)
                    WORK(2,4) = X(KH,1) * D(1,2) + X(KH,2) * D(2,2) + X(KH,3) * D(3,2)
                    WORK(3,3) = X(K,1)  * D(1,3) + X(K,2)  * D(2,3) + X(K,3)  * D(3,3) + X(K,4)  * D(4,3)
                    WORK(3,4) = X(KH,1) * D(1,3) + X(KH,2) * D(2,3) + X(KH,3) * D(3,3) + X(KH,4) * D(4,3)
                    WORK(4,3) = X(K,1)  * D(1,4) + X(K,2)  * D(2,4) + X(K,3)  * D(3,4) + X(K,4)  * D(4,4)
                    WORK(4,4) = X(KH,1) * D(1,4) + X(KH,2) * D(2,4) + X(KH,3) * D(3,4) + X(KH,4) * D(4,4)

                    X(KE:M, 1) = X(KE:M, 1) - (A(K,KE:M) * WORK(1,1) + A(KH, KE:M) * WORK(1,2)) &
                        &   - SGN * (C(K,KE:M) * WORK(1,3)  + C(KH, KE:M) * WORK(1,4))
                    X(KE:M, 2) = X(KE:M, 2) - (A(K,KE:M) * WORK(2,1) + A(KH, KE:M) * WORK(2,2)) &
                        &   - SGN * (C(K,KE:M) * WORK(2,3)  + C(KH, KE:M) * WORK(2,4))
                    X(KE:M, 3) = X(KE:M, 3) - (A(K,KE:M) * WORK(3,1) + A(KH, KE:M) * WORK(3,2)) &
                        &   - SGN * (C(K,KE:M) * WORK(3,3)  + C(KH, KE:M) * WORK(3,4))
                    X(KE:M, 4) = X(KE:M, 4) - (A(K,KE:M) * WORK(4,1) + A(KH, KE:M) * WORK(4,2)) &
                        &   - SGN * (C(K,KE:M) * WORK(4,3)  + C(KH, KE:M) * WORK(4,4))
                END SELECT

            END IF
        END IF

        K = KH + 1
    END DO

END SUBROUTINE



