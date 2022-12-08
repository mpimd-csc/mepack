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


!> \brief Level-2 Bartels-Stewart Algorithm for the generalized Sylvester equation (unoptimized)
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TGSYLV_L2_UNOPT ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANSA(1), TRANSB(1)
!           REAL SGN, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
!           REAL A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
!           REAL WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TGSYLV_L2_UNOPT solves a generalized Sylvester equation of the following forms
!>
!>    op1(A) * X * op2(B) + op1(C) * X * op2(D) = SCALE * Y                              (1)
!>
!> or
!>
!>    op1(A) * X * op2(B) - op1(C) * X * op2(D) = SCALE * Y                              (2)
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by SGGES from LAPACK.
!> \endverbatim
!> \attention The algorithm is implemented using BLAS level 2 operations without further optimizations. For a faster implementation
!> see SLA_TGSYLV_L2 .
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A and C:
!>          == 'N':  op1(A) = A, op1(C) = C (No transpose for A and C)
!>          == 'T':  op1(A) = A**T, op1(C) = C **T (Transpose A and C)
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to B and D:
!>          == 'N':  op2(B) = B, op2(D) = D (No transpose for B and D)
!>          == 'T':  op2(B) = B**T, op2(D) = D **T (Transpose B and D)
!> \endverbatim
!>
!> \param[in] SGN
!> \verbatim
!>          SGN is REAL, allowed values: +/-1
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
!>          could not be solved correctly, 0 < SCALE <= 1 holds true.
!> \endverbatim
!>
!> \param[in] WORK
!> \verbatim
!>          WORK is REAL array, dimension LWORK
!>          Workspace for the algorithm.
!>          The workspace needs to queried before the running the computation.
!>          The query is performed by calling the subroutine with INFO == -1 on input.
!>          The required workspace is then returned in INFO.
!> \endverbatim
!>
!> \param[inout] INFO
!> \verbatim
!>          INFO is INTEGER
!>
!>          On input:
!>            == -1 : Perform a workspace query
!>            <> -1: normal operation
!>
!>          On exit, workspace query:
!>            < 0 :  if INFO = -i, the i-th argument had an illegal value
!>            >= 0:  The value of INFO is the required number of elements in the workspace.
!>
!>          On exit, normal operation:
!>            == 0:  successful exit
!>            < 0:  if INFO = -i, the i-th argument had an illegal value
!>            > 0:  The equation is not solved correctly. One of the arising inner
!>                  system got singular.
!> \endverbatim
!>
!>
!> \par Optimizations:
!       ==============
!>
!> \li Nothing.
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Dezember 2022
!> \ingroup sgltgsylv
!
SUBROUTINE SLA_TGSYLV_L2_ROWWISE ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    REAL SGN, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
    REAL A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
    REAL WORK(*)


    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE, IT
    INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
    REAL MAT(4,4), RHS(4)
    REAL A11, A12, A21, A22, C11, C12, C22
    REAL B11, B12, B21, B22, D11, D12, D21, D22
    REAL SCAL

    LOGICAL BTRANSA, BTRANSB
    REAL ZERO , ONE
    PARAMETER (ZERO = 0.0, ONE = 1.0)

    ! EXTERNAL FUNCTIONS
    EXTERNAL SGEMM
    EXTERNAL SGESC2
    EXTERNAL SGETC2X
    EXTERNAL SLACPY
    EXTERNAL SLA_QTRMM2
    EXTERNAL SLA_SMALL_SOLVE4
    EXTERNAL STRMM
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
    ELSE IF ( LDC .LT. MAX(1, M)) THEN
        INFO = -11
    ELSE IF ( LDD .LT. MAX(1, N)) THEN
        INFO = -13
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -15
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_TGSYLV_L2_ROWWISE', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1) THEN
        INFO = 2*MAX(M,N)
        RETURN
    END IF

    ! Quick Return
    SCALE = ONE
    IF ( M .EQ. 0 .OR. N .EQ. 0 ) THEN
        RETURN
    END IF

    ! Prepare
    INFO = 0

    IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB ) THEN
        ! ( N, N )
        K = M
        DO WHILE (K .GT. 0)
            KB = 1
            KH = K
            IF ( K .GT. 1) THEN
                IF (A(K,K-1) .NE. ZERO ) THEN
                    KB = 2
                    K = K - 1
                END IF
            END IF
            KE = K - 1

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
                    ! Solve inner system
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
                        MAT(1,2) = B11*A12 + D11*C12
                        MAT(2,1) = B11*A21
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
                        MAT(1,2) = B11*A12 + D11*C12
                        MAT(1,3) = B21*A11 + D21*C11
                        MAT(1,4) = B21*A12 + D21*C12

                        MAT(2,1) = B11*A21
                        MAT(2,2) = B11*A22 + D11*C22
                        MAT(2,3) = B21*A21
                        MAT(2,4) = B21*A22 + D21*C22

                        MAT(3,1) = B12*A11 + D12*C11
                        MAT(3,2) = B12*A12 + D12*C12
                        MAT(3,3) = B22*A11 + D22*C11
                        MAT(3,4) = B22*A12 + D22*C12

                        MAT(4,1) = B12*A21
                        MAT(4,2) = B12*A22 + D12*C22
                        MAT(4,3) = B22*A21
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

                ! Update January 2021
                ! in current row
                IF ( LH .LT. N ) THEN
                    CALL SGEMM('N','N', KB, LB, KB, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, MAT, KB )
                    CALL SGEMM('N','N', KB, N-LH, LB, -ONE, MAT, KB, B(L,LH+1), LDB, ONE, X(K,LH+1), LDX)

                    CALL SGEMM('N','N', KB, LB, KB, ONE, C(K,K), LDC, X(K,L), LDX, ZERO, MAT, KB )
                    CALL SGEMM('N','N', KB, N-LH, LB, -ONE*SGN, MAT, KB, D(L,LH+1), LDD, ONE, X(K,LH+1), LDX)

                END IF

                L = LH + 1
            END DO

            ! Update January 2021
            IF ( K .GT. 1 ) THEN
                CALL SLACPY('All', KB, N, X(K,1), LDX,  WORK, 2)
                CALL STRMM('Right', 'Upper', 'NoTrans', 'NoUnit', KB, N, ONE, B, LDB, WORK, 2)
                CALL SLA_QTRMM2('Right','NoTrans',KB, N, ONE, B, LDB, X(K,1), LDX, WORK, 2)
                CALL SGEMM('N', 'N', KE, N, KB, -ONE, A(1,K), LDA, WORK, 2, ONE, X(1,1), LDX)

                CALL SLACPY('All', KB, N, X(K,1), LDX, WORK, 2)
                CALL STRMM('Right', 'Upper', 'NoTrans', 'NoUnit', KB, N, ONE, D, LDD, WORK, 2)
                CALL SLA_QTRMM2('Right','NoTrans',KB, N, ONE, D, LDD, X(K,1), LDX, WORK, 2)
                CALL SGEMM('N', 'N', KE, N, KB, -ONE*SGN, C(1,K), LDC, WORK, 2, ONE, X(1,1), LDX)
            END IF

            K = K - 1
        END DO

    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N,T)
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

            L = N
            DO WHILE ( L .GT. 0 )
                LB = 1
                LH = L
                IF ( L .GT. 1 ) THEN
                    IF ( (B(L,L-1) .NE. ZERO .OR. D(L,L-1) .NE. ZERO) ) THEN
                        LB = 2
                        L = L - 1
                    END IF
                END IF
                LE = L - 1

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
                        MAT(1,2) = B12*A11 + D12*C11
                        MAT(2,1) = B21*A11 + D21*C11
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
                        MAT(1,2) = B11*A12 + D11*C12
                        MAT(2,1) = B11*A21
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
                        MAT(1,2) = B11*A12 + D11*C12
                        MAT(1,3) = B12*A11 + D12*C11
                        MAT(1,4) = B12*A12 + D12*C12

                        MAT(2,1) = B11*A21
                        MAT(2,2) = B11*A22 + D11*C22
                        MAT(2,3) = B12*A21
                        MAT(2,4) = B12*A22 + D12*C22

                        MAT(3,1) = B21*A11 + D21*C11
                        MAT(3,2) = B21*A12 + D21*C12
                        MAT(3,3) = B22*A11 + D22*C11
                        MAT(3,4) = B22*A12 + D22*C12

                        MAT(4,1) = B21*A21
                        MAT(4,2) = B21*A22 + D21*C22
                        MAT(4,3) = B22*A21
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

                ! Update January 2021
                IF ( L .GT. 0 ) THEN
                    CALL SGEMM('N','N', KB, LB, KB, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, MAT, KB )
                    CALL SGEMM('N','T', KB, LE, LB, -ONE, MAT, KB, B(1,L), LDB, ONE, X(K,1), LDX)

                    CALL SGEMM('N','N', KB, LB, KB, ONE, C(K,K), LDC, X(K,L), LDX, ZERO, MAT, KB )
                    CALL SGEMM('N','T', KB, LE, LB, -ONE*SGN, MAT, KB, D(1,L), LDD, ONE, X(K,1), LDX)
                END IF

                L = L - 1
            END DO

            ! Update January 2021
            IF ( K .GT. 1 ) THEN
                CALL SLACPY('All', KB, N, X(K,1), LDX,  WORK, 2)
                CALL STRMM('Right', 'Upper', 'Trans', 'NoUnit', KB, N, ONE, B, LDB, WORK, 2)
                CALL SLA_QTRMM2('Right','Trans',KB, N, ONE, B, LDB, X(K,1), LDX, WORK, 2)
                CALL SGEMM('N', 'N', KE, N, KB, -ONE, A(1,K), LDA, WORK, 2, ONE, X(1,1), LDX)

                CALL SLACPY('All', KB, N, X(K,1), LDX, WORK, 2)
                CALL STRMM('Right', 'Upper', 'Trans', 'NoUnit', KB, N, ONE, D, LDD, WORK, 2)
                CALL SLA_QTRMM2('Right','Trans',KB, N, ONE, D, LDD, X(K,1), LDX, WORK, 2)
                CALL SGEMM('N', 'N', KE, N, KB, -ONE*SGN, C(1,K), LDC, WORK, 2, ONE, X(1,1), LDX)
            END IF

            K = K - 1
        END DO

    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
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
                    CALL SGEMM('T','N', KB, LB, KB, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, MAT, KB )
                    CALL SGEMM('N','N', KB, N-LH, LB, -ONE, MAT, KB, B(L,LH+1), LDB, ONE, X(K,LH+1), LDX)

                    CALL SGEMM('T','N', KB, LB, KB, ONE, C(K,K), LDC, X(K,L), LDX, ZERO, MAT, KB )
                    CALL SGEMM('N','N', KB, N-LH, LB, -ONE*SGN, MAT, KB, D(L,LH+1), LDD, ONE, X(K,LH+1), LDX)

                END IF


                L = LH + 1
            END DO

            IF ( KH .LT. M ) THEN
                CALL SLACPY('All', KB, N, X(K,1), LDX,  WORK, 2)
                CALL STRMM('Right', 'Upper', 'NoTrans', 'NoUnit', KB, N, ONE, B, LDB, WORK, 2)
                CALL SLA_QTRMM2('Right','NoTrans',KB, N, ONE, B, LDB, X(K,1), LDX, WORK, 2)
                CALL SGEMM('T', 'N', M-KH, N, KB, -ONE, A(K,KH+1), LDA, WORK, 2, ONE, X(KH+1,1), LDX)

                CALL SLACPY('All', KB, N, X(K,1), LDX,  WORK, 2)
                CALL STRMM('Right', 'Upper', 'NoTrans', 'NoUnit', KB, N, ONE, D, LDD, WORK, 2)
                CALL SLA_QTRMM2('Right','NoTrans',KB, N, ONE, D, LDD, X(K,1), LDX, WORK, 2)
                CALL SGEMM('T', 'N', M-KH, N, KB, -ONE*SGN, C(K,KH+1), LDC, WORK, 2, ONE, X(KH+1,1), LDX)

            END IF

            K = KH + 1
        END DO


    ELSE IF ( BTRANSA .AND. BTRANSB) THEN
        ! (T, T)
        K = 1
        DO WHILE ( K .LE. M )
            KB = 1
            IF ( K .LT. M ) THEN
                IF ( ( A(K+1,K) .NE. ZERO )) THEN
                    KB = 2
                END IF
            END IF
            KH = K + KB - 1

            L = N
            DO WHILE ( L .GT. 0 )
                LB = 1
                LH = L
                IF ( L .GT. 1 ) THEN
                    IF (B(L,L-1) .NE. ZERO .OR. D(L,L-1) .NE. ZERO)  THEN
                        LB = 2
                        L = L - 1
                    END IF
                END IF
                LE = L - 1

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
                        MAT(1,2) = B12*A11 + D12*C11
                        MAT(2,1) = B21*A11 + D21*C11
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
                        MAT(1,3) = B12*A11 + D12*C11
                        MAT(1,4) = B12*A21 ! + D12*C21

                        MAT(2,1) = B11*A12 + D11*C12
                        MAT(2,2) = B11*A22 + D11*C22
                        MAT(2,3) = B12*A12 + D12*C12
                        MAT(2,4) = B12*A22 + D12*C22

                        MAT(3,1) = B21*A11 + D21*C11
                        MAT(3,2) = B21*A21
                        MAT(3,3) = B22*A11 + D22*C11
                        MAT(3,4) = B22*A21

                        MAT(4,1) = B21*A12 + D21*C12
                        MAT(4,2) = B21*A22 + D21*C22
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
                IF ( L .GT. 1  ) THEN
                    CALL SGEMM('T','N', KB, LB, KB, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, MAT, KB )
                    CALL SGEMM('N','T', KB, LE, LB, -ONE, MAT, KB, B(1,L), LDB, ONE, X(K,1), LDX)

                    CALL SGEMM('T','N', KB, LB, KB, ONE, C(K,K), LDC, X(K,L), LDX, ZERO, MAT, KB )
                    CALL SGEMM('N','T', KB, LE, LB, -ONE*SGN, MAT, KB, D(1,L), LDD, ONE, X(K,1), LDX)

                END IF

                L = L - 1
            END DO

            ! Update January 2021
            IF ( KH .LT. M ) THEN
                CALL SLACPY('All', KB, N, X(K,1), LDX,  WORK, 2)
                CALL STRMM('Right', 'Upper', 'Trans', 'NoUnit', KB, N, ONE, B, LDB, WORK, 2)
                CALL SLA_QTRMM2('Right','Trans',KB, N, ONE, B, LDB, X(K,1), LDX, WORK, 2)
                CALL SGEMM('T', 'N', M-KH, N, KB, -ONE, A(K,KH+1), LDA, WORK, 2, ONE, X(KH+1,1), LDX)

                CALL SLACPY('All', KB, N, X(K,1), LDX,  WORK, 2)
                CALL STRMM('Right', 'Upper', 'Trans', 'NoUnit', KB, N, ONE, D, LDD, WORK, 2)
                CALL SLA_QTRMM2('Right','Trans',KB, N, ONE, D, LDD, X(K,1), LDX, WORK, 2)
                CALL SGEMM('T', 'N', M-KH, N, KB, -ONE*SGN, C(K,KH+1), LDC, WORK, 2, ONE, X(KH+1,1), LDX)

            END IF

            K = KH + 1
        END DO

    END IF
END SUBROUTINE



SUBROUTINE SLA_TGSYLV_L2_UNOPT ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    REAL SGN, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
    REAL A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
    REAL WORK(*)

    LOGICAL BTRANSA, BTRANSB
    REAL ONE

    EXTERNAL SLA_TGSYLV_L2_ROWWISE
    EXTERNAL LSAME
    EXTERNAL XERROR_HANDLER
    LOGICAL LSAME
    PARAMETER ( ONE = 1.0)

    IF (INFO.EQ.-1) THEN
        CALL SLA_TGSYLV_L2_ROWWISE ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, INFO)
        RETURN
    END IF
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
    ELSE IF ( LDC .LT. MAX(1, M)) THEN
        INFO = -11
    ELSE IF ( LDD .LT. MAX(1, N)) THEN
        INFO = -13
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -15
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_TGSYLV_L2_UNOPT', -INFO)
        RETURN
    END IF

    CALL SLA_TGSYLV_L2_ROWWISE ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, INFO)

END SUBROUTINE
