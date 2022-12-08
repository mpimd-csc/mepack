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

!> \brief Level-2 Bartels-Stewart Algorithm for the generalized Stein equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TGSTEIN_L2 ( TRANS, M,  A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANS(1)
!           DOUBLE PRECISION SCALE
!           INTEGER M, N, LDA, LDB, LDX, INFO
!           DOUBLE PRECISION A(LDA, *), B(LDB, *) , X(LDX, *)
!           DOUBLE PRECISION WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLA_TGSTEIN_L2 solves a generalized Lyapunov  equation of the following forms
!>
!>    A * X * A^T - B * X * B^T = SCALE * Y                                              (1)
!>
!> or
!>
!>    A^T * X * A - B^T * X * B =  SCALE * Y                                             (2)
!>
!> where A is a M-by-M quasi upper triangular matrix, B is a M-by-M upper triangular,
!> and X and Y are symmetric  M-by-M matrices.
!> Typically the matrix pencil (A,B) is created by DGGES from LAPACK.
!> \endverbatim
!> \attention The algorithm is implemented using BLAS level 2 operations.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A, B, and C:
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
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,M)
!>          The matrix B must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
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
!>          WORK is DOUBLE PRECISION array, dimension LWORK
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
!
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Dezember 2022
!> \ingroup dbltglyap
!
SUBROUTINE DLA_TGSTEIN_L2 ( TRANS, M, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANS(1)
    DOUBLE PRECISION SCALE
    INTEGER M, LDA, LDB, LDX, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *), X(LDX, *)
    DOUBLE PRECISION WORK(*)

    ! Local Variables
    INTEGER ININFO

    INTEGER K, L, KB, LB, KH, LE, LH, KE, IT, II, JJ
    INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
    DOUBLE PRECISION MAT(4,4), RHS(4)
    DOUBLE PRECISION A11, A12, A21, A22, C11, C12, C22
    DOUBLE PRECISION B11, B12, B21, B22, D11, D12, D21, D22
    DOUBLE PRECISION SCAL


    LOGICAL BTRANSA
    DOUBLE PRECISION ZERO , ONE, HALF
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0)

    ! DOUBLE PRECISION WORKA(256), WORKB(256), WORKC(256), WORKD(256)
    !    !dir$ attributes align: 64:: WORKA, WORKB, WORKC, WORKD


    INTEGER IONE, I
    PARAMETER ( IONE = 1)

    ! EXTERNAL FUNCTIONS
    EXTERNAL DGEMV
    EXTERNAL DGER
    EXTERNAL DGESC2
    EXTERNAL DGETC2X
    EXTERNAL DLA_SMALL_SOLVE4
    EXTERNAL DSYR2
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME


    ! Check Input
    ININFO = INFO
    INFO = 0
    BTRANSA = LSAME(TRANS, 'T')

    IF ( .NOT. BTRANSA .AND. .NOT. LSAME(TRANS, 'N')) THEN
        INFO = -1
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -2
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -4
    ELSE IF ( LDB .LT. MAX(1, M)) THEN
        INFO = -6
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -8
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('DLA_TGSTEIN_L2', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1) THEN
        INFO = 4*M
        RETURN
    END IF

    ! Quick Return
    SCALE = ONE
    IF ( M .EQ. 0 ) THEN
        RETURN
    END IF

    ! Prepare
    INFO = 0

    IF ( .NOT. BTRANSA ) THEN
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
                    IF (  A(K,K-1) .NE. ZERO ) THEN
                        KB = 2
                        K = K - 1
                    END IF
                END IF
                KE = K - 1

                IF (K .EQ. L .AND. KB .EQ. 2 ) THEN
                    X(KH, K ) = X(K,KH)
                END IF


                IT = 0
                DO
                    ! Solve Inner System
                    IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 1

                        MAT(1,1) = A(L,L) * A(K,K) - B(L,L) * B(K,K)
                        RHS(1)   = X(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = A(K,K)
                        C11 = - B(K,K)

                        D11 = B(L,L)
                        D12 = B(L,LH)
                        D21 = B(LH,L)
                        D22 = B(LH,LH)

                        B11 = A(L,L)
                        B21 = A(LH,L)
                        B12 = A(L,LH)
                        B22 = A(LH,LH)


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
                        C11 = B(K,K)
                        C12 = B(K,KH)
                        C22 = B(KH,KH)

                        D11 = - B(L,L)
                        B11 = A(L,L)

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
                        C11 = - B(K,K)
                        C12 = - B(K,KH)
                        C22 =  -B(KH,KH)

                        D11 = B(L,L)
                        D12 = B(L,LH)
                        D21 = B(LH,L)
                        D22 = B(LH,LH)

                        B11 = A(L,L)
                        B12 = A(L,LH)
                        B21 = A(LH, L)
                        B22 = A(LH,LH)

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
                        CALL DLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                        IF ( INFO1 .EQ. 0 ) EXIT
                        IT = 1
                    ELSE
                        CALL DGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                        IF ( INFO1 .NE. 0 ) INFO = INFO1
                        CALL DGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
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

                IF ( K.EQ.L .AND. KB.EQ.2) THEN
                    X(K,KH) = HALF*(RHS(2)+RHS(3))
                    X(KH,K) = X(K,KH)
                ELSE
                    ! Transpose causes a temporary here
                    ! X(L:LH,K:KH) = TRANSPOSE(X(K:KH,L:LH))
                      DO II = K, KH
                            DO JJ = L, LH
                              X(JJ, II) = X(II,JJ)
                          END DO
                    END DO
                END IF

                IF ( K .GT. 1 ) THEN
                    IF ( KB .EQ. IONE .AND. LB .EQ. IONE ) THEN
                        X(1:KE, L ) = X(1:KE, L) - (A(L,L) * X(K,L))* A(1:KE, K) + (B(L,L)*X(K,L)) * B(1:KE,K)
                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                        X(1:KE,L) = X(1:KE,L) - (X(K,L) * A(L,L)  + X(K,LH) * A(L,LH))*A(1:KE,K) &
                            & + ((X(K,L) * B(L,L) + X(K,LH) * B(L,LH)))* B(1:KE,K)
                        X(1:KE,LH) = X(1:KE,LH) - (X(K,L) * A(LH,L)  + X(K,LH) * A(LH,LH))*A(1:KE,K) &
                            & + ((X(K,L) * B(LH,L) + X(K,LH) * B(LH,LH)))* B(1:KE,K)
                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                        X(1:KE, L) = X(1:KE,L) - A(1:KE,K) *(X(K,L) * A(L,L)) - A(1:KE,KH) *(X(KH,L) * A(L,L)) &
                            &  + B(1:KE,K) *(X(K,L) * B(L,L) ) + B(1:KE,KH) *(X(KH,L) * B(L,L) )
                    ELSE
                        MAT(1,1) = X(K,L) *A(L,L)  + X(K,LH)  * A(L,LH)
                        MAT(2,1) = X(KH,L)*A(L,L)  + X(KH,LH) * A(L,LH)
                        MAT(3,1) = X(K,L) *A(LH,L) + X(K,LH)  * A(LH,LH)
                        MAT(4,1) = X(KH,L)*A(LH,L) + X(KH,LH) * A(LH,LH)

                        MAT(1,2) = X(K,L) *B(L,L)  + X(K,LH)  * B(L,LH)
                        MAT(2,2) = X(KH,L)*B(L,L)  + X(KH,LH) * B(L,LH)
                        MAT(3,2) = X(K,L) *B(LH,L) + X(K,LH)  * B(LH,LH)
                        MAT(4,2) = X(KH,L)*B(LH,L) + X(KH,LH) * B(LH,LH)

                        X(1:KE,L) = X(1:KE,L) - MAT(1,1) * A(1:KE,K) - MAT(2,1) * A(1:KE,KH) &
                            & + MAT(1,2) * B(1:KE,K) + MAT(2,2) * B(1:KE,KH)
                        X(1:KE,LH) = X(1:KE,LH) - MAT(3,1) * A(1:KE,K) - MAT(4,1) * A(1:KE,KH) &
                            & + MAT(3,2) * B(1:KE,K) + MAT(4,2) * B(1:KE,KH)
                    END IF

                END IF


                K = K - 1
            END DO

            ! Update January 2021
            IF ( LE .GT. 0 ) THEN
                IF (LB.EQ.1) THEN
                    CALL DGEMV("N", LE, LE, ONE, A(1,1),  LDA, X(1,L), IONE, ZERO, WORK, IONE)
                    CALL DSYR2("U", LE, -ONE, WORK,IONE, A(1, L), IONE, X(1,1), LDX)

                    WORK(1:LE) = X(L,L) * A(1:LE,L)
                    CALL DGER(LE, LE, -ONE, WORK, IONE, A(1,L),IONE, X(1,1), LDX)

                    CALL DGEMV("N", LE, LE, ONE, B(1,1),  LDB, X(1,L), IONE, ZERO, WORK, IONE)
                    CALL DSYR2("U", LE, ONE, WORK,IONE, B(1, L), IONE, X(1,1), LDX)

                    WORK(1:LE) = X(L,L) * B(1:LE,L)
                    CALL DGER(LE, LE, ONE, WORK, IONE, B(1,L),IONE, X(1,1), LDX)

                ELSE
                    CALL DGEMV("N", LE, LE, ONE, A(1,1),  LDA, X(1,L), IONE, ZERO, WORK, IONE)
                    CALL DSYR2("U", LE, -ONE, WORK,IONE, A(1, L), IONE, X(1,1), LDX)

                    CALL DGEMV("N", LE, LE, ONE, A(1,1),  LDA, X(1,LH), IONE, ZERO, WORK, IONE)
                    CALL DSYR2("U", LE, -ONE, WORK,IONE, A(1, LH), IONE, X(1,1), LDX)

                    WORK(1:LE) = X(L,L) * A(1:LE,L) + X(L,LH) * A(1:LE, LH)
                    CALL DGER(LE, LE, -ONE, WORK, IONE, A(1,L),IONE, X(1,1), LDX)

                    WORK(1:LE) = X(L,LH) * A(1:LE,L) + X(LH,LH) * A(1:LE, LH)
                    CALL DGER(LE, LE, -ONE, WORK, IONE, A(1,LH),IONE, X(1,1), LDX)

                    CALL DGEMV("N", LE, LE, ONE, B(1,1),  LDB, X(1,L), IONE, ZERO, WORK, IONE)
                    CALL DSYR2("U", LE, ONE, WORK,IONE, B(1, L), IONE, X(1,1), LDX)

                    CALL DGEMV("N", LE, LE, ONE, B(1,1),  LDB, X(1,LH), IONE, ZERO, WORK, IONE)
                    CALL DSYR2("U", LE, ONE, WORK,IONE, B(1, LH), IONE, X(1,1), LDX)

                    WORK(1:LE) = X(L,L) * B(1:LE,L) + X(L,LH) * B(1:LE, LH)
                    CALL DGER(LE, LE, ONE, WORK, IONE, B(1,L),IONE, X(1,1), LDX)

                    WORK(1:LE) = X(L,LH) * B(1:LE,L) + X(LH,LH) * B(1:LE, LH)
                    CALL DGER(LE, LE, ONE, WORK, IONE, B(1,LH),IONE, X(1,1), LDX)

                END IF
            END IF
            L = L - 1
        END DO


    ELSE
        !     ! (T, N)
        L = 1
        DO WHILE ( L .LE. M )
            LB = 1
            IF ( L .LT. M ) THEN
                IF ( A(L+1,L) .NE. ZERO ) THEN
                    LB =2
                END IF
            END IF
            LH = L + LB - 1

            K = L
            DO WHILE ( K .LE. M )
                KB = 1
                IF ( K .LT. M ) THEN
                    IF ( ( A(K+1,K) .NE. ZERO )) THEN
                        KB = 2
                    END IF
                END IF
                KH = K + KB - 1

                IF (K .EQ. L .AND. KB .EQ. 2 ) THEN
                    X(K, KH ) = X(KH,K)
                END IF

                IT = 0
                DO

                    ! Solve Inner System
                    IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 1

                        MAT(1,1) = A(L,L) * A(K,K) - B(L,L) * B(K,K)
                        RHS(1)   = X(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = A(K,K)
                        C11 = -B(K,K)

                        D11 = B(L,L)
                        D12 = B(L,LH)
                        D21 = B(LH,L)
                        D22 = B(LH,LH)

                        B11 = A(L,L)
                        B21 = A(LH,L)
                        B12 = A(L,LH)
                        B22 = A(LH,LH)


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
                        C11 = B(K,K)
                        C12 = B(K,KH)
                        C22 = B(KH,KH)

                        D11 = - B(L,L)
                        B11 = A(L,L)

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
                        C11 = - B(K,K)
                        C12 = - B(K,KH)
                        C22 = - B(KH,KH)

                        D11 = B(L,L)
                        D12 = B(L,LH)
                        D21 = B(LH,L)
                        D22 = B(LH,LH)

                        B11 = A(L,L)
                        B12 = A(L,LH)
                        B21 = A(LH, L)
                        B22 = A(LH,LH)

                        MAT(1,1) = B11*A11 + D11*C11
                        MAT(1,2) = B11*A21 !+ D11*C21
                        MAT(1,3) = B21*A11 + D21*C11
                        MAT(1,4) = B21*A21 !+ D21*C21

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
                        CALL DLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                        IF ( INFO1 .EQ. 0 ) EXIT
                        IT = 1
                    ELSE
                        CALL DGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                        IF ( INFO1 .NE. 0 ) INFO = INFO1
                        CALL DGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
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

                IF ( K.EQ.L .AND. KB.EQ.2) THEN
                    X(K,KH) = HALF*(RHS(2)+RHS(3))
                    X(KH,K) = X(K,KH)
                ELSE
                    ! Transpose causes a temporary here
                    ! X(L:LH,K:KH) = TRANSPOSE(X(K:KH,L:LH))
                    DO II = K, KH
                        DO JJ = L, LH
                           X(JJ, II) = X(II,JJ)
                        END DO
                       END DO
                END IF

                IF ( KH .LT. M ) THEN
                    IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                        X(KH+1:M, L) = X(KH+1:M, L) - A(K,KH+1:M) * (X(K,L) * A(L,L))  +  B(K,KH+1:M) * (X(K,L) * B(L,L))

                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                        X(KH+1:M, L) = X(KH+1:M, L) - A(K,KH+1:M) * (X(K,L) * A(L,L) + X(K,LH) * A(LH,L) ) &
                            & +  B(K,KH+1:M) * (X(K,L) * B(L,L)+X(K,LH)*B(LH,L))
                        X(KH+1:M, LH) = X(KH+1:M, LH) - A(K,KH+1:M) * (X(K,L) * A(L,LH) + X(K,LH) * A(LH,LH) ) &
                            & +  B(K,KH+1:M) * (X(K,L) * B(L,LH)+X(K,LH)*B(LH,LH))
                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                        X(KH+1:M,L) = X(KH+1:M,L) - A(K,KH+1:M)*(X(K,L)*A(L,L)) - A(KH,KH+1:M)*(X(KH,L)*A(L,L))
                        X(KH+1:M,L) = X(KH+1:M,L) + B(K,KH+1:M)*(X(K,L)*B(L,L)) + B(KH,KH+1:M)*(X(KH,L)*B(L,L))
                    ELSE
                        MAT(1,1) = X(K,L)*A(L,L) + X(K,LH) * A(LH,L)
                        MAT(2,1) = X(KH,L)*A(L,L) + X(KH,LH) * A(LH,L)
                        MAT(3,1) = X(K,L)*A(L,LH) + X(K,LH) * A(LH,LH)
                        MAT(4,1) = X(KH,L)*A(L,LH) + X(KH,LH) * A(LH,LH)

                        MAT(1,2) = X(K,L)*B(L,L) + X(K,LH) * B(LH,L)
                        MAT(2,2) = X(KH,L)*B(L,L) + X(KH,LH) * B(LH,L)
                        MAT(3,2) = X(K,L)*B(L,LH) + X(K,LH) * B(LH,LH)
                        MAT(4,2) = X(KH,L)*B(L,LH) + X(KH,LH) * B(LH,LH)


                        X(KH+1:M,L) = X(KH+1:M,L) - MAT(1,1)*A(K,KH+1:M) - MAT(2,1) * A(KH,KH+1:M) &
                            & +  MAT(1,2)*B(K,KH+1:M) + MAT(2,2) * B(KH,KH+1:M)
                        X(KH+1:M,LH) = X(KH+1:M,LH) - MAT(3,1)*A(K,KH+1:M) - MAT(4,1) * A(KH,KH+1:M) &
                            & +  MAT(3,2)*B(K,KH+1:M) + MAT(4,2) * B(KH,KH+1:M)
                    END IF
                END IF

                K = KH + 1
            END DO

            IF ( LH .LT. M ) THEN
                IF (LB .EQ. 1 ) THEN
                    CALL DGEMV("T", M-LH, M-LH, ONE, A(LH+1,LH+1), LDA, X(L,LH+1), LDX, ZERO, WORK, IONE)

                    WORK(M+1:M+M-LH) = A(L,LH+1:M)
                    CALL DSYR2("L", M-LH, -ONE, WORK(M+1),IONE, WORK, IONE, X(LH+1,LH+1), LDX)

                    ! This causes a temporary copy
                    ! WORK(1:M-LH) = X(L,L)* WORK(M+1:M+M-LH)
                    DO II = 1, M-LH
                        WORK(II) = X(L,L) * WORK(M+II)
                    END DO
                    CALL DGER(M-LH, M-LH,  -ONE, WORK(M+1:M+M-LH), IONE, WORK, IONE, X(LH+1,LH+1), LDX)

                    CALL DGEMV("T", M-LH, M-LH, ONE, B(LH+1,LH+1), LDB, X(L,LH+1), LDX, ZERO, WORK, IONE)

                    WORK(M+1:M+M-LH) = B(L,LH+1:M)
                    CALL DSYR2("L", M-LH, ONE, WORK(M+1),IONE, WORK, IONE, X(LH+1,LH+1), LDX)

                    ! WORK(1:M-LH) = X(L,L)* WORK(M+1:M+M-LH)
                    DO II = 1, M-LH
                        WORK(II) = X(L,L) * WORK(M+II)
                    END DO

                    CALL DGER(M-LH, M-LH,  ONE, WORK(M+1:M+M-LH), IONE, WORK, IONE, X(LH+1,LH+1), LDX)

                ELSE
                    I = M - LH
                    ! Need to be checked if it is really faster
                    ! IF ( I .LE. 256) THEN
                    !     WORKA(1:I) =  A(L,LH+1:M)
                    !     WORKB(1:I) =  A(LH,LH+1:M)
                    !
                    !     CALL DGEMV("T", M-LH, M-LH, ONE, A(LH+1,LH+1), LDA, X(L,LH+1), LDX, ZERO, WORK, IONE)
                    !     CALL DSYR2("L", M-LH, -ONE, WORKA, IONE, WORK, IONE, X(LH+1,LH+1), LDX)
                    !     CALL DGEMV("T", M-LH, M-LH, ONE, A(LH+1,LH+1), LDA, X(LH,LH+1), LDX, ZERO, WORK, IONE)
                    !     CALL DSYR2("L", M-LH, -ONE, WORKB, IONE, WORK, IONE, X(LH+1,LH+1), LDX)
                    !
                    !
                    !     WORKC(1:I) = X(L,L) * WORKA(1:I) + X(LH, L) * WORKB(1:I)
                    !     WORKD(1:I) = X(LH, L) * WORKA(1:I) + X(LH, LH) * WORKB(1:I)
                    !
                    !     CALL DGER(I, I, -ONE, WORKA, IONE, WORKC, IONE, X(LH+1, LH+1), LDX)
                    !     CALL DGER(I, I, -ONE, WORKB, IONE, WORKD, IONE, X(LH+1, LH+1), LDX)
                    !
                    !     WORKA(1:I) =  B(L,LH+1:M)
                    !     WORKB(1:I) =  B(LH,LH+1:M)
                    !
                    !     CALL DGEMV("T", M-LH, M-LH, ONE, B(LH+1,LH+1), LDB, X(L,LH+1), LDB, ZERO, WORK, IONE)
                    !     CALL DSYR2("L", M-LH, ONE, WORKA, IONE, WORK, IONE, X(LH+1,LH+1), LDX)
                    !     CALL DGEMV("T", M-LH, M-LH, ONE, B(LH+1,LH+1), LDB, X(LH,LH+1), LDB, ZERO, WORK, IONE)
                    !     CALL DSYR2("L", M-LH, ONE, WORKB, IONE, WORK, IONE, X(LH+1,LH+1), LDX)
                    !
                    !
                    !     WORKC(1:I) = X(L,L) * WORKA(1:I) + X(LH, L) * WORKB(1:I)
                    !     WORKD(1:I) = X(LH, L) * WORKA(1:I) + X(LH, LH) * WORKB(1:I)
                    !
                    !     CALL DGER(I, I, ONE, WORKA, IONE, WORKC, IONE, X(LH+1, LH+1), LDX)
                    !     CALL DGER(I, I, ONE, WORKB, IONE, WORKD, IONE, X(LH+1, LH+1), LDX)
                    !
                    !
                    ! ELSE
                        CALL DGEMV("T", M-LH, M-LH, ONE, A(LH+1,LH+1), LDA, X(L,LH+1), LDX, ZERO, WORK, IONE)
                        WORK(M+1:M+I) = A(L,LH+1:M)
                        CALL DSYR2("L", M-LH, -ONE, WORK(M+1), IONE, WORK, IONE, X(LH+1,LH+1), LDX)
                        CALL DGEMV("T", M-LH, M-LH, ONE, A(LH+1,LH+1), LDA, X(LH,LH+1), LDX, ZERO, WORK, IONE)
                        WORK(M+1:M+I) = A(LH,LH+1:M)
                        CALL DSYR2("L", M-LH, -ONE, WORK(M+1), IONE, WORK, IONE, X(LH+1,LH+1), LDX)


                        WORK(1:I) = X(L,L) * A(L,LH+1:M) + X(LH, L) * A(LH, LH+1:M)
                        WORK(M+1:M+I) = X(LH, L) * A(L,LH+1:M) + X(LH, LH) * A(LH, LH+1:M)

                        CALL DGER(I, I, -ONE, A(L,LH+1), LDA, WORK(1), IONE, X(LH+1, LH+1), LDX)
                        CALL DGER(I, I, -ONE, A(LH,LH+1), LDA, WORK(M+1), IONE, X(LH+1, LH+1), LDX)

                        CALL DGEMV("T", M-LH, M-LH, ONE, B(LH+1,LH+1), LDB, X(L,LH+1), LDB, ZERO, WORK, IONE)
                        WORK(M+1:M+I) = B(L,LH+1:M)
                        CALL DSYR2("L", M-LH, ONE, WORK(M+1), IONE, WORK, IONE, X(LH+1,LH+1), LDX)
                        CALL DGEMV("T", M-LH, M-LH, ONE, B(LH+1,LH+1), LDB, X(LH,LH+1), LDB, ZERO, WORK, IONE)
                        WORK(M+1:M+I) = B(LH,LH+1:M)
                        CALL DSYR2("L", M-LH, ONE, WORK(M+1), IONE, WORK, IONE, X(LH+1,LH+1), LDX)


                        WORK(1:I) = X(L,L) * B(L,LH+1:M) + X(LH, L) * B(LH, LH+1:M)
                        WORK(M+1:M+I) = X(LH, L) * B(L,LH+1:M) + X(LH, LH) * B(LH, LH+1:M)

                        CALL DGER(I, I, ONE, B(L,LH+1), LDB, WORK(1), IONE, X(LH+1, LH+1), LDX)
                        CALL DGER(I, I, ONE, B(LH,LH+1), LDB, WORK(M+1), IONE, X(LH+1, LH+1), LDX)

                    ! END IF
                END IF
            END IF
            L = LH + 1
        END DO

    END IF

END SUBROUTINE



