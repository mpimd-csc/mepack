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
!
!> \brief Level-2 Bartels-Stewart Algorithm for the Sylvester equation (optimized)
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TRSYLV_L2 ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANSA(1), TRANSB(1)
!           REAL SGN, SCALE
!           INTEGER M, N, LDA, LDB, LDX, INFO
!           REAL A(LDA, *), B(LDB, *) , X(LDX, *)
!           REAL WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TRSYLV_L2 solves a Sylvester equation of the following forms
!>
!>    op1(A) * X  +  X * op2(B) = SCALE * Y                              (1)
!>
!> or
!>
!>    op1(A) * X  -  X * op2(B) = SCALE * Y                              (2)
!>
!> where A is a M-by-M quasi upper triangular matrix, B is a N-by-N quasi upper triangular
!> matrix. The right hand side Y and the solution X are M-by-N matrices.  Typically the matrices
!> A and B are generated via SGEES form LAPACK.
!> \endverbatim
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
!> \li Replaced level-3 by level-2 and level-1 calls,
!> \li Replaced SAXPY operation by Fortran intrinsic
!> \li Replaced all BLAS calls except of STRMV and SGER by Fortran intrinsic,
!
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup sgltrsylv
!
SUBROUTINE SLA_TRSYLV_L2 ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    REAL SGN, SCALE
    INTEGER M, N, LDA, LDB, LDX, INFO
    REAL A(LDA, *), B(LDB, *) ,  X(LDX, *)
    REAL WORK(*)


    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE, IT
    INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
    REAL MAT(4,4), RHS(4)
    REAL A11, A12, A21, A22
    REAL B11, B12, B21, B22
    REAL SCAL

    LOGICAL BTRANSA, BTRANSB
    REAL ZERO , ONE
    PARAMETER (ZERO = 0.0, ONE = 1.0)
    INTEGER IONE
    PARAMETER (IONE = 1)

    ! EXTERNAL FUNCTIONS
    EXTERNAL SGER
    EXTERNAL SGESC2
    EXTERNAL SGETC2X
    EXTERNAL SLA_SMALL_SOLVE4
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
        CALL XERROR_HANDLER('SLA_TRSYLV_L2', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1) THEN
        INFO = 1
        RETURN
    END IF


    ! Quick Return
    SCALE = ONE
    IF ( M .EQ. 0 .OR. N .EQ. 0 ) THEN
        RETURN
    END IF

    ! Prepare
    INFO = 0
    WORK(1)=WORK(1)

    IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB ) THEN
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
                ! First half AXB
                IF ( KB .EQ. 1 ) THEN
                    CALL SGER(KE, N, -ONE, A(1,K) , IONE, X(K,1), LDX, X(1,1), LDX)
                ELSE
                    CALL SGER(KE, N, -ONE, A(1,K) , IONE, X(K,1), LDX, X(1,1), LDX)
                    CALL SGER(KE, N, -ONE, A(1,KH) , IONE, X(KH,1), LDX, X(1,1), LDX)
                ENDIF
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
                    IF (B(L,L-1) .NE. ZERO) THEN
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

                        MAT(1,1) = A(K,K) + SGN * B(L,L)
                        RHS(1)   = X(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = A(K,K)

                        B11 = SGN*B(L,L)
                        B21 = SGN*B(LH,L)
                        B12 = SGN*B(L,LH)
                        B22 = SGN*B(LH,LH)


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

                        B11 = SGN * B(L,L)

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

                        B11 = SGN*B(L,L)
                        B12 = SGN*B(L,LH)
                        B21 = SGN*B(LH, L)
                        B22 = SGN*B(LH,LH)

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
                    IF ( KB .EQ. IONE .AND. LB .EQ. IONE ) THEN
                        X(K,1:LE) =  X(K,1:LE) - SGN*X(K,L)*B(1:LE, L)

                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE ) THEN
                        X(K,1:LE) =  X(K,1:LE) -SGN*X(K,L)*B(1:LE, L) - SGN*X(K,LH)*B(1:LE, LH)

                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE ) THEN
                        X(K, 1:LE)  = X(K,1:LE) - SGN* X(K,L)  *  B(1:LE, L)
                        X(KH, 1:LE) = X(KH,1:LE) - SGN* X(KH,L) * B(1:LE, L)

                    ELSE
                        X(K, 1:LE) = X(K,1:LE)   - SGN * X(K,L) * B(1:LE, L) - SGN * X(K,LH) * B(1:LE, LH)
                        X(KH, 1:LE) = X(KH,1:LE) - SGN * X(KH,L)* B(1:LE, L) - SGN * X(KH,LH) * B(1:LE, LH)
                    END IF
                END IF

                L = L - 1
            END DO

            ! Update January 2021
            IF ( K .GT. 1 ) THEN
                IF ( KB .EQ. IONE ) THEN
                    CALL SGER(KE, N, -ONE, A(1,K), IONE, X(K,1), LDX, X(1,1),LDX)
                ELSE
                    CALL SGER(KE, N, -ONE, A(1,K), IONE, X(K,1), LDX, X(1,1),LDX)
                    CALL SGER(KE, N, -ONE, A(1,KH), IONE, X(KH,1), LDX, X(1,1),LDX)
                END IF

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
                    IF ( B(L+1,L) .NE. ZERO ) THEN
                        LB =2
                    END IF
                END IF
                LH = L + LB - 1

                IT = 0
                DO

                    ! Solve Inner System
                    IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 1

                        MAT(1,1) = A(K,K) + SGN * B(L,L)
                        RHS(1)   = X(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = A(K,K)


                        B11 = SGN * B(L,L)
                        B21 = SGN * B(LH,L)
                        B12 = SGN * B(L,LH)
                        B22 = SGN * B(LH,LH)


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

                        B11 = SGN * B(L,L)

                        MAT(1,1) = A11 + B11
                        MAT(1,2) = A21
                        MAT(2,1) = A12
                        MAT(2,2) = A22 + B11

                        RHS(1) = X(K,L)
                        RHS(2) = X(KH,L)

                    ELSE
                        DIMMAT = 4

                        A11 = A(K,K)
                        A12 = A(K,KH)
                        A21 = A(KH,K)
                        A22 = A(KH,KH)


                        B11 = SGN * B(L,L)
                        B12 = SGN * B(L,LH)
                        B21 = SGN * B(LH, L)
                        B22 = SGN * B(LH,LH)

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
                    IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                        X(K,LH+1:N) = X(K,LH+1:N) - SGN * X(K,L) * B(L,LH+1:N)
                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                        X(K,LH+1:N) = X(K,LH+1:N) - SGN*X(K,L) * B(L,LH+1:N) - SGN*X(K,LH) * B(LH,LH+1:N)
                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                        X(K,LH+1:N) = X(K,LH+1:N) - SGN * X(K,L)   * B(L,LH+1:N)
                        X(KH,LH+1:N) = X(KH,LH+1:N)  - SGN * X(KH,L)  * B(L,LH+1:N)

                    ELSE
                        X(K,LH+1:N) = X(K,LH+1:N) - SGN * X(K,L) * B(L,LH+1:N) - SGN * X(K,LH) * B(LH, LH+1:N)
                        X(KH,LH+1:N) = X(KH,LH+1:N) - SGN * X(KH,L) * B(L,LH+1:N) - SGN * X(KH,LH) * B(LH, LH+1:N)

                    END IF
                END IF

                L = LH + 1
            END DO

            IF ( KH .LT. M ) THEN
                IF (KB .EQ. IONE ) THEN
                    CALL SGER(M-KH, N, -ONE, A(K,KH+1), LDA, X(K,1), LDX, X(KH+1, 1), LDX)
                ELSE
                    CALL SGER(M-KH, N, -ONE, A(K,KH+1), LDA, X(K,1), LDX, X(KH+1, 1), LDX)
                    CALL SGER(M-KH, N, -ONE, A(KH,KH+1), LDA, X(KH,1), LDX, X(KH+1, 1), LDX)
                END IF
            END IF
            K = KH + 1
        END DO


    ELSE IF ( BTRANSA .AND. BTRANSB) THEN
        ! (T, T)
        K = 1
        DO WHILE ( K .LE. M )
            KB = 1
            IF ( K .LT. M) THEN
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
                    IF (B(L,L-1) .NE. ZERO) THEN
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

                        MAT(1,1) = A(K,K) + SGN * B(L,L)
                        RHS(1)   = X(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = A(K,K)

                        B11 = SGN * B(L,L)
                        B21 = SGN * B(LH,L)
                        B12 = SGN * B(L,LH)
                        B22 = SGN * B(LH,LH)


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

                        B11 = SGN * B(L,L)

                        MAT(1,1) = A11 + B11
                        MAT(1,2) = A21
                        MAT(2,1) = A12
                        MAT(2,2) = A22 + B11

                        RHS(1) = X(K,L)
                        RHS(2) = X(KH,L)

                    ELSE
                        DIMMAT = 4

                        A11 = A(K,K)
                        A12 = A(K,KH)
                        A21 = A(KH,K)
                        A22 = A(KH,KH)

                        B11 = SGN * B(L,L)
                        B12 = SGN * B(L,LH)
                        B21 = SGN * B(LH, L)
                        B22 = SGN * B(LH,LH)

                        MAT(1,1) = A11 + B11
                        MAT(1,2) = A21
                        MAT(1,3) =         + B12
                        MAT(1,4) = ZERO

                        MAT(2,1) = A12
                        MAT(2,2) = A22 + B11
                        MAT(2,3) =  ZERO
                        MAT(2,4) =  B12

                        MAT(3,1) =       B21
                        MAT(3,2) = ZERO
                        MAT(3,3) = A11 + B22
                        MAT(3,4) = A21

                        MAT(4,1) = ZERO
                        MAT(4,2) =       B21
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
                IF ( L .GT. 0  ) THEN
                    IF ( KB .EQ. IONE .AND. LB .EQ. IONE ) THEN
                        X(K,1:LE) =  X(K,1:LE) - SGN*X(K,L)*B(1:LE, L)

                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE ) THEN
                        X(K,1:LE) =  X(K,1:LE) -SGN*X(K,L)*B(1:LE, L) - SGN*X(K,LH)*B(1:LE, LH)
                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE ) THEN
                        X(K, 1:LE) = X(K,1:LE) - SGN*X(K,L) * B(1:LE, L)
                        X(KH, 1:LE) = X(KH,1:LE) - SGN*X(KH,L) * B(1:LE, L)

                    ELSE
                        MAT(1,2) = X(K,L)
                        MAT(2,2) = X(KH,L)
                        MAT(3,2) = X(K,LH)
                        MAT(4,2) = X(KH,LH)


                        X(K, 1:LE) = X(K,1:LE) - SGN * X(K,L) * B(1:LE, L) - SGN * X(K,LH) * B(1:LE, LH)
                        X(KH, 1:LE) = X(KH,1:LE)  - SGN * X(KH,L) * B(1:LE, L) - SGN * X(KH,LH) * B(1:LE, LH)
                    END IF
                END IF

                L = L - 1
            END DO

            ! Update January 2021
            IF ( KH .LT. M ) THEN
                IF ( KB .EQ. IONE ) THEN
                    CALL SGER(M-KH, N, -ONE, A(K,KH+1), LDA, X(K,1), LDX, X(KH+1, 1), LDX)
                ELSE
                    CALL SGER(M-KH, N, -ONE, A(K,KH+1), LDA, X(K,1), LDX, X(KH+1, 1), LDX)
                    CALL SGER(M-KH, N, -ONE, A(KH,KH+1), LDA, X(KH,1), LDX, X(KH+1, 1), LDX)
                END IF
            END IF

            K = KH + 1
        END DO

    END IF
END SUBROUTINE



