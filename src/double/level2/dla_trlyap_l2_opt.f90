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

!> \brief Level-2 Bartels-Stewart Algorithm for the Lyapunov Equation (Optimized)
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TRLYAP_L2_OPT ( TRANSA, M, A, LDA, X, LDX, SCALE, WORK, INFO)
!           CHARACTER(1) TRANSA
!           DOUBLE PRECISION SCALE
!           INTEGER M, LDA, LDX, INFO
!           DOUBLE PRECISION A(LDA, *), X(LDX, *)
!           DOUBLE PRECISION WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLA_TRLYAP_L2 solves a Lyapunov equation of the following forms
!>
!>    A * X  +  X * A**T = SCALE * Y                              (1)
!>
!> or
!>
!>    A **T * X  -  X * A = SCALE * Y                              (2)
!>
!> where A is a M-by-M quasi upper triangular matrix.
!> The right hand side Y and the solution X are M-by-N matrices.  Typically the matrix A
!> is generated by DGEES from LAPACK.
!> \endverbatim
!> \remark The algorithm is implemented using BLAS level 2 operations.
!> \remark The transposed case (2) is optimized w.r.t. to the usage of the DSYR2 operation.
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
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and C.  M >= 0.
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
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,N)
!>          On input, the matrix X contains the right hand side Y.
!>          On output, the matrix X contains the solution of Equation (1) or (2)
!>          Right hand side Y and the solution X are symmetric M-by-M matrices.
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
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup dbltrlyap
!
SUBROUTINE DLA_TRLYAP_L2_OPT( TRANSA, M, A, LDA, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1)
    DOUBLE PRECISION SCALE
    INTEGER M, LDA, LDX, INFO
    DOUBLE PRECISION A(LDA, *), X(LDX, *)
    DOUBLE PRECISION WORK(*)


    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE, IT,II,JJ
    INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
    DOUBLE PRECISION MAT(4,4), RHS(4)
    DOUBLE PRECISION A11, A12, A21, A22
    DOUBLE PRECISION B11, B12, B21, B22
    DOUBLE PRECISION SCAL
    DOUBLE PRECISION WORKI(64), WORKII(64)

    DOUBLE PRECISION AL(64,64), XL(64,64)
    !dir$ attributes align: 64:: WORKI
    !dir$ attributes align: 64:: MAT
    !dir$ attributes align: 64:: RHS
    !dir$ attributes align: 64:: IPIV
    !dir$ attributes align: 64:: JPIV
    !dir$ attributes align: 64:: AL
    !dir$ attributes align: 64:: XL

    LOGICAL BTRANSA
    DOUBLE PRECISION ZERO , ONE
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
    INTEGER IONE
    PARAMETER (IONE = 1)
    INTEGER :: J

    ! EXTERNAL FUNCTIONS
    EXTERNAL DGESC2
    EXTERNAL DGETC2X
    EXTERNAL DLA_SMALL_SOLVE4
    EXTERNAL DLA_TRLYAP_L2
    EXTERNAL DSYR2
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME


    IF ( M.GT.64) THEN
        CALL DLA_TRLYAP_L2( TRANSA, M, A, LDA, X, LDX, SCALE, WORK, INFO)
        RETURN
    END IF

    ! Check Input
    ININFO = INFO
    INFO = 0
    BTRANSA = LSAME(TRANSA, 'T')

    IF ( .NOT. BTRANSA .AND. .NOT. LSAME(TRANSA, 'N')) THEN
        INFO = -1
    ELSE IF ( M .LT. 0 .OR. M .GT. 64) THEN
        INFO = -2
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -4
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -6
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('DLA_TRLYAP_L2_OPT', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1) THEN
        INFO = 1
        RETURN
    END IF

    ! Quick Return
    SCALE = ONE
    IF ( M .EQ. 0 ) THEN
        RETURN
    END IF

    ! Prepare
    INFO = 0
    WORK(1)=WORK(1)

    IF ( .NOT. BTRANSA  ) THEN
        ! (N,T)
        XL(1:M,1:M) = X(1:M,1:M)
        AL(1:M,1:M) = A(1:M,1:M)

        L = M
        DO WHILE ( L .GT. 0 )
            LB = 1
            LH = L
            IF ( L .GT. 1 ) THEN
                IF (AL(L,L-1) .NE. ZERO) THEN
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
                    IF ( AL(K,K-1) .NE. ZERO ) THEN
                        KB = 2
                        K = K - 1
                    END IF
                END IF
                KE = K - 1

                IT = 0

                DO

                    ! Solve Inner System
                    IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 1

                        MAT(1,1) = AL(K,K) + AL(L,L)
                        RHS(1)   = XL(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)

                        B11 = AL(L,L)
                        B21 = AL(LH,L)
                        B12 = AL(L,LH)
                        B22 = AL(LH,LH)


                        MAT(1,1) = A11 + B11
                        MAT(1,2) =       B12
                        MAT(2,1) =       B21
                        MAT(2,2) = A11 + B22

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(K,LH)

                    ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)
                        A12 = AL(K,KH)
                        A21 = AL(KH,K)
                        A22 = AL(KH,KH)

                        B11 = AL(L,L)

                        MAT(1,1) = A11 + B11
                        MAT(1,2) = A12
                        MAT(2,1) = A21
                        MAT(2,2) = A22 + B11

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)

                    ELSE
                        DIMMAT = 4
                        IF ( K.EQ. L) XL(KH, K ) = XL(K,KH)

                        A11 = AL(K,K)
                        A12 = AL(K,KH)
                        A21 = AL(KH,K)
                        A22 = AL(KH,KH)

                        B11 = AL(L,L)
                        B12 = AL(L,LH)
                        B21 = AL(LH, L)
                        B22 = AL(LH,LH)

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

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)
                        RHS(3) = XL(K, LH)
                        RHS(4) = XL(KH,LH)

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
                    XL(K,L) = RHS(1)
                    XL(L,K) = RHS(1)
                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    XL(K,L)  = RHS(1)
                    XL(K,LH) = RHS(2)
                    XL(L,K)  = RHS(1)
                    XL(LH,K) = RHS(2)
                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    XL(K,L) = RHS(1)
                    XL(KH,L) = RHS(2)
                    XL(L,K) = RHS(1)
                    XL(L,KH) = RHS(2)
                ELSE
                    XL(K,L)   = RHS(1)
                    XL(KH,L)  = RHS(2)
                    XL(K, LH) = RHS(3)
                    XL(KH,LH) = RHS(4)

                    XL(L, K)  = RHS(1)
                    XL(L, KH) = RHS(2)
                    XL(LH,K)  = RHS(3)
                    XL(LH,KH) = RHS(4)

                END IF

                IF ( K .GT. 1 ) THEN
                    IF ( KB .EQ. IONE .AND. LB .EQ. IONE ) THEN
                        XL(1:KE, L ) = XL(1:KE, L) - XL(K,L)* AL(1:KE, K)
                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                        XL(1:KE,L) = XL(1:KE,L) - XL(K,L) * AL(1:KE,K)
                        XL(1:KE,LH) = XL(1:KE,LH) - XL(K,LH) *AL(1:KE,K)
                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                        XL(1:KE, L) = XL(1:KE,L) - AL(1:KE,K) *XL(K,L) - AL(1:KE,KH) *XL(KH,L)
                    ELSE
                        XL(1:KE,L) = XL(1:KE,L) -   XL(K,L)  * AL(1:KE,K) - XL(KH,L) * AL(1:KE,KH)
                        XL(1:KE,LH) = XL(1:KE,LH) - XL(K,LH) * AL(1:KE,K) - XL(KH,LH) * AL(1:KE,KH)
                    END IF

                END IF

                K = K - 1
            END DO

            ! Update January 2021
            IF ( L .GT. 1 ) THEN
                IF ( LB .EQ. IONE ) THEN
                    DO J = 1, L-1
                       XL(1:J,J) = XL(1:J,J) - XL(1:J,L) * AL(J,L) - AL(1:J,L) * XL(J,L)
                    END DO
                    ! CALL DSYR2("U", L-1, -ONE, XL(1,L), IONE, AL(1,L), IONE, XL(1,1), 64)
                ELSE
                    DO J = 1, L-1
                       XL(1:J,J) = XL(1:J,J) - XL(1:J,L) * AL(J,L) - AL(1:J,L) * XL(J,L)  &
                           & - XL(1:J,LH) * AL(J,LH) - AL(1:J,LH) * XL(J,LH)
                    END DO
                    ! CALL DSYR2("U", L-1, -ONE, XL(1,L), IONE, AL(1,L), IONE, XL(1,1), 64)
                    ! CALL DSYR2("U", L-1, -ONE, XL(1,LH), IONE, AL(1,LH), IONE, XL(1,1), 64)
                ENDIF
            END IF
            L = L - 1
        END DO

        X(1:M,1:M) = XL(1:M,1:M)
    ELSE
        ! (T, N)
        XL(1:M,1:M) = X(1:M,1:M)
        AL(1:M,1:M) = A(1:M,1:M)


        L = 1
        DO WHILE ( L .LE. M )
            LB = 1
            IF ( L .LT. M ) THEN
                IF ( AL(L+1,L) .NE. ZERO ) THEN
                    LB =2
                END IF
            END IF
            LH = L + LB - 1

            K = L
            DO WHILE ( K .LE. M )
                KB = 1
                IF ( K .LT. M ) THEN
                    IF ( ( AL(K+1,K) .NE. ZERO )) THEN
                        KB = 2
                    END IF
                END IF
                KH = K + KB - 1

                IF (K.EQ.L .AND. KB.EQ.2) THEN
                    XL(K,KH) = XL(KH,K)
                END IF

                IT = 0
                DO

                    ! Solve Inner System
                    IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 1

                        MAT(1,1) = AL(K,K) + AL(L,L)
                        RHS(1)   = XL(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)


                        B11 = AL(L,L)
                        B21 = AL(LH,L)
                        B12 = AL(L,LH)
                        B22 = AL(LH,LH)


                        MAT(1,1) = A11 + B11
                        MAT(1,2) =       B21
                        MAT(2,1) =       B12
                        MAT(2,2) = A11 + B22

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(K,LH)

                    ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)
                        A12 = AL(K,KH)
                        A21 = AL(KH,K)
                        A22 = AL(KH,KH)

                        B11 = AL(L,L)

                        MAT(1,1) = A11 + B11
                        MAT(1,2) = A21
                        MAT(2,1) = A12
                        MAT(2,2) = A22 + B11

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)

                    ELSE
                        DIMMAT = 4

                        A11 = AL(K,K)
                        A12 = AL(K,KH)
                        A21 = AL(KH,K)
                        A22 = AL(KH,KH)


                        B11 = AL(L,L)
                        B12 = AL(L,LH)
                        B21 = AL(LH, L)
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

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)
                        RHS(3) = XL(K, LH)
                        RHS(4) = XL(KH,LH)

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
                            XL(1:M, 1:M) = SCAL*XL(1:M, 1:M)
                            SCALE = SCALE * SCAL
                            INFO = 1
                        END IF
                        EXIT
                    END IF

                END DO

                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    XL(K,L) = RHS(1)
                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    XL(K,L) = RHS(1)
                    XL(K,LH) = RHS(2)
                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    XL(K,L) = RHS(1)
                    XL(KH,L) = RHS(2)
                ELSE
                    XL(K,L) = RHS(1)
                    XL(KH,L) = RHS(2)
                    XL(K, LH) = RHS(3)
                    XL(KH,LH) = RHS(4)
                END IF

                IF ( K.NE.L) THEN
                    ! Transpose causes a temporary here
                    ! XL(L:LH,K:KH) = TRANSPOSE(XL(K:KH,L:LH))
                     DO II = K, KH
                        DO JJ = L, LH
                        XL(JJ, II) = XL(II,JJ)
                END DO
            END DO

                END IF

                ! in current row

                IF ( KH .LT. M ) THEN
                    IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                        XL(KH+1:M, L) = XL(KH+1:M, L) - AL(K,KH+1:M) * XL(K,L)
                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                        XL(KH+1:M, L) = XL(KH+1:M, L) - AL(K,KH+1:M) * XL(K,L)
                        XL(KH+1:M, LH) = XL(KH+1:M, LH) - AL(K,KH+1:M) * XL(K,LH)

                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                        XL(KH+1:M,L) = XL(KH+1:M,L) - AL(K,KH+1:M)*XL(K,L) - AL(KH,KH+1:M)*XL(KH,L)
                    ELSE
                        XL(KH+1:M,L) = XL(KH+1:M,L) -  XL(K,L) *AL(K,KH+1:M) - XL(KH,L) * AL(KH,KH+1:M)
                        XL(KH+1:M,LH) = XL(KH+1:M,LH) - XL(K,LH)*AL(K,KH+1:M) - XL(KH,LH) * AL(KH,KH+1:M)

                    END IF
                END IF



                K = KH + 1
            END DO

            IF ( LH .LT. M ) THEN
                IF ( LB .EQ. IONE ) THEN
                    WORKI(LH+1:M) = AL(L,LH+1:M)

                    ! CALL DSYR2("L", M-LH, -ONE, XL(LH+1,L), IONE, WORKI, IONE, XL(LH+1, LH+1), 64)

                    DO J = LH+1,M
                        XL(J:M,J) = XL(J:M,J) - XL(J:M,L) * WORKI(J) - WORKI(J:M) * XL(J,L)
                    END DO

                ELSE
                    WORKI(LH+1:M) = AL(L,LH+1:M)
                    WORKII(LH+1:M) = AL(LH,LH+1:M)

                    DO J = LH+1,M
                        XL(J:M,J) = XL(J:M,J) - XL(J:M,L) * WORKI(J) - WORKI(J:M) * XL(J,L) &
                            & - XL(J:M,LH) * WORKII(J) - WORKII(J:M) * XL(J,LH)
                    END DO


                    ! WORKI(1:M-LH) = AL(L,LH+1:M)
                    ! CALL DSYR2("L", M-LH, -ONE, XL(LH+1,L), IONE, WORKI, IONE, XL(LH+1, LH+1), 64)
                    ! WORKI(1:M-LH) = AL(LH,LH+1:M)
                    ! CALL DSYR2("L", M-LH, -ONE, XL(LH+1,LH),IONE, WORKI, IONE, XL(LH+1, LH+1), 64)
                END IF
            END IF

            L = LH + 1
        END DO

        X(1:M,1:M) = XL(1:M,1:M)

    END IF
END SUBROUTINE



