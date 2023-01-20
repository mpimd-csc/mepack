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


!> \brief Level-3 Bartels-Stewart Algorithm for the standard Lyapunov Equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TRLYAP_L3 ( TRANS, M, A, LDA, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANS
!           DOUBLE PRECISION SCALE
!           INTEGER M, LDA, LDX, INFO
!           DOUBLE PRECISION A(LDA, *), X(LDX, *)
!           DOUBLE PRECISION WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> DLA_TRLYAP_L3 solves a Lyapunov equation of the following forms
!>
!>    A  * X  +  X * A**T = SCALE * Y                              (1)
!>
!> or
!>
!>    A ** T * X  +  X * A = SCALE * Y                              (2)
!>
!> where A is a M-by-M quasi upper triangular matrix.
!> The right hand side Y and the solution X are M-by-N matrices.  Typically the matrix A
!> is generated by DGEES from LAPACK.
!> \endverbatim
!> \remark The algorithm is implemented using BLAS level 3 operations.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER
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
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Januar 2023
!> \ingroup dbltrlyap
!
SUBROUTINE DLA_TRLYAP_L3 ( TRANS, M, A, LDA, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    USE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
    USE MEPACK_OPTIONS_ISOLVER
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANS(1)
    DOUBLE PRECISION SCALE
    INTEGER M, LDA, LDX, INFO
    DOUBLE PRECISION A(LDA, *), X(LDX, *)
    DOUBLE PRECISION WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE
    INTEGER INFO1
    DOUBLE PRECISION SCAL
    INTEGER MB, ISOLVER

    LOGICAL BTRANS
    CHARACTER TRANSA, TRANSB
    DOUBLE PRECISION ZERO , ONE, SGN, HALF
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, SGN = 1.0D0, HALF = 0.5D0 )
    INTEGER IZERO, IONE, I, J
    PARAMETER(IZERO = 0, IONE = 1 )

    ! EXTERNAL FUNCTIONS
    EXTERNAL DCOPY
    EXTERNAL DGEMM
    EXTERNAL DLACPY
    EXTERNAL DLA_ITRANSPOSE
    EXTERNAL DLA_TRLYAP_L2
    EXTERNAL DLA_TRLYAP_L2_OPT

    EXTERNAL DLA_TRLYAP_RECURSIVE
    EXTERNAL DLA_TRSYLV_L2
    EXTERNAL DLA_TRSYLV_L2_LOCAL_COPY
    EXTERNAL DLA_TRSYLV_L2_LOCAL_COPY_128
    EXTERNAL DLA_TRSYLV_L2_LOCAL_COPY_32
    EXTERNAL DLA_TRSYLV_L2_LOCAL_COPY_64
    EXTERNAL DLA_TRSYLV_L2_LOCAL_COPY_96
    EXTERNAL DLA_TRSYLV_L2_REORDER
    EXTERNAL DLA_TRSYLV_L2_UNOPT
    EXTERNAL DLA_TRSYLV_RECURSIVE
    EXTERNAL DSYR2K
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME

    MB = TRLYAP_BLOCKSIZE(M)
    ISOLVER = TRLYAP_ISOLVER()

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
    ELSE IF ( MB .LT. 2 ) THEN
        INFO = -30
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('DLA_TRLYAP_L3', -INFO)
        RETURN
    END IF


    IF ( ININFO .EQ. -1 ) THEN
        INFO =  MB*MB
        RETURN
    END IF



    ! Quick Return
    SCALE = ONE
    SCAL = ONE
    INFO = 0
    INFO1 = 0

    IF ( BTRANS ) THEN
        TRANSA = 'N'
        TRANSB = 'T'
    ELSE
        TRANSA = 'T'
        TRANSB = 'N'
    END IF

    IF ( M .EQ. 0 ) THEN
        RETURN
    END IF

    IF ( M .LE. MB ) THEN
        KB = M
        LB = M
        CALL DLA_TRLYAP_L2(TRANSA, KB, A, LDA,  X, LDX, SCAL, WORK, INFO1)

        IF ( INFO1.NE. 0 ) THEN
            INFO = INFO1
        END IF
        IF ( SCAL .NE. ONE) THEN
            SCALE = SCAL
        END IF

        RETURN
    END IF

    ! Prepare
    INFO = 0

    IF ( BTRANS  ) THEN
        ! (N,T)

        L = M
        DO WHILE ( L .GT. 0 )
            IF (L .EQ. M .AND. MOD(M,MB) .NE. IZERO) THEN
                LH = L
                L = (L/MB)*MB + 1
            ELSE
                LH = L
                L = MAX(1, LH - MB + 1)
            END IF

            LB = LH - L + 1
            IF ( L .GT. 1 ) THEN
                IF ( A(L,L-1) .NE. ZERO ) THEN
                    IF ( LB .EQ. 1 ) THEN
                        LB = LB + 1
                        L = L - 1
                    ELSE
                        LB = LB - 1
                        L = L + 1
                    END IF
                END IF
            END IF
            LE = L - 1


            K = LH
            DO WHILE (K .GT. 0)
                IF (K .EQ. M .AND. MOD(M,MB) .NE. IZERO) THEN
                    KH = K
                    K = (M/MB)*MB + 1
                ELSE
                    KH = K
                    K = MAX(1, KH - MB + 1)
                END IF

                KB = KH - K + 1
                IF ( K .GT. 1 ) THEN
                    IF ( A(K,K-1) .NE. ZERO ) THEN
                        IF ( KB.EQ.1) THEN
                            KB = KB + 1
                            K = K - 1
                        ELSE
                            KB = KB - 1
                            K = K + 1
                        END IF
                    END IF
                END IF
                KE = K - 1


                IF ( K .EQ. L ) THEN

                    IF ( ISOLVER .EQ. 6 ) THEN
                        DO I = K,KH-1
                            CALL DCOPY(KH-I,X(I,I+1),LDX,X(I+1,I),IONE)
                        END DO

                        CALL DLA_TRLYAP_RECURSIVE(TRANSA, KB, A(K,K), LDA,  X(K,K), LDX, SCAL, WORK, INFO1)
                    ELSE
                        IF ( KB .LE. 64 ) THEN
                            CALL DLA_TRLYAP_L2_OPT(TRANSA, KB, A(K,K), LDA,  X(K,K), LDX, SCAL, WORK, INFO1)
                        ELSE
                            CALL DLA_TRLYAP_L2(TRANSA, KB, A(K,K), LDA,  X(K,K), LDX, SCAL, WORK, INFO1)
                        END IF
                        ! CALL DLA_TRLYAP_L2_OPT(TRANSA, KB, A(K,K), LDA,  X(K,K), LDX, SCAL, WORK, INFO1)
                    END IF

                ELSE
                    SELECT CASE(ISOLVER)
                    CASE(1)
                        IF ( MB .LE. 32 ) THEN
                            CALL DLA_TRSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( MB .LE. 64) THEN
                            CALL DLA_TRSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( MB .LE. 96) THEN
                            CALL DLA_TRSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( MB .LE. 128) THEN
                            CALL DLA_TRSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE
                            CALL DLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        END IF
                    CASE(2)
                        IF ( KB .GT. 128 .OR. LB .GT. 128) THEN
                            CALL DLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE
                            CALL DLA_TRSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        END IF
                    CASE(3)
                        CALL DLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            & X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(4)
                        CALL DLA_TRSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            & X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(5)
                        CALL DLA_TRSYLV_L2_UNOPT(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            & X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(6)
                        CALL DLA_TRSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            & X(K,L), LDX, SCAL, WORK, INFO1)
                    END SELECT

                    DO J = K,KH
                        DO  I = L, LH
                            X(I,J) = X(J,I)
                        END DO
                    END DO
                    ! X(L:LH, K:KH) = TRANSPOSE(X(K:KH,L:LH))

                    ! CALL DLA_ITRANSPOSE(K,KH,L,LH, X, LDX)
                END IF

                IF ( INFO1 .LT. 0 ) THEN
                    INFO = 100
                END IF

                IF ( SCAL .NE. ONE ) THEN
                    CALL DLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:M) = X(1:M, 1:M) * SCAL
                    SCALE = SCALE * SCAL
                    CALL DLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF

                ! Update January 2021
                IF ( K .GT. 1 ) THEN
                    CALL DGEMM('N', 'N', KE, LB, KB, -ONE, A(1,K), LDA, X(K,L), LDX, ONE, X(1,L), LDX)
                END IF

                K = K - 1
            END DO



            ! Update January 2021
            IF ( L .GT. 1 ) THEN
                CALL DSYR2K('Upper', 'NoTrans', LE, LB, -ONE, A(1,L), LDA, X(1, L), LDX, ONE, X(1,1), LDX)
!                !$omp parallel do private(KB)
!                DO K = 1, LE, 256
!                    KB = MIN(256, LE - K + 1)
!
!                    CALL DGEMM("N", "T", LE, KB, LB, -ONE, A(1,L), LDA, X(K,L), LDX, ONE, X(1,K), LDX)
!                    CALL DGEMM("N", "T", LE, KB, LB, -ONE, X(1,L), LDX, A(K,L), LDA, ONE, X(1,K), LDX)
!                END DO
!                !$omp end parallel do

                    ! CALL DGEMM("N", "T", LE, LE, LB, -ONE, A(1,L), LDA, X(1,L), LDX, ONE, X(1,1), LDX)
                    ! CALL DGEMM("N", "T", LE, LE, LB, -ONE, X(1,L), LDX, A(1,L), LDA, ONE, X(1,1), LDX)

            END IF

            L = L - 1
        END DO
    ELSE
        ! (T, N)
        K = 1
        KH = 1
        DO WHILE ( KH .LE. M )
            K = KH
            KH = MIN(M, K + MB - 1)
            KB = KH - K + 1
            IF ( KH .LT. M ) THEN
                IF ( A(KH+1, KH) .NE. ZERO ) THEN
                    KB = KB -1
                    KH = KH -1
                END IF
            END IF


            L = K
            LH = K
            DO WHILE ( LH .LE. M )
                L = LH
                LH = MIN(M, L + MB - 1)
                LB = LH - L + 1
                IF ( LH .LT. M ) THEN
                    IF ( ( A(LH+1,LH) .NE. ZERO ) ) THEN
                        LB = LB - 1
                        LH = LH - 1
                    END IF
                END IF

                IF ( L.EQ.K) THEN
                    ! Restore the symmetry in the diagonal block
                    DO I = K,KH-1
                        CALL DCOPY(KH-I,X(I,I+1),LDX,X(I+1,I),IONE)
                    END DO

                    IF ( ISOLVER .EQ. 6 ) THEN

                        CALL DLA_TRLYAP_RECURSIVE(TRANSA,  KB, A(K,K), LDA, X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE
                        IF ( KB .LE. 64 ) THEN
                            CALL DLA_TRLYAP_L2_OPT(TRANSA, KB, A(K,K), LDA,  X(K,K), LDX, SCAL, WORK, INFO1)
                        ELSE
                            CALL DLA_TRLYAP_L2(TRANSA,  KB, A(K,K), LDA, X(K,L), LDX, SCAL, WORK, INFO1)
                        END IF
                    END IF
                ELSE

                    SELECT CASE(ISOLVER)
                    CASE(1)
                        IF ( MB .LE. 32 ) THEN
                            CALL DLA_TRSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( MB .LE. 64) THEN
                            CALL DLA_TRSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( MB .LE. 96) THEN
                            CALL DLA_TRSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( MB .LE. 128) THEN
                            CALL DLA_TRSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE
                            CALL DLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        END IF
                    CASE(2)
                        IF ( KB .GT. 128 .OR. LB .GT. 128) THEN
                            CALL DLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE
                            CALL DLA_TRSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                & X(K,L), LDX, SCAL, WORK, INFO1)
                        END IF
                    CASE(3)
                        CALL DLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            & X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(4)
                        CALL DLA_TRSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            & X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(5)
                        CALL DLA_TRSYLV_L2_UNOPT(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            & X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(6)
                        CALL DLA_TRSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            & X(K,L), LDX, SCAL, WORK, INFO1)
                    END SELECT

                    CALL DLA_ITRANSPOSE(K,KH,L,LH,X,LDX)
                END IF

                IF ( SCAL .NE. ONE ) THEN
                    CALL DLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:M) = X(1:M, 1:M) * SCAL
                    SCALE = SCALE * SCAL
                    CALL DLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF


                ! in current row
                IF ( LH .LT. M ) THEN
                    CALL DGEMM('N','N', KB, M-LH, LB, -ONE, X(K,L), LDX, A(L,LH+1), LDA, ONE, X(K,LH+1), LDX)

                END IF


                LH = LH + 1
            END DO

            IF ( KH .LT. M ) THEN
                CALL DSYR2K("U", "T", M-KH, KB, -ONE, A(K,KH+1), LDA , X(K,KH+1) , LDX, ONE, X(KH+1,KH+1), LDX)
            END IF

            KH = KH + 1
        END DO

    END IF
END SUBROUTINE



