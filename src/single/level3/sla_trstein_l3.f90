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


!> \brief Level-3 Bartels-Stewart Algorithm for the Stein equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TRSTEIN_L3 ( TRANS, M,  A, LDA, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANS
!           REAL SCALE
!           INTEGER M, N, LDA, LDX, INFO
!           REAL A(LDA, *), X(LDX, *)
!           REAL WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TRSTEIN_L3 solves a standard Stein  equation of the following forms
!>
!>    A * X * A^T - X  = SCALE * Y                                              (2)
!>
!> or
!>
!>    A^T * X * A - X  =  SCALE * Y                                             (1)
!>
!> where A is a M-by-M quasi upper triangular matrix.
!> and X and Y are symmetric  M-by-M matrices.
!> Typically the matrix A is created by SGEES from LAPACK.
!> \endverbatim
!> \attention The algorithm is implemented using BLAS level 3 operations.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER
!>          Specifies the form of the system of equations with respect to A and C:
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
!>          The leading dimension of the array X.  LDX >= max(1,M).
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
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup sgltrlyap
!
SUBROUTINE SLA_TRSTEIN_L3 ( TRANS, M, A, LDA, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE MEPACK_OPTIONS_ISOLVER
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANS(1)
    REAL SCALE
    INTEGER M, LDA, LDX, INFO
    REAL A(LDA, *), X(LDX, *)
    REAL WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE
    INTEGER INFO1
    REAL SCAL
    INTEGER NB, ISOLVER

    LOGICAL BTRANS
    REAL ZERO , ONE, HALF
    PARAMETER (ZERO = 0.0, ONE = 1.0, HALF = 0.5)
    INTEGER IZERO, IONE
    PARAMETER(IZERO = 0, IONE = 1 )
    CHARACTER TRANSA(1), TRANSB(1)
    REAL SGN
    PARAMETER(SGN = -ONE)

    ! EXTERNAL FUNCTIONS
    EXTERNAL SGEMM
    EXTERNAL SLACPY
    EXTERNAL SLA_ITRANSPOSE
    EXTERNAL SLA_TRSTEIN_L2
    EXTERNAL SLA_TRSTEIN_RECURSIVE
    EXTERNAL SLA_TRSYLV2_L2
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_128
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_32
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_64
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_96
    EXTERNAL SLA_TRSYLV2_L2_REORDER
    EXTERNAL SLA_TRSYLV2_L2_UNOPT
    EXTERNAL SLA_TRSYLV2_RECURSIVE
    EXTERNAL SSYMM
    EXTERNAL SSYR2K
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME

    NB = TRSTEIN_BLOCKSIZE(M)
    ISOLVER = TRSTEIN_ISOLVER()


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
    ELSE IF ( NB .LT. 2 ) THEN
        INFO = -30
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_TRSTEIN_L3', -INFO)
        RETURN
    END IF


    IF ( ININFO .EQ. -1 ) THEN
        INFO = MAX(4*NB, NB*NB, M*NB)
        RETURN
    END IF



    ! Quick Return
    SCALE = ONE
    SCAL  = ONE
    INFO = 0
    INFO1 = 0
    IF ( M .EQ. 0 ) THEN
        RETURN
    END IF

    IF ( BTRANS ) THEN
        TRANSA = 'N'
        TRANSB = 'T'
    ELSE
        TRANSA = 'T'
        TRANSB = 'N'
    END IF

    IF ( M .LE. NB ) THEN
        KB = M
        LB = M
        IF ( ISOLVER .LT. 6 ) THEN
            CALL SLA_TRSTEIN_L2(TRANS,M,A,LDA,X,LDX, SCAL, WORK, INFO1)
        ELSE IF (ISOLVER .EQ. 6 ) THEN
            CALL SLA_TRSTEIN_RECURSIVE(TRANS,M,A,LDA,X,LDX, SCAL, WORK, INFO1)
        END IF

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

    IF ( BTRANS ) THEN
        ! (N,T)

        L = M
        DO WHILE ( L .GT. 0 )
            IF (L .EQ. M .AND. MOD(M,NB) .NE. IZERO) THEN
                LH = L
                L = (L/NB)*NB + 1
            ELSE
                LH = L
                L = MAX(1, LH - NB + 1)
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
                IF (K .EQ. M .AND. MOD(M,NB) .NE. IZERO) THEN
                    KH = K
                    K = (M/NB)*NB + 1
                ELSE
                    KH = K
                    K = MAX(1, KH - NB + 1)
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


                IF ( L.EQ. K ) THEN
                    ! DO I = K,KH-1
                    !     CALL SCOPY(KH-I,X(I,I+1),LDX,X(I+1,I),IONE)
                    ! END DO
                    IF (ISOLVER .EQ. 6 ) THEN
                        CALL SLA_TRSTEIN_RECURSIVE(TRANS,KB,A(K,K),LDA,X(K,K),LDX, SCAL, WORK, INFO1)
                    ELSE
                        CALL SLA_TRSTEIN_L2(TRANS,KB,A(K,K),LDA,X(K,K),LDX, SCAL, WORK, INFO1)
                    END IF
                    IF ( INFO1 .LT. 0 ) THEN
                        INFO = 100
                    END IF

                ELSE
                    SELECT CASE(ISOLVER)
                    CASE(1)
                        IF ( NB .LE. 32 ) THEN
                            CALL SLA_TRSYLV2_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( NB .LE. 64) THEN
                            CALL SLA_TRSYLV2_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( NB .LE. 96) THEN
                            CALL SLA_TRSYLV2_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( NB .LE. 128) THEN
                            CALL SLA_TRSYLV2_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE
                            CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        END IF
                    CASE(2)
                        IF ( KB .GT. 128 .OR. LB .GT. 128) THEN
                            CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE
                            CALL SLA_TRSYLV2_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        END IF

                    CASE(3)
                        CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            &   X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(4)
                        CALL SLA_TRSYLV2_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            &   X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(5)
                        CALL SLA_TRSYLV2_L2_UNOPT(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            &   X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(6)
                        CALL SLA_TRSYLV2_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            &   X(K,L), LDX, SCAL, WORK, INFO1)
                    END SELECT
                    IF ( INFO1 .LT. 0 ) THEN
                        INFO = -100
                    END IF

                    CALL SLA_ITRANSPOSE(K,KH,L,LH,X,LDX)

                END IF

                IF ( SCAL .NE. ONE ) THEN
                    CALL SLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:M) = X(1:M, 1:M) * SCAL
                    SCALE = SCALE * SCAL
                    CALL SLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF

                ! Update January 2021
                IF ( K .GT. 1 ) THEN
                    CALL SGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, A(L,L), LDA, ZERO, WORK, KB)
                    CALL SGEMM('N', 'N', KE, LB, KB, -ONE, A(1,K), LDA, WORK, KB, ONE, X(1,L), LDX)
                END IF

                K = K - 1
            END DO


            ! Update January 2021
            IF ( LE .GT. 0 ) THEN

                CALL SGEMM("N", "N", LE, LB, LE, ONE, A(1,1), LDA, X(1,L), LDX, ZERO, WORK, LE)
                CALL SSYR2K("U", "N", LE, LB, -ONE, WORK, LE, A(1,L), LDA, ONE, X(1,1), LDX)


                CALL SSYMM("Right", "Upper", LE, LB, ONE, X(L,L), LDX, A(1,L), LDA, ZERO, WORK, LE)
                CALL SGEMM("N", "T", LE, LE, LB, -ONE, WORK, LE, A(1,L), LDA, ONE, X(1,1), LDX)
            END IF

            L = L - 1
        END DO

    ELSE
        ! (T, N)
        L = 1
        LH = 1
        DO WHILE ( LH .LE. M )
            L = LH
            LH = MIN(M, L + NB - 1)
            LB = LH - L + 1
            IF ( LH .LT. M ) THEN
                IF ( A(LH+1,LH) .NE. ZERO ) THEN
                    LB = LB - 1
                    LH = LH - 1
                END IF
            END IF

            K = L
            KH = L
            DO WHILE ( KH .LE. M )
                K = KH
                KH = MIN(M, K + NB - 1)
                KB = KH - K + 1
                IF ( KH .LT. M ) THEN
                    IF ( A(KH+1, KH) .NE. ZERO ) THEN
                        KB = KB - 1
                        KH = KH - 1
                    END IF
                END IF

                IF ( L.EQ.K) THEN
                    IF (ISOLVER .EQ. 6 ) THEN
                        CALL SLA_TRSTEIN_RECURSIVE(TRANS,KB,A(K,K),LDA,X(K,K),LDX, SCAL, WORK, INFO1)
                    ELSE
                        CALL SLA_TRSTEIN_L2(TRANS,KB,A(K,K),LDA,X(K,K),LDX, SCAL, WORK, INFO1)
                    END IF
                    IF ( INFO1 .LT. 0 ) THEN
                        INFO = 100
                    END IF

                ELSE

                    SELECT CASE(ISOLVER)
                    CASE(1)
                        IF ( NB .LE. 32 ) THEN
                            CALL SLA_TRSYLV2_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( NB .LE. 64) THEN
                            CALL SLA_TRSYLV2_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( NB .LE. 96) THEN
                            CALL SLA_TRSYLV2_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE IF ( NB .LE. 128) THEN
                            CALL SLA_TRSYLV2_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE
                            CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        END IF
                    CASE(2)
                        IF ( KB .GT. 128 .OR. LB .GT. 128) THEN
                            CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        ELSE
                            CALL SLA_TRSYLV2_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                                &   X(K,L), LDX, SCAL, WORK, INFO1)
                        END IF

                    CASE(3)
                        CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            &   X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(4)
                        CALL SLA_TRSYLV2_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            &   X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(5)
                        CALL SLA_TRSYLV2_L2_UNOPT(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            &   X(K,L), LDX, SCAL, WORK, INFO1)
                    CASE(6)
                        CALL SLA_TRSYLV2_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                            &   X(K,L), LDX, SCAL, WORK, INFO1)
                    END SELECT

                    IF ( INFO1 .LT. 0 ) THEN
                        INFO = 100
                    END IF

                    CALL SLA_ITRANSPOSE(K,KH,L,LH,X, LDX)
                END IF

                IF ( SCAL .NE. ONE ) THEN
                    CALL SLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:M) = X(1:M, 1:M) * SCAL
                    SCALE = SCALE * SCAL
                    CALL SLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF


                ! in current Column
                IF ( KH  .LT. M ) THEN
                    CALL SGEMM("N", "N", KB, LB, LB, ONE, X(K,L), LDX, A(L,L), LDA, ZERO, WORK, KB)
                    CALL SGEMM("T", "N", M-KH, LB, KB, -ONE, A(K,KH+1), LDA, WORK, KB, ONE, X(KH+1, L), LDX)

                END IF

                KH = KH + 1
            END DO


            ! Update January 2021
            IF ( LH .LT. M ) THEN
                CALL SGEMM("N", "N", LB, M-LH, M-LH, ONE, X(L,LH+1), LDX, A(LH+1,LH+1), LDA, ZERO, WORK, NB)
                CALL SSYR2K("L", "T", M-LH, LB, -ONE, A(L,LH+1), LDA, WORK, NB, ONE, X(LH+1,LH+1), LDX)


                CALL SSYMM("Left", "Lower", LB, M-LH, ONE, X(L,L), LDX, A(L,LH+1), LDA, ZERO, WORK, LB)
                CALL SGEMM("T", "N", M-LH, M-LH, LB, -ONE, A(L,LH+1), LDA, WORK, LB, ONE, X(LH+1,LH+1), LDX)

            END IF

            LH = LH + 1
        END DO

    END IF
END SUBROUTINE



