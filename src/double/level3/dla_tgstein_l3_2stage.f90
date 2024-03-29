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


!> \brief Level-3 Bartels-Stewart Algorithm for the generalized Stein equation (Two Stage)
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TGSTEIN_L3_2S ( TRANS, M,  A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
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
!> DLA_TGSTEIN_L3_2S solves a generalized Lyapunov  equation of the following forms
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
!> \attention The algorithm is implemented using BLAS level 3 operations with two level blocking.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A, B,  and C:
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
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup dbltglyap
!
SUBROUTINE DLA_TGSTEIN_L3_2S ( TRANS, M, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    USE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
    USE MEPACK_OPTIONS_ISOLVER
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANS(1)
    DOUBLE PRECISION SCALE
    INTEGER M, LDA, LDB, LDX, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *) , X(LDX, *)
    DOUBLE PRECISION WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER I, K, L, KB, LB, KH, LE, LH, KE
    INTEGER INFO1
    DOUBLE PRECISION SCAL
    INTEGER NB, ISOLVER
    INTEGER DAG_INFO, DAG_INFO2

    LOGICAL BTRANS
    DOUBLE PRECISION ZERO , ONE, HALF
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0)
    INTEGER IZERO, IONE, BIGNB
    PARAMETER(IZERO = 0, IONE = 1 )
    CHARACTER TRANSA(1), TRANSB(1)
    DOUBLE PRECISION SGN
    PARAMETER(SGN = -ONE)

    ! EXTERNAL FUNCTIONS
    EXTERNAL DAXPY
    EXTERNAL DCOPY
    EXTERNAL DGEMM
    EXTERNAL DLA_ITRANSPOSE
    EXTERNAL DLA_TGSTEIN_DAG
    EXTERNAL DLA_TGSYLV_DAG
    EXTERNAL DLA_TGSYLV_L2
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_128
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_32
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_64
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_96
    EXTERNAL DLA_TGSYLV_L2_REORDER
    EXTERNAL DLA_TGSYLV_L2_UNOPT
    EXTERNAL DLA_TGSYLV_RECURSIVE
    EXTERNAL DSCAL
    EXTERNAL DSYMM
    EXTERNAL DSYR2K
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME



    NB = TGSTEIN_BLOCKSIZE(M)
    ISOLVER = TGSTEIN_ISOLVER()
    BIGNB = TGSTEIN_BLOCKSIZE_2STAGE(M)


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
    ELSE IF ( LDB .LT. MAX(1, M)) THEN
        INFO = -6
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -8
    ELSE IF ( NB .LT. 2 ) THEN
        INFO = -30
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32
    ELSE IF ( BIGNB .LT. NB ) THEN
        INFO = -33
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('DLA_TGSTEIN_L3_2S', -INFO)
        RETURN
    END IF

    IF ( BTRANS ) THEN
        TRANSA = 'N'
        TRANSB = 'T'
    ELSE
        TRANSA = 'T'
        TRANSB = 'N'
    END IF

    IF ( ININFO .EQ. -1 ) THEN
        ! MEMORY DAG
        DAG_INFO = -1
        DAG_INFO2 = -1
        CALL DLA_TGSTEIN_DAG ( TRANS, BIGNB, A, BIGNB, B, BIGNB,  X, BIGNB, SCALE, WORK, DAG_INFO)
        CALL DLA_TGSYLV_DAG(TRANSA, TRANSB, ONE, BIGNB, BIGNB,  A, BIGNB, A, BIGNB, &
                        & B, BIGNB, B, BIGNB,  X, BIGNB, SCAL, WORK, DAG_INFO2)
        INFO = MAX( DAG_INFO, DAG_INFO2, BIGNB**2, M*BIGNB, 4*NB)

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

    IF ( M .LE. NB ) THEN
        KB = M
        LB = M
        SELECT CASE(ISOLVER)
        CASE(1)
            IF ( NB .LE. 32 ) THEN
                CALL DLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                    & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( NB .LE. 64) THEN
                CALL DLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                    & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( NB .LE. 96) THEN
                CALL DLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                    & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( NB .LE. 128) THEN
                CALL DLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                    & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
            ELSE
                CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                    & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
            END IF
        CASE(2)
            IF ( M .GT. 128 ) THEN
                CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                    & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
            ELSE
                CALL DLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                    & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
            END IF
        CASE(3)
            CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
        CASE(4)
            CALL DLA_TGSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
        CASE(5)
            CALL DLA_TGSYLV_L2_UNOPT(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
        CASE(6)
            CALL DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, A, LDA, &
                & B, LDB, B, LDB,  X, LDX, SCAL, WORK, INFO1)
        END SELECT

        ! Symmetries
        DO K = 1,M-1
            CALL DSCAL(M-K,HALF,X(K+1,K),IONE)
            CALL DAXPY(M-K,HALF,X(K,K+1),LDX,X(K+1,K),IONE)
            CALL DCOPY(M-K,X(K+1,K),IONE,X(K,K+1),LDX)
        END DO


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
            IF (L .EQ. M .AND. MOD(M,BIGNB) .NE. IZERO) THEN
                LH = L
                L = (L/BIGNB)*BIGNB + 1
            ELSE
                LH = L
                L = MAX(1, LH - BIGNB + 1)
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
                IF (K .EQ. M .AND. MOD(M,BIGNB) .NE. IZERO) THEN
                    KH = K
                    K = (M/BIGNB)*BIGNB + 1
                ELSE
                    KH = K
                    K = MAX(1, KH - BIGNB + 1)
                END IF

                KB = KH - K + 1
                IF ( K .GT. 1 ) THEN
                    IF( A(K,K-1) .NE. ZERO ) THEN
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

                IF ( L.EQ.K) THEN
                    DO I = K,KH-1
                        CALL DCOPY(KH-I,X(I,I+1),LDX,X(I+1,I),IONE)
                    END DO
                    CALL DLA_TGSTEIN_DAG ( TRANS, KB, A(K,K) , LDA, B(K,K), LDB, X(K,K), LDX, SCALE, WORK, INFO)
                    IF ( INFO1 .LT. 0 ) THEN
                        INFO = 100
                    END IF

                ELSE
                    CALL DLA_TGSYLV_DAG(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                        & B(K,K), LDB, B(L,L), LDB,  X(K,L), LDX, SCAL, WORK, INFO1)

                    IF ( INFO1 .LT. 0 ) THEN
                        INFO = 100
                    END IF

                    CALL DLA_ITRANSPOSE(K,KH,L,LH, X, LDX)
                END IF


                ! Update January 2021
                IF ( K .GT. 1 ) THEN
                    CALL DGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, A(L,L), LDA, ZERO, WORK, KB)
                    CALL DGEMM('N', 'N', KE, LB, KB, -ONE, A(1,K), LDA, WORK, KB, ONE, X(1,L), LDX)

                    CALL DGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK, KB)
                    CALL DGEMM('N', 'N', KE, LB, KB, ONE, B(1,K), LDB, WORK, KB, ONE, X(1,L), LDX)
                END IF

                K = K - 1
            END DO


            ! Update January 2021
            IF ( LE .GT. 0 ) THEN

                CALL DGEMM("N", "N", LE, LB, LE, ONE, A(1,1), LDA, X(1,L), LDX, ZERO, WORK, LE)
                CALL DSYR2K("U", "N", LE, LB, -ONE, WORK, LE, A(1,L), LDA, ONE, X(1,1), LDX)

                CALL DGEMM("N", "N", LE, LB, LE, ONE, B(1,1), LDB, X(1,L), LDX, ZERO, WORK, LE)
                CALL DSYR2K("U", "N", LE, LB, ONE, WORK, LE, B(1,L), LDB, ONE, X(1,1), LDX)

                CALL DSYMM("Right", "Upper", LE, LB, ONE, X(L,L), LDX, A(1,L), LDA, ZERO, WORK, LE)
                CALL DGEMM("N", "T", LE, LE, LB, -ONE, WORK, LE, A(1,L), LDA, ONE, X(1,1), LDX)
                CALL DSYMM("Right", "Upper", LE, LB, ONE, X(L,L), LDX, B(1,L), LDB, ZERO, WORK, LE)
                CALL DGEMM("N", "T", LE, LE, LB, ONE, WORK, LE, B(1,L), LDB, ONE, X(1,1), LDX)

           END IF

            L = L - 1
        END DO

    ELSE
        ! (T, N)
        L = 1
        LH = 1
        DO WHILE ( LH .LE. M )
            L = LH
            LH = MIN(M, L + BIGNB - 1)
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
                KH = MIN(M, K + BIGNB - 1)
                KB = KH - K + 1
                IF ( KH .LT. M ) THEN
                    IF ( A(KH+1, KH) .NE. ZERO ) THEN
                    KB = KB - 1
                    KH = KH - 1
                END IF
                END IF

                IF ( L.EQ.K) THEN
                    DO I = K,KH-1
                        CALL DCOPY(KH-I,X(I+1,I),IONE,X(I,I+1),LDX)
                    END DO
                    CALL DLA_TGSTEIN_DAG ( TRANS, KB, A(K,K) , LDA, B(K,K), LDB, X(K,K), LDX, SCALE, WORK, INFO)
                    IF ( INFO1 .LT. 0 ) THEN
                        INFO = 100
                    END IF

                ELSE
                    CALL DLA_TGSYLV_DAG(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                        & B(K,K), LDB, B(L,L), LDB,  X(K,L), LDX, SCAL, WORK, INFO1)

                    IF ( INFO1 .LT. 0 ) THEN
                        INFO = 100
                    END IF

                    CALL DLA_ITRANSPOSE(K,KH,L,LH, X, LDX)
                END IF

                ! in current Column
                IF ( KH  .LT. M ) THEN
                    CALL DGEMM("N", "N", KB, LB, LB, ONE, X(K,L), LDX, A(L,L), LDA, ZERO, WORK, KB)
                    CALL DGEMM("T", "N", M-KH, LB, KB, -ONE, A(K,KH+1), LDA, WORK, KB, ONE, X(KH+1, L), LDX)

                    CALL DGEMM("N", "N", KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK, KB)
                    CALL DGEMM("T", "N", M-KH, LB, KB, ONE, B(K,KH+1), LDB, WORK, KB, ONE, X(KH+1, L), LDX)
                END IF

                KH = KH + 1
            END DO

            ! Update January 2021
            IF ( LH .LT. M ) THEN
                CALL DGEMM("N", "N", LB, M-LH, M-LH, ONE, X(L,LH+1), LDX, A(LH+1,LH+1), LDA, ZERO, WORK, LB)
                CALL DSYR2K("L", "T", M-LH, LB, -ONE, A(L,LH+1), LDA, WORK, LB, ONE, X(LH+1,LH+1), LDX)

                CALL DGEMM("N", "N", LB, M-LH, M-LH, ONE, X(L,LH+1), LDX, B(LH+1,LH+1), LDB, ZERO, WORK, LB)
                CALL DSYR2K("L", "T", M-LH, LB, ONE, B(L,LH+1), LDB, WORK, LB, ONE, X(LH+1,LH+1), LDX)

                CALL DSYMM("Left", "Lower", LB, M-LH, ONE, X(L,L), LDX, A(L,LH+1), LDA, ZERO, WORK, LB)
                CALL DGEMM("T", "N", M-LH, M-LH, LB, -ONE, A(L,LH+1), LDA, WORK, LB, ONE, X(LH+1,LH+1), LDX)

                CALL DSYMM("Left", "Lower", LB, M-LH, ONE, X(L,L), LDX, B(L,LH+1), LDB, ZERO, WORK, LB)
                CALL DGEMM("T", "N", M-LH, M-LH, LB,  ONE, B(L,LH+1), LDB, WORK, LB, ONE, X(LH+1,LH+1), LDX)


            END IF

            LH = LH + 1
        END DO
    END IF
END SUBROUTINE



