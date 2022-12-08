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


!> \brief DAG Scheduled Bartels-Stewart Algorithm for the standard Lyapunov Equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TRLYAP_DAG ( TRANS, M, A, LDA, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANS(1)
!           REAL SCALE
!           INTEGER M, LDA, LDX, INFO
!           REAL A(LDA, *), X(LDX, *)
!           REAL WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> SLA_TRLYAP_DAG solves a Lyapunov equation of the following forms
!>
!>    A  * X  +  X * A**T = SCALE * Y                              (1)
!>
!> or
!>
!>    A ** T * X  +  X * A = SCALE * Y                              (2)
!>
!> where A is a M-by-M quasi upper triangular matrix.
!> The right hand side Y and the solution X are M-by-N matrices.  Typically the matrix A
!> is generated by SGEES from LAPACK.
!> \endverbatim
!> \remark The algorithm is implemented using DAG Scheduling
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER(1)
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
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Dezember 2022
!> \ingroup sgltrlyap
!
SUBROUTINE SLA_TRLYAP_DAG ( TRANS, M, A, LDA, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANS(1)
    REAL SCALE
    INTEGER M, LDA, LDX, INFO
    REAL A(LDA, M), X(LDX, M)
    REAL WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, I, KOLD, LOLD
    INTEGER INFO1
    REAL SCAL
    INTEGER MB, ISOLVER, MAXNB
    INTEGER J,JH,JB, IB, IH

    LOGICAL BTRANS
    CHARACTER(1) TRANSA, TRANSB
    REAL ZERO , ONE, SGN, HALF
    PARAMETER (ZERO = 0.0, ONE = 1.0, SGN = 1.0, HALF = 0.5 )
    INTEGER IZERO, IONE
    PARAMETER(IZERO = 0, IONE = 1 )

    ! EXTERNAL FUNCTIONS
    EXTERNAL SCOPY
    EXTERNAL SGEMM
    EXTERNAL SLA_TRLYAP_L2
    EXTERNAL SLA_TRLYAP_RECURSIVE
    EXTERNAL SLA_TRSYLV_L2
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_128
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_32
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_64
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_96
    EXTERNAL SLA_TRSYLV_L2_REORDER
    EXTERNAL SLA_TRSYLV_L2_UNOPT
    EXTERNAL SLA_TRSYLV_RECURSIVE
    EXTERNAL ISLA_TRLYAP_SOLVE_BLOCK_T
    EXTERNAL ISLA_TRLYAP_SOLVE_BLOCK_N

    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME


    MB = TRLYAP_BLOCKSIZE(M)
    ISOLVER = TRLYAP_ISOLVER()
    MAXNB = MB


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
        CALL XERROR_HANDLER('SLA_TRLYAP_DAG', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1) THEN
        INFO = 1
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
        CALL SLA_TRLYAP_L2(TRANSA, KB, A, LDA,  X, LDX, SCAL, WORK, INFO1)

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

        !$omp parallel default(shared)
        !$omp master
        KOLD = 1
        LOLD = 1
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

                IF ( KH.EQ.M .AND. LH.EQ.M) THEN
                    !$omp task firstprivate(TRANSA,TRANSB,M,K,KB,KH,L,LB,LH,INFO1,SCAL) depend(out:X(K,L)) default(shared)
                    CALL ISLA_TRLYAP_SOLVE_BLOCK_N(TRANSA, TRANSB, M, K, KB, KH, L, LB, LH, A, LDA, X, LDX, SCAL, INFO1)
                    IF ( INFO1 .NE. 0 ) THEN
                        !$omp atomic
                        INFO = INFO + 1
                    END IF
                    IF ( SCAL .NE. ONE ) THEN
                        !$omp atomic
                        SCALE = SCALE * SCAL
                    END IF

                    !$omp end task
                ELSE IF ( KH.NE.M .AND. LH.EQ.M ) THEN
                    !$omp task firstprivate(TRANSA,TRANSB,M,K,KB,KH,L,LB,LH) default(shared) &
                    !$omp& depend(in:X(KOLD,L)) depend(out:X(K,L))
                    CALL ISLA_TRLYAP_SOLVE_BLOCK_N(TRANSA, TRANSB, M, K, KB, KH, L, LB, LH, A, LDA, X, LDX, SCAL, INFO1)
                    IF ( INFO1 .NE. 0 ) THEN
                        !$omp atomic
                        INFO = INFO + 1
                    END IF
                    IF ( SCAL .NE. ONE ) THEN
                        !$omp atomic
                        SCALE = SCALE * SCAL
                    END IF

                    !$omp end task
                ELSE IF ( KH.EQ.LH ) THEN
                    !$omp task firstprivate(TRANSA,TRANSB,M,K,KB,KH,L,LB,LH) default(shared) &
                    !$omp& depend(in:X(K,LOLD)) depend(out:X(K,L))
                    CALL ISLA_TRLYAP_SOLVE_BLOCK_N(TRANSA, TRANSB, M, K, KB, KH, L, LB, LH, A, LDA, X, LDX, SCAL, INFO1)
                    IF ( INFO1 .NE. 0 ) THEN
                        !$omp atomic
                        INFO = INFO + 1
                    END IF
                    IF ( SCAL .NE. ONE ) THEN
                        !$omp atomic
                        SCALE = SCALE * SCAL
                    END IF

                    !$omp end task
                ELSE
                    !$omp task firstprivate(TRANSA,TRANSB,M,K,KB,KH,L,LB,LH) default(shared) &
                    !$omp& depend(in:X(K,LOLD),X(KOLD,L)) depend(out:X(K,L))
                    CALL ISLA_TRLYAP_SOLVE_BLOCK_N(TRANSA, TRANSB, M, K, KB, KH, L, LB, LH, A, LDA, X, LDX, SCAL, INFO1)
                    IF ( INFO1 .NE. 0 ) THEN
                        !$omp atomic
                        INFO = INFO + 1
                    END IF
                    IF ( SCAL .NE. ONE ) THEN
                        !$omp atomic
                        SCALE = SCALE * SCAL
                    END IF

                    !$omp end task
                END IF


                ! Update January 2021
                IF ( K .GT. 1 ) THEN
                    J = K-1
                    DO WHILE (J .GT. 0)
                        JH = J
                        J = MAX(1, JH - MB + 1)
                        JB = JH - J + 1
                        IF ( J .GT. 1 ) THEN
                            IF (A(J,J-1) .NE. ZERO ) THEN
                                JB = JB -1
                                J = J + 1
                            END IF
                        END IF

                        !$omp task firstprivate(JB,KB,LB,J,K,L) default(shared) depend(in:X(K,L)) depend(inout: X(J,L))
                        CALL SGEMM('N', 'N', JB, LB, KB, -ONE, A(J,K), LDA, X(K,L), LDX, ONE, X(J,L), LDX)
                        !$omp end task
                        J = J - 1
                    END DO
                END IF

                KOLD = K
                K = K - 1
            END DO



            ! Update January 2021
            IF ( L .GT. 1 ) THEN
                J = L-1
                DO WHILE (J .GT. 0)
                    JH = J
                    J = MAX(1, JH - MB + 1)
                    JB = JH - J + 1
                    IF ( J .GT. 1 ) THEN
                        IF ( A(J,J-1) .NE. ZERO ) THEN
                            JB = JB -1
                            J = J + 1
                        END IF
                    END IF

                    I = JH
                    DO WHILE (I .GT. 0)
                        IH = I
                        I = MAX(1, IH - MB + 1)
                        IB = IH - I + 1
                        IF ( I .GT. 1 ) THEN
                            IF ( A(I,I-1) .NE. ZERO ) THEN
                                IB = IB - 1
                                I = I + 1
                            END IF
                        END IF

                        !$omp task firstprivate(I,IB,JB,LB,J,L) default(shared) depend(in:X(J,L),X(I,L)) depend(inout: X(I,J))
                        CALL SGEMM("N", "T", IB, JB, LB, -ONE, A(I,L), LDA, X(J,L), LDX, ONE, X(I,J), LDX)
                        !!$omp end task
                        !!$omp task firstprivate(I,IB,JB,LB,J,L) default(shared) depend(in:X(I,L)) depend(out: X(I,J))
                        CALL SGEMM("N", "T", IB, JB, LB, -ONE, X(I,L), LDX, A(J,L), LDA, ONE, X(I,J), LDX)
                        !$omp end task
                        I = I - 1
                    END DO
                    J = J - 1
                END DO
            END IF

            LOLD = L
            L = L - 1
        END DO
        !$omp end master
        !$omp taskwait
        !$omp end parallel
    ELSE
        ! (T, N)
        !$omp parallel default(shared)
        !$omp master

        KOLD = 1
        LOLD = 1
        K = 1
        KH = 1
        DO WHILE ( KH .LE. M )
            K = KH
            KH = MIN(M, K + MB - 1)
            KB = KH - K + 1
            IF ( KH .LT. M ) THEN
                IF ( A(KH+1, KH) .NE. ZERO ) THEN
                    KB = KB - 1
                    KH = KH - 1
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

                IF ( K.EQ.1 .AND. L.EQ. 1 ) THEN
                    !$omp task firstprivate(TRANSA,TRANSB,M,K,KB,KH,L,LB,LH,INFO1,SCAL) depend(out:X(K,L)) default(shared)
                    CALL ISLA_TRLYAP_SOLVE_BLOCK_T(TRANSA,TRANSB,M,K,KH,KB,L,LH,LB,A,LDA,X,LDX,SCAL,INFO1)
                    IF ( INFO1 .NE. 0 ) THEN
                        !$omp atomic
                        INFO = INFO + 1
                    END IF
                    IF ( SCAL .NE. ONE ) THEN
                        !$omp atomic
                        SCALE = SCALE * SCAL
                    END IF

                    !$omp end task
                ELSE IF (K.EQ.1 .AND. L.GT.1 ) THEN
                    !$omp task firstprivate(TRANSA,TRANSB,M,K,KB,KH,L,LB,LH) default(shared) &
                    !$omp& depend(in:X(K,LOLD)) depend(out:X(K,L))
                    CALL ISLA_TRLYAP_SOLVE_BLOCK_T(TRANSA,TRANSB,M,K,KH,KB,L,LH,LB,A,LDA,X,LDX,SCAL,INFO1)
                    IF ( INFO1 .NE. 0 ) THEN
                        !$omp atomic
                        INFO = INFO + 1
                    END IF
                    IF ( SCAL .NE. ONE ) THEN
                        !$omp atomic
                        SCALE = SCALE * SCAL
                    END IF

                    !$omp end task
                ELSE IF (K.EQ.L) THEN
                    !$omp task firstprivate(TRANSA,TRANSB,M,K,KB,KH,L,LB,LH) default(shared) &
                    !$omp& depend(in:X(KOLD,L)) depend(out:X(K,L))
                    CALL ISLA_TRLYAP_SOLVE_BLOCK_T(TRANSA,TRANSB,M,K,KH,KB,L,LH,LB,A,LDA,X,LDX,SCAL,INFO1)
                    IF ( INFO1 .NE. 0 ) THEN
                        !$omp atomic
                        INFO = INFO + 1
                    END IF
                    IF ( SCAL .NE. ONE ) THEN
                        !$omp atomic
                        SCALE = SCALE * SCAL
                    END IF

                    !$omp end task
                ELSE
                    !$omp task firstprivate(TRANSA,TRANSB,M,K,KB,KH,L,LB,LH) default(shared) &
                    !$omp& depend(in:X(K,LOLD),X(KOLD,L)) depend(out:X(K,L))
                    CALL ISLA_TRLYAP_SOLVE_BLOCK_T(TRANSA,TRANSB,M,K,KH,KB,L,LH,LB,A,LDA,X,LDX,SCAL,INFO1)
                    IF ( INFO1 .NE. 0 ) THEN
                        !$omp atomic
                        INFO = INFO + 1
                    END IF
                    IF ( SCAL .NE. ONE ) THEN
                        !$omp atomic
                        SCALE = SCALE * SCAL
                    END IF
                    !$omp end task
                END IF



                ! in current row
                IF ( LH .LT. M ) THEN
                    J = LH+1
                    JH = LH +1
                    DO WHILE ( JH .LE. M )
                        J = JH
                        JH = MIN(M, J + MB - 1)
                        JB = JH - J + 1
                        IF ( JH .LT. M ) THEN
                            IF ( ( A(JH+1,JH) .NE. ZERO ) ) THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF

                        !$omp task firstprivate(JB,KB,LB,J,K,L) default(shared) depend(in:X(K,L)) depend(inout: X(K,J))
                        CALL SGEMM('N','N', KB, JB, LB, -ONE, X(K,L), LDX, A(L,J), LDA, ONE, X(K,J), LDX)
                        !$omp end task
                        JH = JH + 1
                    END DO

                END IF


                LOLD = L
                LH = LH + 1
            END DO

            IF ( KH .LT. M ) THEN
                J = KH+1
                JH = KH +1
                DO WHILE ( JH .LE. M )
                    J = JH
                    JH = MIN(M, J + MB - 1)
                    JB = JH - J + 1
                    IF ( JH .LT. M ) THEN
                        IF ( ( A(JH+1,JH) .NE. ZERO ) ) THEN
                            JB = JB - 1
                            JH = JH - 1
                        END IF
                    END IF

                    I = KH+1
                    IH = KH +1
                    DO WHILE ( IH .LE. JH )
                        I = IH
                        IH= MIN(JH,I + MB - 1)
                        IB = IH- I + 1
                        IF ( IH .LT. M ) THEN
                            IF ( ( A(IH+1,IH) .NE. ZERO ) ) THEN
                                IB = IB - 1
                                IH = IH - 1
                            END IF
                        END IF

                        !$omp task firstprivate(I,IB,JB,KB,J,K) default(shared) depend(in:X(K,J),X(K,I)) depend(inout: X(I,J))
                        CALL SGEMM("T", "N", IB, JB, KB, -ONE, A(K,I), LDA, X(K,J), LDX, ONE, X(I, J), LDX)
                        CALL SGEMM("T", "N", IB, JB, KB, -ONE, X(K,I), LDX, A(K,J), LDA,  ONE, X(I, J), LDX)
                        !$omp end task
                        IH = IH +1
                    END DO
                    JH = JH + 1
                END DO

            END IF

            KOLD = K
            KH = KH + 1
        END DO

        !$omp end master
        !$omp taskwait
        !$omp end parallel

    END IF
    RETURN

END SUBROUTINE

SUBROUTINE ISLA_TRLYAP_SOLVE_BLOCK_N(TRANSA, TRANSB, M, K, KB, KH, L, LB, LH, A, LDA, X, LDX, SCAL, INFO)
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    IMPLICIT NONE
    CHARACTER TRANSA(1), TRANSB(1)
    INTEGER M, K, KB, KH, L, LB, LH, LDA, LDX, INFO
    REAL A(LDA, M), X(LDX, M), SCAL

    EXTERNAL SCOPY
    EXTERNAL SGEMM
    EXTERNAL SLA_TRLYAP_L2
    EXTERNAL SLA_TRLYAP_RECURSIVE
    EXTERNAL SLA_TRSYLV_L2
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_128
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_32
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_64
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_96
    EXTERNAL SLA_TRSYLV_L2_REORDER
    EXTERNAL SLA_TRSYLV_L2_UNOPT
    EXTERNAL SLA_TRSYLV_RECURSIVE


    INTEGER I, ISOLVER, IONE, MAXNB, II, JJ
    REAL WORK(10), SGN

    PARAMETER(IONE = 1, SGN = 1.0)

    ISOLVER = TRLYAP_ISOLVER()
    MAXNB = MAX(KB, LB)

    IF ( K .EQ. L ) THEN

        DO I = K,KH-1
            CALL SCOPY(KH-I,X(I,I+1),LDX,X(I+1,I),IONE)
        END DO

        IF ( ISOLVER .EQ. 6 ) THEN
            CALL SLA_TRLYAP_RECURSIVE(TRANSA, KB, A(K,K), LDA,  X(K,K), LDX, SCAL, WORK, INFO)

        ELSE
            CALL SLA_TRLYAP_L2(TRANSA, KB, A(K,K), LDA,  X(K,K), LDX, SCAL, WORK, INFO)
        END IF

    ELSE
        SELECT CASE(ISOLVER)
        CASE(1)
            IF ( MAXNB .LE. 32 ) THEN
                CALL SLA_TRSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            ELSE IF ( MAXNB .LE. 64) THEN
                CALL SLA_TRSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            ELSE IF ( MAXNB .LE. 96) THEN
                CALL SLA_TRSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            ELSE IF ( MAXNB .LE. 128) THEN
                CALL SLA_TRSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            ELSE
                CALL SLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            END IF
        CASE(2)
            IF ( KB .GT. 128 .OR. LB .GT. 128) THEN
                CALL SLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            ELSE
                CALL SLA_TRSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            END IF
        CASE(3)
            CALL SLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                & X(K,L), LDX, SCAL, WORK, INFO)
        CASE(4)
            CALL SLA_TRSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                & X(K,L), LDX, SCAL, WORK, INFO)
        CASE(5)
            CALL SLA_TRSYLV_L2_UNOPT(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                & X(K,L), LDX, SCAL, WORK, INFO)
        CASE(6)
            CALL SLA_TRSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                & X(K,L), LDX, SCAL, WORK, INFO)
        END SELECT
    END IF

    IF ( K .NE. L ) THEN
        ! Transpose causes a temporary here
        ! X(L:LH,K:KH) = TRANSPOSE(X(K:KH,L:LH))
        DO II = K, KH
            DO JJ = L, LH
                X(JJ, II) = X(II,JJ)
            END DO
        END DO
    END IF

END SUBROUTINE

SUBROUTINE ISLA_TRLYAP_SOLVE_BLOCK_T(TRANSA,TRANSB,M,K,KH,KB,L,LH,LB,A,LDA,X,LDX,SCAL,INFO)
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    IMPLICIT NONE
    CHARACTER TRANSA(1), TRANSB(1)
    INTEGER M, K, KB, KH, L, LB, LH, LDA, LDX, INFO
    REAL A(LDA, M), X(LDX, M), SCAL

    EXTERNAL SCOPY
    EXTERNAL SGEMM
    EXTERNAL SLA_TRLYAP_L2
    EXTERNAL SLA_TRLYAP_RECURSIVE
    EXTERNAL SLA_TRSYLV_L2
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_128
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_32
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_64
    EXTERNAL SLA_TRSYLV_L2_LOCAL_COPY_96
    EXTERNAL SLA_TRSYLV_L2_REORDER
    EXTERNAL SLA_TRSYLV_L2_UNOPT
    EXTERNAL SLA_TRSYLV_RECURSIVE


    INTEGER I, ISOLVER, IONE, MAXNB, II,JJ
    REAL WORK(10), SGN

    PARAMETER(IONE = 1, SGN = 1.0)

    ISOLVER = TRLYAP_ISOLVER()
    MAXNB = MAX(KB, LB)


    IF ( L.EQ.K) THEN
        ! Restore the symmetry in the diagonal block
        DO I = K,KH-1
            CALL SCOPY(KH-I,X(I,I+1),LDX,X(I+1,I),IONE)
        END DO
        IF ( ISOLVER .EQ. 6 ) THEN
            CALL SLA_TRLYAP_RECURSIVE(TRANSA,  KB, A(K,K), LDA, X(K,L), LDX, SCAL, WORK, INFO)
        ELSE
            CALL SLA_TRLYAP_L2(TRANSA,  KB, A(K,K), LDA, X(K,L), LDX, SCAL, WORK, INFO)
        END IF
    ELSE

        SELECT CASE(ISOLVER)
        CASE(1)
            IF ( MAXNB .LE. 32 ) THEN
                CALL SLA_TRSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            ELSE IF ( MAXNB .LE. 64) THEN
                CALL SLA_TRSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            ELSE IF ( MAXNB .LE. 96) THEN
                CALL SLA_TRSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            ELSE IF ( MAXNB .LE. 128) THEN
                CALL SLA_TRSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            ELSE
                CALL SLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            END IF
        CASE(2)
            IF ( KB .GT. 128 .OR. LB .GT. 128) THEN
                CALL SLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            ELSE
                CALL SLA_TRSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                    & X(K,L), LDX, SCAL, WORK, INFO)
            END IF
        CASE(3)
            CALL SLA_TRSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                & X(K,L), LDX, SCAL, WORK, INFO)
        CASE(4)
            CALL SLA_TRSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                & X(K,L), LDX, SCAL, WORK, INFO)
        CASE(5)
            CALL SLA_TRSYLV_L2_UNOPT(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                & X(K,L), LDX, SCAL, WORK, INFO)
        CASE(6)
            CALL SLA_TRSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, A(L,L), LDA, &
                & X(K,L), LDX, SCAL, WORK, INFO)
        END SELECT

    END IF

    ! Copy Block to its transposed place
    IF ( K.NE.L) THEN
        ! Transpose causes a temporary here
        ! X(L:LH,K:KH) = TRANSPOSE(X(K:KH,L:LH))

        DO II = K, KH
            DO JJ = L, LH
                X(JJ, II) = X(II,JJ)
            END DO
        END DO
    END IF
END SUBROUTINE




