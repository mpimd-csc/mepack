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


!> \brief Level-3 Bartels-Stewart Algorithm for the generalized Lyapunov equation with DAG scheduling
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TGLYAP_DAG ( TRANS, M,  A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANS
!           REAL SCALE
!           INTEGER M, N, LDA, LDB, LDX, INFO
!           REAL A(LDA, *), B(LDB, *) , X(LDX, *)
!           REAL WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TGLYAP_DAG solves a generalized Lyapunov  equation of the following forms
!>
!>    A * X * B^T + A * X * B^T = SCALE * Y                                              (2)
!>
!> or
!>
!>    A^T * X * B + A^T * X * B =  SCALE * Y                                             (1)
!>
!> where A is a M-by-M quasi upper triangular matrix, B is a M-by-M upper triangular,
!> and X and Y are symmetric  M-by-M matrices.
!> Typically the matrix pencil (A,B) is created by SGGES from LAPACK.
!> \endverbatim
!> \attention The algorithm is implemented using BLAS level 3 operations and DAG .
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER
!>          Specifies the form of the system of equations with respect to A and B:
!>          == 'N':  Equation (1) is solved.
!>          == 'T':  Equation (2) is solved.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and B.  M >= 0.
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
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB,M)
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
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date June 2023
!> \ingroup sgltglyap
!
SUBROUTINE SLA_TGLYAP_DAG ( TRANS, M, A, LDA, B, LDB,  X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE OMP_LIB
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANS(1)
    REAL SCALE
    INTEGER M, LDA, LDX, LDB, INFO
    REAL A(LDA, *), X(LDX, *), B(LDB,*)
    REAL WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE
    INTEGER J, JB, JH, KOLD, LOLD
    INTEGER INFO1
    REAL SCAL
    INTEGER MB, ISOLVER, MAXNB, WORKOFFSET
    INTEGER NTH

    LOGICAL BTRANS
    REAL ZERO , ONE, HALF
    PARAMETER (ZERO = 0.0, ONE = 1.0, HALF = 0.5)
    INTEGER IZERO, IONE
    PARAMETER(IZERO = 0, IONE = 1 )
    CHARACTER TRANSA, TRANSB
    REAL SGN
    PARAMETER(SGN = ONE)

    ! EXTERNAL FUNCTIONS
    EXTERNAL SGEMM
    EXTERNAL SLA_TGLYAP_L2
    EXTERNAL SLA_TGSYLV_L2
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_128
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_32
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_64
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_96
    EXTERNAL SLA_TGSYLV_L2_REORDER
    EXTERNAL SLA_TGSYLV_L2_COLWISE
    EXTERNAL SLA_TGSYLV_RECURSIVE
    EXTERNAL ISLA_TGLYAP_SOLVE_BLOCK_N
    EXTERNAL ISLA_TGLYAP_SOLVE_BLOCK_T

    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME


    MB = TGLYAP_BLOCKSIZE(M)
    ISOLVER = TGLYAP_ISOLVER()
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
    ELSE IF ( LDB .LT. MAX(1, M)) THEN
        INFO = -6
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -8
    ELSE IF ( MB .LT. 2 ) THEN
        INFO = -30
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_TGLYAP_DAG', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1) THEN
        NTH = OMP_GET_MAX_THREADS()
        INFO = 4*M*M + NTH*MAX(MB*MB, 4*MB)
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

    IF ( M .LE. MB ) THEN
        KB = M
        LB = M

        CALL SLA_TGLYAP_L2 (TRANSA, KB, A, LDA, B, LDB, X, LDX, SCAL, WORK, INFO1)

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
                IF (A(L,L-1) .NE. ZERO ) THEN
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

                    !$omp task firstprivate(TRANSA,TRANSB,M,K,KB,KH,L,LB,LH,INFO1,SCAL) depend(inout:X(K,L)) default(shared)
                    CALL ISLA_TGLYAP_SOLVE_BLOCK_N(TRANSA, TRANSB, M, MB,K, KB, KH, L, LB, LH, A, LDA, B, LDB, X, LDX, WORK, &
                        & SCAL, INFO)

                    !$omp end task


                ! Prepare Update January 2021
                IF ( K .GT. 1 ) THEN
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                    WORKOFFSET = (L-1)*M+(K-1)*LB+1+2*M*M
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, A(L,L), LDA, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                END IF

                IF ( L.NE.K .AND. L.GT.1) THEN
                    WORKOFFSET = (K-1)*M+(L-1)*KB+1
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM('T', 'T', LB, KB, KB, ONE, X(K,L), LDX, B(K,K), LDB, ZERO, WORK(WORKOFFSET), LB)
                    !$omp end task

                    WORKOFFSET = (K-1)*M+(L-1)*KB+1+2*M*M
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM('T', 'T', LB, KB, KB, ONE, X(K,L), LDX, A(K,K), LDA, ZERO, WORK(WORKOFFSET), LB)
                    !$omp end task

                END IF



                IF ( K .LT. L) THEN
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1+M*M
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM("N", "N", KB, LB, M-K+1, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                    WORKOFFSET = (L-1)*M+(K-1)*LB+1+3*M*M
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM("N", "N", KB, LB, M-K+1, ONE, B(K,K), LDB, X(K,L), LDX, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                END IF

                !
                ! ! Update January 2021
                IF ( K .GT. 1 ) THEN
                    J = K-1
                    DO WHILE (J .GT. 0)
                        JH = J
                        J = MAX(1, JH - MB + 1)

                        JB = JH - J + 1
                        IF ( J .GT. 1 ) THEN
                            IF ( A(J,J-1) .NE. ZERO ) THEN
                                JB = JB - 1
                                J = J + 1
                            END IF
                        END IF

                        WORKOFFSET = (L-1)*M+(K-1)*LB+1
                        !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,L))
                        CALL SGEMM('N', 'N', JB, LB, KB, -ONE, A(J,K), LDA, WORK(WORKOFFSET), KB, ONE, X(J,L), LDX)
                        !$omp end task


                        WORKOFFSET = (L-1)*M+(K-1)*LB+1 + 2*M*M
                        !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,L))
                        CALL SGEMM('N', 'N', JB, LB, KB, -ONE, B(J,K), LDB, WORK(WORKOFFSET), KB, ONE, X(J,L), LDX)
                        !$omp end task

                        J = J - 1
                    END DO
                END IF

                IF ( L .GT. 1 )  THEN
                    IF ( L.NE.K) THEN
                        J = KH
                        DO WHILE (J .GT. 0)
                            JH = J
                            J = MAX(1, JH - MB + 1)

                            JB = JH - J + 1
                            IF ( J .GT. 1 ) THEN
                                IF ( A(J,J-1) .NE. ZERO ) THEN
                                    JB = JB - 1
                                    J = J + 1
                                END IF
                            END IF


                            WORKOFFSET = (K-1)*M+(L-1)*KB+1
                            !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                            !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,K))
                            CALL SGEMM('N', 'N', JB, KB, LB, -ONE, A(J,L), LDA, WORK(WORKOFFSET), LB, ONE, X(J,K), LDX)
                            !$omp end task

                            WORKOFFSET = (K-1)*M+(L-1)*KB+1+2*M*M
                            !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                            !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,K))
                            CALL SGEMM('N', 'N', JB, KB, LB, -ONE, B(J,L), LDB, WORK(WORKOFFSET), LB, ONE, X(J,K), LDX)
                            !$omp end task

                            J = J - 1
                        END DO
                    END IF

                    IF ( K .LT. L) THEN
                        J = L -1
                        DO WHILE ( J .GE. K )
                            JH = J
                            J = MAX(K, JH - MB + 1)
                            JB = JH - J + 1
                            IF ( J .GT. 1 ) THEN
                                IF ( A(J,J-1) .NE. ZERO ) THEN
                                    JB = JB - 1
                                    J = J + 1
                                END IF
                            END IF

                            WORKOFFSET = (L-1)*M+(K-1)*LB+1+M*M
                            !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                            !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(K, J))
                            CALL SGEMM('N', 'T', KB, JB , LB, -ONE, WORK(WORKOFFSET), KB, B(J,L), LDB, ONE, X(K,J), LDX)
                            !$omp end task

                            WORKOFFSET = (L-1)*M+(K-1)*LB+1+3*M*M
                            !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                            !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(K, J))
                            CALL SGEMM('N', 'T', KB, JB , LB, -ONE, WORK(WORKOFFSET), KB, A(J,L), LDA, ONE, X(K,J), LDX)
                            !$omp end task

                            J = J - 1
                        END DO
                    END IF

                END IF

                KOLD = K
                K = K - 1
            END DO

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

        L = 1
        LH = 1
        DO WHILE ( LH .LE. M )
            L = LH
            LH = MIN(M, L + MB - 1)
            LB = LH - L + 1
            IF ( LH .LT. M ) THEN
                IF (( A(LH+1,LH) .NE. ZERO  ) ) THEN
                    LB = LB - 1
                    LH = LH - 1
                END IF
            END IF

            K = L
            KH = L
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

                !$omp task firstprivate(TRANSA,TRANSB,M,K,KB,KH,L,LB,LH,INFO1,SCAL) depend(inout:X(K,L)) default(shared)
                CALL ISLA_TGLYAP_SOLVE_BLOCK_T(TRANSA, TRANSB, M, MB, K, KB, KH, L, LB, LH, A, LDA, B, LDB, X, LDX, &
                    & WORK, SCAL, INFO)
                !$omp end task

                ! Prepare update January 2021
                IF ( KH  .LT. M ) THEN
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1
                    !$omp task firstprivate(KB, LB, K, L, KH,WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM("N", "N", KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1+2*M*M
                    !$omp task firstprivate(KB, LB, K, L, KH,WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM("N", "N", KB, LB, LB, ONE, X(K,L), LDX, A(L,L), LDA, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                END IF

                IF ( L.NE.K ) THEN
                    WORKOFFSET = (K-1)*M+(L-1)*KB+1
                    !$omp task firstprivate(KB, LB, K, L, KH,WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM("T", "N", LB, KB, KB, ONE, X(K,L), LDX, B(K,K), LDB, ZERO, WORK(WORKOFFSET), LB)
                    !$omp end task
                    WORKOFFSET = (K-1)*M+(L-1)*KB+1+2*M*M
                    !$omp task firstprivate(KB, LB, K, L, KH,WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM("T", "N", LB, KB, KB, ONE, X(K,L), LDX, A(K,K), LDA, ZERO, WORK(WORKOFFSET), LB)
                    !$omp end task

                END IF

                IF ( LH .LT. M ) THEN
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1+M*M
                    !$omp task firstprivate(KB, LB, K, L, KH,WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM("T", "N", KB, LB, KH, ONE, A(1,K), LDA, X(1,L), LDX, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                    WORKOFFSET = (L-1)*M+(K-1)*LB+1+3*M*M
                    !$omp task firstprivate(KB, LB, K, L, KH,WORKOFFSET) default(shared) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET))
                    CALL SGEMM("T", "N", KB, LB, KH, ONE, B(1,K), LDB, X(1,L), LDX, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                END IF



                ! in current Column
                IF ( KH  .LT. M ) THEN

                    JH = KH+1
                    DO WHILE ( JH .LE. M )
                        J = JH
                        JH = MIN(M, J + MB - 1)
                        JB = JH - J + 1
                        IF ( JH .LT. M ) THEN
                            IF ( A(JH+1, JH) .NE. ZERO ) THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF

                        WORKOFFSET = (L-1)*M+(K-1)*LB+1
                        !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,L))
                        CALL SGEMM("T", "N", JB, LB, KB, -ONE, A(K,J), LDA, WORK(WORKOFFSET), KB, ONE, X(J, L), LDX)
                        !$omp end task

                        WORKOFFSET = (L-1)*M+(K-1)*LB+1+2*M*M
                        !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,L))
                        CALL SGEMM("T", "N", JB, LB, KB, -ONE, B(K,J), LDB, WORK(WORKOFFSET), KB, ONE, X(J, L), LDX)
                        !$omp end task

                        JH = JH + 1
                    END DO
                END IF

                IF ( L.NE.K ) THEN

                    JH = K
                    DO WHILE ( JH .LE. M )
                        J = JH
                        JH = MIN(M, J + MB - 1)
                        JB = JH - J + 1
                        IF ( JH .LT. M ) THEN
                            IF ( A(JH+1, JH) .NE. ZERO ) THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF

                        WORKOFFSET = (K-1)*M+(L-1)*KB+1
                        !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,K))
                        CALL SGEMM("T", "N", JB, KB, LB, -ONE, A(L,J), LDA, WORK(WORKOFFSET), LB, ONE, X(J, K), LDX)
                        !$omp end task

                        WORKOFFSET = (K-1)*M+(L-1)*KB+1+2*M*M
                        !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,K))
                        CALL SGEMM("T", "N", JB, KB, LB, -ONE, B(L,J), LDB, WORK(WORKOFFSET), LB, ONE, X(J, K), LDX)
                        !$omp end task

                        JH = JH + 1
                    END DO

                END IF

                IF ( LH .LT. M ) THEN
                    JH = LH+1
                    DO WHILE ( JH .LE. KH )
                        J = JH
                        JH = MIN(KH, J + MB - 1)
                        JB = JH - J + 1
                        IF ( JH .LT. M ) THEN
                            IF (( A(JH+1,JH) .NE. ZERO  ) ) THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF


                        WORKOFFSET = (L-1)*M+(K-1)*LB+1+M*M
                        !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(K,J))
                        CALL SGEMM("N", "N", KB, JB, LB, -ONE, WORK(WORKOFFSET), KB, B(L,J), LDB, ONE, X(K,J), LDX)
                        !$omp end task

                        WORKOFFSET = (L-1)*M+(K-1)*LB+1+3*M*M
                        !$omp task firstprivate(JB,LB,KB,WORKOFFSET,J,K,L) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(K,J))
                        CALL SGEMM("N", "N", KB, JB, LB, -ONE, WORK(WORKOFFSET), KB, A(L,J), LDA, ONE, X(K,J), LDX)
                        !$omp end task

                        JH = JH +1
                    END DO
                END IF

                KOLD = K
                KH = KH + 1
            END DO

            LOLD = L
            LH = LH + 1
        END DO
        !$omp end master
        !$omp taskwait
        !$omp end parallel


    END IF
    RETURN

END SUBROUTINE

SUBROUTINE ISLA_TGLYAP_SOLVE_BLOCK_N(TRANSA, TRANSB, M, MB, K, KB, KH, L, LB, LH, A, LDA, B, LDB, X, LDX, WORK, SCAL, INFO)
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE OMP_LIB
    IMPLICIT NONE
    CHARACTER TRANSA(1), TRANSB(1)
    INTEGER M, K, KB, KH, L, LB, LH, LDA, LDX, INFO, LDB, MB
    REAL A(LDA, M), X(LDX, M), SCAL, B(LDB,M)
    INTEGER ISOLVER, IONE, MAXNB
    REAL WORK(*), SGN, ONE, ZERO, HALF
    PARAMETER(IONE = 1, SGN = 1.0, ONE = 1.0, ZERO=0.0, HALF=0.5)
    INTEGER WORKOFFSET
    INTEGER INFO1, II,JJ
    INTEGER TID

    EXTERNAL SLA_TGLYAP_L2
    EXTERNAL SLA_TGSYLV_L2
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_128
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_32
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_64
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_96
    EXTERNAL SLA_TGSYLV_L2_REORDER
    EXTERNAL SLA_TGSYLV_L2_COLWISE
    EXTERNAL SLA_TGSYLV_RECURSIVE


    MAXNB = MAX(KB, LB)
    ISOLVER = TGLYAP_ISOLVER()
    INFO1 = 0

    TID =  OMP_GET_THREAD_NUM()
    WORKOFFSET = 4*M*M + (TID)*MAX(MB*MB,4*MB) +1

    IF ( L.EQ.K) THEN
        CALL SLA_TGLYAP_L2 (TRANSA, KB, A(K,K), LDA, B(K,K), LDB, X(K,K), LDX, SCAL, WORK(WORKOFFSET), INFO1)
    ELSE
        SELECT CASE(ISOLVER)
        CASE(1)
            IF ( MAXNB .LE. 32 ) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            ELSE IF ( MAXNB .LE. 64) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            ELSE IF ( MAXNB .LE. 96) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            ELSE IF ( MAXNB .LE. 128) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            ELSE
                CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            END IF
        CASE(2)
            IF ( KB .GT. 128 .OR. LB .GT. 128) THEN
                CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            ELSE
                CALL SLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            END IF

        CASE(3)
            CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        CASE(4)
            CALL SLA_TGSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        CASE(5)
            CALL SLA_TGSYLV_L2_COLWISE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        CASE(6)
            CALL SLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        END SELECT

        IF ( L .NE. K ) THEN
            ! Transpose causes a temporary here
            ! X(L:LH,K:KH) = TRANSPOSE(X(K:KH,L:LH))
            DO II = K, KH
                DO JJ = L, LH
                    X(JJ, II) = X(II,JJ)
                END DO
            END DO
        END IF

    END IF

    IF ( INFO1 .LT. 0 ) THEN
        INFO = 100
    END IF



END SUBROUTINE

SUBROUTINE ISLA_TGLYAP_SOLVE_BLOCK_T(TRANSA, TRANSB, M, MB, K, KB, KH, L, LB, LH, A, LDA, B, LDB, X, LDX, WORK, SCAL, INFO)
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE OMP_LIB
    IMPLICIT NONE
    CHARACTER TRANSA(1), TRANSB(1)
    INTEGER M, MB, K, KB, KH, L, LB, LH, LDA, LDX, INFO, LDB
    REAL A(LDA, M), X(LDX, M), SCAL, B(LDB,M)
    INTEGER ISOLVER, IONE, MAXNB
    REAL WORK(*), SGN, ONE, ZERO, HALF
    PARAMETER(IONE = 1, SGN = 1.0, ONE = 1.0, ZERO=0.0, HALF=0.5)
    INTEGER WORKOFFSET, TID
    INTEGER INFO1, II, JJ

    EXTERNAL SLA_TGLYAP_L2
    EXTERNAL SLA_TGSYLV_L2
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_128
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_32
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_64
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_96
    EXTERNAL SLA_TGSYLV_L2_REORDER
    EXTERNAL SLA_TGSYLV_L2_COLWISE
    EXTERNAL SLA_TGSYLV_RECURSIVE


    MAXNB = MAX(KB, LB)
    ISOLVER = TGLYAP_ISOLVER()
    INFO1=0

    TID =  OMP_GET_THREAD_NUM()
    WORKOFFSET = 4*M*M + TID*MAX(MB*MB,4*MB) +1

    IF ( L.EQ.K) THEN
        CALL SLA_TGLYAP_L2 (TRANSA, KB, A(K,K), LDA, B(K,K), LDB, X(K,K), LDX, SCAL, WORK(WORKOFFSET), INFO1)
    ELSE
        SELECT CASE(ISOLVER)
        CASE(1)
            IF ( MAXNB .LE. 32 ) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            ELSE IF ( MAXNB .LE. 64) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            ELSE IF ( MAXNB .LE. 96) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            ELSE IF ( MAXNB .LE. 128) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            ELSE
                CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            END IF
        CASE(2)
            IF ( KB .GT. 128 .OR. LB .GT. 128) THEN
                CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            ELSE
                CALL SLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
            END IF

        CASE(3)
            CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        CASE(4)
            CALL SLA_TGSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        CASE(5)
            CALL SLA_TGSYLV_L2_COLWISE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        CASE(6)
            CALL SLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & B(K,K), LDB, A(L,L), LDA,  X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        END SELECT

        ! X(L:LH,K:KH) = TRANSPOSE(X(K:KH,L:LH))

        DO II = K, KH
            DO JJ = L, LH
                X(JJ, II) = X(II,JJ)
            END DO
        END DO
    END IF
    IF ( INFO1 .LT. 0 ) THEN
        INFO = 100
    END IF

END SUBROUTINE





