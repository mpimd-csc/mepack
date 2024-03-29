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


!> \brief Level-3 Bartels-Stewart Algorithm for the discrete time Sylvester equation with DAG based parallelization.
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TRSYLV2_DAG ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
!           CHARACTER(1) TRANSA, TRANSB
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
!> SLA_TRSYLV2_DAG solves a discrete time Sylvester equation of the following forms
!>
!>    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)
!>
!> or
!>
!>    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)
!>
!> where A is a M-by-M quasi upper triangular matrix and B is a N-by-N upper
!> quasi triangular matrices. The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrices A and B are created by SGEES from LAPACK.
!> \endverbatim
!> \remark The algorithm is implemented using BLAS level 3 operations and OpenMP 4.0 DAG scheduling.
!
!> \attention Due to the parallel nature of the algorithm the scaling is not applied to the right hand. If the problem
!> is ill-posed and scaling appears you have to solve the equation again with a solver with complete scaling support
!> like the SLA_TRSYLV2_L3 routine.
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
!>          If SCALE .NE. 1 the problem is no solved correctly in this case
!>          one have to use an other solver.
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
!> \ingroup sgltrsylv
!
SUBROUTINE SLA_TRSYLV2_DAG ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE OMP_LIB
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    REAL SGN, SCALE
    INTEGER M, N, LDA, LDB, LDX, INFO
    REAL A(LDA, M), B(LDB, N) , X(LDX, N)
    REAL WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KOLD, LOLD
    INTEGER INFO1, WORKOFFSET
    REAL SCAL
    INTEGER MB, NB, ISOLVER, MAXNB
    INTEGER J, JH, JB
    INTEGER NTH, BLO

    LOGICAL BTRANSA, BTRANSB
    REAL ZERO , ONE
    PARAMETER (ZERO = 0.0, ONE = 1.0)
    INTEGER IZERO, IONE
    PARAMETER(IZERO = 0, IONE = 1 )

    ! EXTERNAL FUNCTIONS
    EXTERNAL SGEMM
    EXTERNAL SLA_TRSYLV2_L2
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_128
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_32
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_64
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_96
    EXTERNAL SLA_TRSYLV2_L2_REORDER
    EXTERNAL SLA_TRSYLV2_L2_UNOPT
    EXTERNAL SLA_TRSYLV2_RECURSIVE
    EXTERNAL ISLA_TRSYLV2_SOLVE_BLOCK
    EXTERNAL LSAME
    EXTERNAL XERROR_HANDLER
    LOGICAL LSAME

    MB = TRSYLV2_BLOCKSIZE_MB(M,N)
    NB = TRSYLV2_BLOCKSIZE_NB(M,N)
    ISOLVER = TRSYLV2_ISOLVER()

    MAXNB = MAX(MB, NB)

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
    ELSE IF ( MB .LT. 2 ) THEN
        INFO = -30
    ELSE IF ( NB .LT. 2 ) THEN
        INFO = -31
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_TRSYLV2_DAG', -INFO)
        RETURN
    END IF

    BLO = MAX(MB*NB, 2*MAX(MB,NB))
    NTH = OMP_GET_MAX_THREADS()
    IF ( ININFO .EQ. -1) THEN
        INFO = 2*M*N + NTH * BLO
        RETURN
    END IF



    INFO = 0
    INFO1 = 0

    ! Quick Return
    IF ( M .EQ. 0 .OR. N .EQ. 0 ) THEN
        RETURN
    END IF

    IF ( M .LE. MB .AND. N .LE. NB ) THEN
        KB = M
        LB = N
        SELECT CASE(ISOLVER)
        CASE(1)
            IF ( MAXNB .LE. 32 ) THEN
                CALL SLA_TRSYLV2_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    &   X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 64) THEN
                CALL SLA_TRSYLV2_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    &   X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 96) THEN
                CALL SLA_TRSYLV2_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    &   X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 128) THEN
                CALL SLA_TRSYLV2_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    &   X, LDX, SCAL, WORK, INFO1)
            ELSE
                CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    &   X, LDX, SCAL, WORK, INFO1)
            END IF
        CASE(2)
            IF ( M .GT. 128 .OR. N .GT. 128) THEN
                CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    &   X, LDX, SCAL, WORK, INFO1)
            ELSE
                CALL SLA_TRSYLV2_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    &   X, LDX, SCAL, WORK, INFO1)
            END IF
        CASE(3)
            CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                &   X, LDX, SCAL, WORK, INFO1)
        CASE(4)
            CALL SLA_TRSYLV2_L2(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                &   X, LDX, SCAL, WORK, INFO1)
        CASE(5)
            CALL SLA_TRSYLV2_L2_UNOPT(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                &   X, LDX, SCAL, WORK, INFO1)
        CASE(6)
            CALL SLA_TRSYLV2_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                &   X, LDX, SCAL, WORK, INFO1)
        END SELECT
        IF ( INFO1 .NE. 0 ) THEN
            !$omp atomic
            INFO = INFO + 1
        END IF
        IF ( SCAL .NE. ONE ) THEN
            !$omp atomic
            SCALE = SCALE * SCAL
        END IF
        RETURN
    END IF

    ! Prepare
    INFO = 0

    IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB ) THEN
        ! ( N, N )

        !$omp parallel default(shared)
        !$omp master
        LOLD = 1
        KOLD = 1
        L = 1
        LH = 1
        DO WHILE ( LH .LE. N )
            L = LH
            LH = MIN(N, L + NB - 1)
            LB = LH - L + 1
            IF ( LH .LT. N ) THEN
                IF ( ( B(LH+1,LH) .NE. ZERO ) ) THEN
                    LB = LB - 1
                    LH = LH - 1
                END IF
            END IF

            K = M
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

                !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                !$omp& depend(inout: X(K,L)) default(shared)
                CALL ISLA_TRSYLV2_SOLVE_BLOCK(TRANSA, TRANSB, M,N,BLO,K,KH,KB,L,LH,LB,A,LDA,B,LDB,&
                    & X,LDX,SGN,SCAL, WORK, INFO1 )
                !$omp end task



                ! Prepare Update January 2021
                IF ( K .GT. 1 ) THEN
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) &
                    !$omp depend(in:X(K,L)) depend(out:WORK(WORKOFFSET)) default(shared)
                    CALL SGEMM('N', 'N', KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

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


                        WORKOFFSET = (L-1)*M+(K-1)*LB+1
                        !$omp task firstprivate(J,K,L,JB,KB,LB,WORKOFFSET) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,L))
                        CALL SGEMM('N', 'N', JB, LB, KB, -ONE, A(J,K), LDA, WORK(WORKOFFSET), KB, ONE, X(J,L), LDX)
                        !$omp end task


                        J = J - 1
                    END DO


                END IF

                IF ( LH .LT. N  ) THEN
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1+M*N
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) &
                    !$omp depend(in:X(K,L)) depend(out:WORK(WORKOFFSET)) default(shared)
                    CALL SGEMM("N", "N", KB, LB, M-K+1, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                    J = LH+1
                    JH = LH +1
                    DO WHILE (JH .LE. N)
                        J = JH
                        JH = MIN(N, J + NB - 1)
                        JB = JH - J + 1

                        IF ( JH .LT. N ) THEN
                            IF (( B(JH+1,JH) .NE. ZERO ) ) THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF


                        WORKOFFSET = (L-1)*M+(K-1)*LB+1+M*N

                        !$omp task firstprivate(J,K,L,JB,KB,LB,WORKOFFSET) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(K,J))
                        CALL SGEMM("N", "N", KB, JB, LB, -ONE, WORK(WORKOFFSET), KB, B(L,J), LDB, ONE, X(K,J), LDX)
                        !$omp end task
                        JH = JH + 1
                    END DO

                END IF

                KOLD = K
                K = K - 1
            END DO


            LOLD = L
            LH = LH + 1
        END DO

        !$omp end master
        !$omp taskwait
        !$omp end parallel

    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N,T)
        !$omp parallel default(shared)
        !$omp master
        KOLD = 1
        LOLD = 1
        L = N
        DO WHILE ( L .GT. 0 )
            IF (L .EQ. N .AND. MOD(N,NB) .NE. IZERO) THEN
                LH = L
                L = (L/NB)*NB + 1
            ELSE
                LH = L
                L = MAX(1, LH - NB + 1)
            END IF

            LB = LH - L + 1
            IF ( L .GT. 1 ) THEN
                IF ((B(L,L-1) .NE. ZERO ) ) THEN
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


            K = M
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

                !$omp task depend(inout: X(K,L)) firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1)&
                !$omp& default(shared)
                CALL ISLA_TRSYLV2_SOLVE_BLOCK(TRANSA, TRANSB, M,N,BLO,K,KH,KB,L,LH,LB,A,LDA,B,LDB,&
                    & X,LDX,SGN,SCAL, WORK, INFO1 )

                !$omp end task

                IF ( K .GT. 1 ) THEN
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) &
                    !$omp depend(in:X(K,L)) depend(out:WORK(WORKOFFSET)) default(shared)
                    CALL SGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                    J = K-1
                    DO WHILE (J .GT. 0)
                        JH = J
                        J = MAX(1, JH - MB + 1)
                        JB = JH - J + 1
                        IF ( J .GT. 1) THEN
                            IF ( A(J,J-1) .NE. ZERO ) THEN
                                JB = JB -1
                                J = J + 1
                            END IF
                        END IF


                        WORKOFFSET = (L-1)*M+(K-1)*LB+1

                        !$omp task firstprivate(J,K,L,JB,KB,LB,WORKOFFSET) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,L))
                        CALL SGEMM('N', 'N', JB, LB, KB, -ONE, A(J,K), LDA, WORK(WORKOFFSET), KB, ONE, X(J,L), LDX)
                        !$omp end task

                        J = J - 1
                    END DO


                END IF


                IF ( L .GT. 1 ) THEN
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1 + M*N

                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) &
                    !$omp depend(in:X(K,L)) depend(out:WORK(WORKOFFSET)) default(shared)
                    CALL SGEMM("N", "N", KB, LB, M-K+1, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                    J = L-1
                    DO WHILE (J .GT. 0)
                        JH = J
                        J = MAX(1, JH - NB + 1)
                        JB = JH - J + 1
                        IF ( J .GT. 1 ) THEN
                            IF (( B(J,J-1) .NE. ZERO )) THEN
                                JB = JB -1
                                J = J + 1
                            END IF
                        END IF

                        WORKOFFSET = (L-1)*M+(K-1)*LB+1 + M*N

                        !$omp task firstprivate(J,K,L,JB,KB,LB,WORKOFFSET) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(K,J))
                        CALL SGEMM('N', 'T', KB, JB, LB, -ONE, WORK(WORKOFFSET), KB, B(J,L), LDB, ONE, X(K,J), LDX)
                        !$omp end task

                        J = J - 1
                    END DO



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

    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
        ! (T, N)

        !$omp parallel default(shared)
        !$omp master
        KOLD = 1
        LOLD = 1

        L = 1
        LH = 1
        DO WHILE ( LH .LE. N )
            L = LH
            LH = MIN(N, L + NB - 1)
            LB = LH - L + 1
            IF ( LH .LT. N ) THEN
                IF (( B(LH+1,LH) .NE. ZERO ) ) THEN
                    LB = LB - 1
                    LH = LH - 1
                END IF
            END IF

            K = 1
            KH = 1
            DO WHILE ( KH .LE. M )
                K = KH
                KH = MIN(M, K + MB - 1)
                KB = KH - K + 1
                IF ( KH .LT. M ) THEN
                    IF (A(KH+1, KH) .NE. ZERO ) THEN
                        KB = KB - 1
                        KH = KH - 1
                    END IF
                END IF

                !$omp task depend(inout: X(K,L)) firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1)&
                !$omp& default(shared)
                CALL ISLA_TRSYLV2_SOLVE_BLOCK(TRANSA, TRANSB, M,N,BLO,K,KH,KB,L,LH,LB,A,LDA,B,LDB,&
                    & X,LDX,SGN,SCAL, WORK, INFO1 )

                !$omp end task

                IF ( KH  .LT. M ) THEN
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) &
                    !$omp& depend(in:X(K,L)) depend(out:WORK(WORKOFFSET)) default(shared)
                    CALL SGEMM("N", "N", KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                    J = KH+1
                    JH = KH +1
                    DO WHILE (JH .LE. M)
                        J = JH
                        JH = MIN(M, J + MB - 1)
                        JB = JH - J + 1
                        IF ( JH .LT. M ) THEN
                            IF (( A(JH+1,JH) .NE. ZERO  ) ) THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF


                        WORKOFFSET = (L-1)*M+(K-1)*LB+1
                        !$omp task firstprivate(J,K,L,JB,KB,LB,WORKOFFSET) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,L))
                        CALL SGEMM("T", "N", JB, LB, KB, -ONE, A(K,J), LDA, WORK(WORKOFFSET), KB, ONE, X(J, L), LDX)
                        !$omp end task

                        JH = JH + 1
                    END DO

                END IF

                IF ( LH .LT. N ) THEN
                    WORKOFFSET = (L-1)*M+(K-1)*LB+1 + M * N
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET, KH) &
                    !$omp depend(in:X(K,L)) depend(out:WORK(WORKOFFSET)) default(shared)
                    CALL SGEMM("T", "N", KB, LB, KH, ONE, A(1,K), LDA, X(1,L), LDX, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                    J = LH+1
                    JH = LH +1
                    DO WHILE (JH .LE. N)
                        J = JH
                        JH = MIN(N, J + NB - 1)
                        JB = JH - J + 1

                        IF ( JH .LT. N ) THEN
                            IF ( ( B(JH+1,JH) .NE. ZERO ) ) THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF


                        WORKOFFSET = (L-1)*M+(K-1)*LB+1 + M * N
                        !$omp task firstprivate(KB,JB,LB,WORKOFFSET, K,J,L) &
                        !$omp& default(shared) depend(in:WORK(WORKOFFSET)) depend(inout:X(K,J))
                        CALL SGEMM("N", "N", KB, JB, LB, -ONE, WORK(WORKOFFSET), KB, B(L,J), LDB, ONE, X(K,J), LDX)
                        !$omp end task

                        JH = JH + 1
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


    ELSE IF ( BTRANSA .AND. BTRANSB) THEN
        ! (T, T)
        !$omp parallel default(shared)
        !$omp master
        KOLD = 1
        LOLD = 1

        L = N
        DO WHILE ( L .GT. 0 )
            IF (L .EQ. N .AND. MOD(N,NB) .NE. IZERO) THEN
                LH = L
                L = (L/NB)*NB + 1
            ELSE
                LH = L
                L = MAX(1, LH - NB + 1)
            END IF

            LB = LH - L + 1
            IF ( L .GT. 1 ) THEN
                IF ( B(L,L-1) .NE. ZERO ) THEN
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

                !$omp task depend(inout: X(K,L)) firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1)&
                !$omp& default(shared)
                CALL ISLA_TRSYLV2_SOLVE_BLOCK(TRANSA, TRANSB, M,N,BLO,K,KH,KB,L,LH,LB,A,LDA,B,LDB,&
                    & X,LDX,SGN,SCAL, WORK, INFO1 )
                !$omp end task


                IF ( KH .LT. M  ) THEN

                    WORKOFFSET = (L-1)*M+(K-1)*LB+1
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET) &
                    !$omp depend(in:X(K,L)) depend(out:WORK(WORKOFFSET)) default(shared)
                    CALL SGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task
                    J = KH+1
                    JH = KH +1
                    DO WHILE (JH .LE. M)
                        J = JH
                        JH = MIN(M, J + MB - 1)
                        JB = JH - J + 1
                        IF ( JH .LT. M ) THEN
                            IF (( A(JH+1,JH) .NE. ZERO  ) ) THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF

                        WORKOFFSET = (L-1)*M+(K-1)*LB+1
                        !$omp task firstprivate(J,K,L,JB,KB,LB,WORKOFFSET) default(shared) &
                        !$omp& depend(in:WORK(WORKOFFSET)) depend(inout:X(J,L))
                        CALL SGEMM('T', 'N', JB, LB, KB, -ONE, A(K, J), LDA, WORK(WORKOFFSET), KB, ONE, X(J,L),LDX)
                        !$omp end task


                        JH = JH + 1
                    END DO

                END IF

                IF ( L .GT. IONE ) THEN

                    WORKOFFSET = (L-1)*M+(K-1)*LB+1 + M*N
                    !$omp task firstprivate(KB, LB, K, L, WORKOFFSET,KH) &
                    !$omp depend(in:X(K,L)) depend(out:WORK(WORKOFFSET)) default(shared)
                    CALL SGEMM('T', 'N', KB, LB, KH, ONE, A(1,K), LDA, X(1,L), LDX, ZERO, WORK(WORKOFFSET), KB)
                    !$omp end task

                    J = L-1
                    DO WHILE (J .GT. 0)
                        JH = J
                        J = MAX(1, JH - NB + 1)
                        JB = JH - J + 1
                        IF ( J .GT. 1 ) THEN
                            IF ( ( B(J,J-1) .NE. ZERO )) THEN
                                JB = JB -1
                                J = J + 1
                            END IF
                        END IF

                        WORKOFFSET = (L-1)*M+(K-1)*LB+1 + M*N
                        !$omp task firstprivate(KB,JB,LB,WORKOFFSET, K,J,L) &
                        !$omp& default(shared) depend(in:WORK(WORKOFFSET)) depend(inout:X(K,J))

                        CALL SGEMM('N', 'T', KB, JB, LB, -ONE, WORK(WORKOFFSET), KB, B(J, L), LDB, ONE, X(K,J), LDX)
                        !$omp end task


                        J = J - 1
                    END DO


                END IF



                KOLD = K
                KH = KH + 1

            END DO


            LOLD = L
            L = L - 1
        END DO
        !$omp end master
        !$omp taskwait
        !$omp end parallel

    END IF
    RETURN

END SUBROUTINE

!> \brief Internal Block
!> \internal
SUBROUTINE ISLA_TRSYLV2_SOLVE_BLOCK(TRANSA, TRANSB, M,N,BLO,K,KH,KB,L,LH,LB,A,LDA,B,LDB,X,LDX,SGN,SCAL, &
        & WORK, INFO1 )
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE OMP_LIB
    IMPLICIT NONE
    CHARACTER TRANSA(1), TRANSB(1)
    INTEGER M,N,BLO,K,KH,KB,L,LH,LB, LDA, LDB, LDX, INFO1
    REAL A(LDA,M), B(LDB,N), X(LDX,N), SGN, SCAL, WORK(*)
    INTEGER MAXNB, ISOLVER
    REAL ONE, ZERO
    PARAMETER (ONE=1.0, ZERO = 0.0)
    INTEGER WORKOFFSET, TID

    EXTERNAL SLA_TRSYLV2_L2
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_128
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_32
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_64
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_96
    EXTERNAL SLA_TRSYLV2_L2_REORDER
    EXTERNAL SLA_TRSYLV2_L2_UNOPT
    EXTERNAL SLA_TRSYLV2_RECURSIVE


    MAXNB = MAX(KB, LB)
    ISOLVER = TRSYLV2_ISOLVER()

    KH = KH
    LH = LH

    TID = OMP_GET_THREAD_NUM()
    WORKOFFSET = 2*M*N+1 + TID * BLO

    SELECT CASE(ISOLVER)
    CASE(1)
        IF ( MAXNB .LE. 32 ) THEN
            CALL SLA_TRSYLV2_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        ELSE IF ( MAXNB .LE. 64) THEN
            CALL SLA_TRSYLV2_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        ELSE IF ( MAXNB .LE. 96) THEN
            CALL SLA_TRSYLV2_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        ELSE IF ( MAXNB .LE. 128) THEN
            CALL SLA_TRSYLV2_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        ELSE
            CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        END IF
    CASE(2)
        IF ( KB .GT. 128 .OR. LB .GT. 128) THEN
            CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        ELSE
            CALL SLA_TRSYLV2_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
        END IF
    CASE(3)
        CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
            & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
    CASE(4)
        CALL SLA_TRSYLV2_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
            & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
    CASE(5)
        CALL SLA_TRSYLV2_L2_UNOPT(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
            & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
    CASE(6)
        CALL SLA_TRSYLV2_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
            & X(K,L), LDX, SCAL, WORK(WORKOFFSET), INFO1)
    END SELECT


END SUBROUTINE
