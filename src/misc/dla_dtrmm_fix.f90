!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) Martin Koehler, 2017-2023

! This is an internal routine of MEPACK and not intended to be used from outside.
! The routine computes
!
!    C <=  BETA * C + ALPHA * op(A) * B
!
! where A is a upper (quasi) triangular matrix. The matrices B and C are general
! ones.
!
! The routine uses OpenMP 4 for acceleration.
!
!
SUBROUTINE DLA_DTRMM_FIX(TRANS, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    USE MEPACK_OPTIONS_OPENMP
    IMPLICIT NONE
    INTEGER, INTENT(IN)  ::  M, N, LDA, LDB, LDC
    DOUBLE PRECISION, INTENT(IN) ::  ALPHA, BETA
    DOUBLE PRECISION, INTENT(IN) ::  A(LDA, *), B(LDB, *)
    DOUBLE PRECISION, INTENT(INOUT) ::  C(LDC,*)
    CHARACTER(1), INTENT(IN) ::  TRANS

    EXTERNAL LSAME, DGEMM, DAXPY
    LOGICAL LSAME

    DOUBLE PRECISION ZERO, DONE
    INTEGER L, LH, LB, NB, K, KB, KH, J, JB
    PARAMETER(ZERO = 0.0D0, DONE = 1.0D0)

    IF ( M .LT. 1536 .OR. N .LT. 1536) THEN
        NB = 64
    ELSE IF ( M .GT. 4000 .AND. N.GT.4000) THEN
        NB = 256
    ELSE
        NB = 128
    END IF

    IF ( OPENMP_ENABLED() .NE. 0 ) THEN
        IF ( LSAME(TRANS, "N") ) THEN
            !$omp parallel do private(LH, LB, K, KB, KH)
            DO L=1, M, NB
                LH = MIN(M, L + NB  - 1)
                LB = LH - L + 1

                DO K = 1, N, NB
                    KB = MIN ( NB, N - K + 1)
                    KH = K + KB - 1
                    IF ( BETA .EQ. ZERO ) THEN
                        C(L:LH, K:KH) = ZERO
                    ELSE
                        C(L:LH, K:KH) = BETA * C(L:LH, K:KH)
                    END IF

                    CALL DGEMM("N", "N", LB, KB, M-L+1, ALPHA, A(L,L), LDA, B(L,K), LDB, DONE, C(L,K), LDC)

                    IF ( L.GT.1 ) THEN
                        IF (A(L,L-1) .NE. ZERO ) THEN
                            C(L,K:KH) = (ALPHA*A(L,L-1)) * B(L-1, K:KH) + C(L,K:KH)
                        END IF
                    END IF
                END DO
            END DO
            !$omp end parallel do

        ELSE
            !$omp parallel do private(LH,LB,K,KB,KH)
            DO L=1, M, NB
                LH = MIN(M, L + NB  - 1)
                LB = LH - L + 1

                DO K = 1, N, NB
                    KB = MIN ( NB, N - K + 1)
                    KH = K + KB - 1

                    IF ( BETA .EQ. ZERO ) THEN
                        C(L:LH, K:KH) = ZERO
                    ELSE
                        C(L:LH, K:KH) = BETA * C(L:LH, K:KH)
                    END IF

                    CALL DGEMM("T", "N", LB, KB, LH, ALPHA, A(1,L), LDA, B(1,K), LDB, DONE, C(L,K), LDC)

                    IF ( LH.LT.M ) THEN
                        IF (A(LH+1,LH) .NE. ZERO ) THEN
                            C(LH,K:KH) = (ALPHA*A(LH+1,LH)) * B(LH+1, K:KH) + C(LH,K:KH)
                        END IF
                    END IF
                END DO
            END DO
            !$omp end parallel do
        END IF
!   This code is faster but leads to random deadlocks using Intel's Compilers
!
!         IF ( LSAME(TRANS, "N") ) THEN
!             !$omp parallel
!             !$omp master
!             DO L=1, M, NB
!                 LH = MIN(M, L + NB  - 1)
!                 LB = LH - L + 1
!
!                 DO K = 1, N, NB
!                     KB = MIN ( NB, N - K + 1)
!                     KH = K + KB - 1
!                     !$omp task firstprivate(L, LH, K, KH) depend(inout:C(L,K))
!                     IF ( BETA .EQ. ZERO ) THEN
!                         C(L:LH, K:KH) = ZERO
!                     ELSE
!                         C(L:LH, K:KH) = BETA * C(L:LH, K:KH)
!                     END IF
!                     !$omp end task
!
!
!                     DO J = L, M, NB
!                         JB = MIN(NB, M - J + 1)
!                         !$omp task firstprivate(LB, KB, L, K, KH, J, JB) depend(inout:C(L,K))
!                         CALL DGEMM("N", "N", LB, KB, JB, ALPHA, A(L,J), LDA, B(J,K), LDB, DONE, C(L,K), LDC)
!                         !$omp end task
!                     END DO
!
!                     IF ( L.GT.1 ) THEN
!                         IF (A(L,L-1) .NE. ZERO ) THEN
!                             !$omp task firstprivate(L,K,KH) depend(inout:C(L,K))
!                             C(L,K:KH) = (ALPHA*A(L,L-1)) * B(L-1, K:KH) + C(L,K:KH)
!                             !$omp end task
!                         END IF
!                     END IF
!                 END DO
!             END DO
!             !$omp end master
!             !$omp taskwait
!             !$omp end parallel
!         ELSE
!             !$omp parallel
!             !$omp master
!             DO L=1, M, NB
!                 LH = MIN(M, L + NB  - 1)
!                 LB = LH - L + 1
!
!                 DO K = 1, N, NB
!                     KB = MIN ( NB, N - K + 1)
!                     KH = K + KB - 1
!
!                     !$omp task firstprivate(L, LH, K, KH) depend(inout:C(L,K))
!                     IF ( BETA .EQ. ZERO ) THEN
!                         C(L:LH, K:KH) = ZERO
!                     ELSE
!                         C(L:LH, K:KH) = BETA * C(L:LH, K:KH)
!                     END IF
!
!                     !$omp end task
!
!                     DO J = 1, LH, NB
!                         JB = MIN(NB, LH - J + 1)
!                         !$omp task firstprivate(LB, KB, L, K, KH, J, JB) depend(inout:C(L,K))
!                         CALL DGEMM("T", "N", LB, KB, JB, ALPHA, A(J,L), LDA, B(J,K), LDB, DONE, C(L,K), LDC)
!                         !$omp end task
!
!                     END DO
!
!                     IF ( LH.LT.M ) THEN
!                         IF (A(LH+1,LH) .NE. ZERO ) THEN
!                             !$omp task firstprivate(L,LH,K,KH) depend(inout:C(L,K))
!                             C(LH,K:KH) = (ALPHA*A(LH+1,LH)) * B(LH+1, K:KH) + C(LH,K:KH)
!                             !$omp end task
!                         END IF
!                     END IF
!                 END DO
!             END DO
!             !$omp end master
!             !$omp taskwait
!             !$omp end parallel
!         END IF

    ELSE
        IF ( LSAME(TRANS, "N") ) THEN
            DO L=1, M, NB
                LH = MIN(M, L + NB  - 1)
                LB = LH - L + 1

                DO K = 1, N, NB
                    KB = MIN ( NB, N - K + 1)
                    KH = K + KB - 1
                    IF ( BETA .EQ. ZERO ) THEN
                        C(L:LH, K:KH) = ZERO
                    ELSE
                        C(L:LH, K:KH) = BETA * C(L:LH, K:KH)
                    END IF

                    DO J = L, M, NB
                        JB = MIN(NB, M - J + 1)
                        CALL DGEMM("N", "N", LB, KB, JB, ALPHA, A(L,J), LDA, B(J,K), LDB, DONE, C(L,K), LDC)
                    END DO

                    IF ( L.GT.1 ) THEN
                        IF (A(L,L-1) .NE. ZERO ) THEN
                            C(L,K:KH) = (ALPHA*A(L,L-1)) * B(L-1, K:KH) + C(L,K:KH)
                        END IF
                    END IF
                END DO
            END DO
        ELSE
            DO L=1, M, NB
                LH = MIN(M, L + NB  - 1)
                LB = LH - L + 1

                DO K = 1, N, NB
                    KB = MIN ( NB, N - K + 1)
                    KH = K + KB - 1

                    IF ( BETA .EQ. ZERO ) THEN
                        C(L:LH, K:KH) = ZERO
                    ELSE
                        C(L:LH, K:KH) = BETA * C(L:LH, K:KH)
                    END IF

                    DO J = 1, LH, NB
                        JB = MIN(NB, LH - J + 1)
                        CALL DGEMM("T", "N", LB, KB, JB, ALPHA, A(J,L), LDA, B(J,K), LDB, DONE, C(L,K), LDC)
                    END DO

                    IF ( LH.LT.M ) THEN
                        IF (A(LH+1,LH) .NE. ZERO ) THEN
                            C(LH,K:KH) = (ALPHA*A(LH+1,LH)) * B(LH+1, K:KH) + C(LH,K:KH)
                        END IF
                    END IF
                END DO
            END DO
        END IF
    END IF
END SUBROUTINE

