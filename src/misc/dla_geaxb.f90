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

! This is an internal MEPACK routine and not intended to be called from outside.
!
! The routine computes
!
!  C <= BETA * C  + ALPHA * opA(A) * X * op(B)
!
! Thereby, it takes the size of A and B into account and reorders the computation
! such that the number of required operations is minimized.
!
SUBROUTINE DLA_GEAXB(TRANSA, TRANSB, FORMA, FORMB, M, N, K, L, ALPHA, A, LDA, X, LDX, B, LDB, BETA, C, LDC, WORK)
    IMPLICIT NONE
    ! Input Arguments

    CHARACTER(1) TRANSA, TRANSB, FORMA, FORMB
    INTEGER M,N,K,L, LDA, LDX, LDB, LDC
    DOUBLE PRECISION ALPHA, BETA
    DOUBLE PRECISION A(LDA,*), B(LDB,*), C(LDC,*), X(LDX,*)
    DOUBLE PRECISION WORK(*)


    ! Local Variables
    DOUBLE PRECISION MULT1A, MULT1B, MULT2A, MULT2B
    DOUBLE PRECISION ONE, ZERO, TWO
    PARAMETER (ONE = 1.0D0, ZERO = 0.0D0, TWO = 2.0D0)

    ! External CODE
    EXTERNAL DGEMM
    EXTERNAL DLA_QTRMM2
    EXTERNAL DTRMM
    EXTERNAL LSAME
    LOGICAL LSAME

    ! Check Input
    ! Skipped for Performance Reasons


    ! Quick Return
    IF (M .EQ. 0 .OR. N.EQ.0) THEN
        RETURN
    ENDIF


    IF (FORMA.eq.'F') THEN
        MULT1A = (TWO*M)*N*K
        MULT2B = (TWO*M)*N*L
    ELSE
        MULT1A = (ONE*M)*M*K+(TWO*M)*(N-M)*K
        MULT2B = (ONE*M)*N*L+(TWO*M)*(N-M)*L
    ENDIF

    IF (FORMB.eq.'F') THEN
        MULT1B = (TWO*M)*K*L
        MULT2A = (TWO*N)*K*L
    ELSE
        MULT1B = (ONE*M)*L*L+(TWO*M)*L*(K-L)
        MULT2A = (ONE*N)*L*L+(TWO*N)*L*(K-L)
    ENDIF

    IF ((MULT1A+MULT1B) .LT. (MULT2A+MULT2B)) THEN

        CALL DGEMM_LEFT (FORMA, TRANSA, M, K, N, ALPHA, A, LDA, X, LDX, ZERO, WORK, M,  WORK(M*K+1))
        CALL DGEMM_RIGHT(FORMB, TRANSB, M, L, K, ONE, WORK, M, B, LDB, BETA, C, LDC, WORK(M*K+1))

    ELSE
        CALL DGEMM_RIGHT(FORMB, TRANSB, N, L, K, ALPHA, X, LDX, B, LDB, ZERO, WORK, N, WORK(N*L+1))
        CALL DGEMM_LEFT (FORMA, TRANSA, M, L, N, ONE, A, LDA, WORK, N, BETA, C, LDC, WORK(N*L+1))
    END IF
    RETURN

CONTAINS

    SUBROUTINE DGEMM_LEFT(FORMA, TRANSA, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, WORK)
        IMPLICIT NONE

        INTEGER M,N,K,LDA,LDB,LDC
        DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *)
        DOUBLE PRECISION WORK(*)
        DOUBLE PRECISION ALPHA, BETA
        CHARACTER(1) TRANSA,FORMA

        ! Local Variables
        DOUBLE PRECISION ONE
        PARAMETER (ONE = 1.0D0)
        INTEGER I,J

        ! External Functions
        EXTERNAL DGEMM, DTRMM, DLA_QTRMM2

        IF (FORMA .EQ.'F' .OR. MAX(MAX(M,N),K) .LT. 40) THEN ! .or.MAX(MAX(M,N),K).lt.MACHINE(4)) THEN
            CALL DGEMM(TRANSA, 'N', M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
            RETURN
        END IF
        IF (TRANSA.EQ.'N') THEN
            IF (K.GE.M) THEN
                !     K     N
                !   AAAaa BBBB CCCC
                ! M  AAaa BBBB CCCC
                !     Aaa BBBB CCCC
                !         bbbb
                !         bbbb
                IF (BETA.EQ.0D0) THEN
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = B(I,J)
                        END DO
                    END DO
                    CALL DTRMM('L','U','N','N', M, N, ALPHA, A, LDA, C, LDC)
                    IF (FORMA .EQ. 'Q') THEN
                        CALL DLA_QTRMM2('L','N', M, N, ALPHA, A, LDA, B, LDB, C, LDC)
                    END IF
                ELSE
                    DO J=1,N
                        DO I=1,M
                            WORK(I+(J-1)*M) = B(I,J)
                        END DO
                    END DO
                    CALL DTRMM('L','U','N','N', M, N, ALPHA, A, LDA, WORK, M)
                    IF (FORMA .EQ.'Q') THEN
                        CALL DLA_QTRMM2('L','N', M, N, ALPHA, A, LDA, B, LDB, WORK, M)
                    END IF
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (K.GT.M) THEN
                    CALL DGEMM('N','N', M, N, K-M, ALPHA, A(1,M+1), LDA, B(M+1,1), LDB, ONE, C, LDC)
                END IF
            ELSE !(K>=M)
                !     K     N
                !   aaa  BBBB cccc
                ! M aaa  BBBB cccc
                !   AAA  BBBB CCCC
                !    AA       CCCC
                !     A       CCCC
                IF (BETA.EQ.0D0) THEN
                    DO J=1,N
                        DO I=1,K
                            C(I+M-K-1,J) = B(I,J)
                        END DO
                    END DO
                    CALL DTRMM('L','U','N','N', K, N, ALPHA, A(M-K+1,1), LDA, C(M-K+1,1), LDC)
                    IF (FORMA .EQ. 'Q') THEN
                        CALL DLA_QTRMM2('L','N', K, N, ALPHA, A(M-K+1,1), LDA, B, LDB, C(M-K+1,1), LDC)
                    END IF
                ELSE
                    DO J=1,N
                        DO I=1,K
                            WORK(I+(J-1)*K) = B(I,J)
                        END DO
                    END DO
                    CALL DTRMM('L','U','N','N', K, N, ALPHA, A(M-K+1,1), LDA, WORK, K)
                    IF (FORMA.eq.'Q') THEN
                        CALL DLA_QTRMM2('L','N', K, N, ALPHA, A(M-K+1,1), LDA, B, LDB, WORK, K)
                    END IF
                    DO J=1,N
                        DO I=1,K
                            C(I+M-K-1,J) = BETA * C(I+M-K-1,J) + WORK(I+J*K)
                        END DO
                    END DO
                END IF
                IF (M.GT.K) THEN
                    CALL DGEMM('N','N', M-K, N, K, ALPHA, A(1,1), LDA, B(1,1), LDB, BETA, C(1,1), LDC )
                END IF
            END IF
        ELSE !TRANSA == 'T'
            IF (M.GE.K) THEN
                !   K
                !  A   BBBB  CCCC
                !M AA  BBBB  CCCC
                !  AAA BBBB  CCCC
                !  aaa       cccc
                !  aaa       cccc
                ! CALL DGEMM(TRANSA, 'N', M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
                IF (BETA.EQ.0D0) THEN
                    DO J=1,N
                        DO I=1,K
                            C(I,J) = B(I,J)
                        END DO
                    END DO
                    CALL DTRMM('L','U','T', 'N', K, N, ALPHA, A, LDA, C, LDC)
                    IF (FORMA .EQ. 'Q') THEN
                        CALL DLA_QTRMM2('L','T', K, N, ALPHA, A, LDA, B, LDB, C, LDC)
                    END IF
                ELSE
                    DO J=1,N
                        DO I=1,K
                            WORK(I+(J-1)*K) = B(I,J)
                        END DO
                    END DO
                    CALL DTRMM('L','U','T', 'N', K, N, ALPHA, A, LDA, WORK, K)
                    IF (FORMA .EQ. 'Q') THEN
                        CALL DLA_QTRMM2('L','T', K, N, ALPHA, A, LDA, B, LDB, WORK, K)
                    END IF
                    DO J=1,N
                        DO I=1,K
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*K)
                        END DO
                    END DO
                END IF
                IF (M.GT.K) THEN
                    CALL DGEMM('T','N', M-K, N, K, ALPHA, A(1,K+1), LDA, B(1,1), LDB, BETA, C(K+1,1), LDC)
                END IF
            ELSE ! M< K
                !   K      N
                !  aaA   bbbb  CCCC
                !M aaAA  bbbb  CCCC
                !  aaAAA BBBB  CCCC
                !        BBBB
                !        BBBB
                ! CALL DGEMM(TRANSA, 'N', M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)

                IF (BETA.EQ.0D0) THEN
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = B(I+K-M,J)
                        END DO
                    END DO
                    CALL DTRMM('L','U','T', 'N', M, N, ALPHA, A(K-M+1,1), LDA, C, LDC)
                    IF (FORMA.eq.'Q') THEN
                        CALL DLA_QTRMM2('L','T', M, N, ALPHA, A(K-M+1,1), LDA, B(K-M+1,1), LDB, C, LDC)
                    END IF
                ELSE
                    DO J=1,N
                        DO I=1,M
                            WORK(I+(J-1)*M) = B(I+K-M,J)
                        END DO
                    END DO
                    CALL DTRMM('L','U','T', 'N', M, N, ALPHA, A(K-M+1,1), LDA, WORK, M)
                    IF (FORMA.eq.'Q') THEN
                        CALL DLA_QTRMM2('L','T', M, N, ALPHA, A(K-M+1,1), LDA, B(K-M+1,1), LDB, WORK, M)
                    END IF
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (K.GT.M) THEN
                    CALL DGEMM('T','N', M, N, K-M, ALPHA, A(1,1), LDA, B(1,1), LDB, ONE, C(1,1), LDC)
                END IF
            END IF
        END IF
    END SUBROUTINE

    SUBROUTINE DGEMM_RIGHT(FORMB, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, WORK)
        IMPLICIT NONE

        INTEGER M,N,K,LDA,LDB,LDC,I,J
        DOUBLE PRECISION A(LDA,*), B(LDB, *), C(LDC,*)
        DOUBLE PRECISION WORK(*)
        DOUBLE PRECISION ALPHA, BETA
        CHARACTER(1) TRANSB,FORMB

        ! Local Variables
        DOUBLE PRECISION ONE
        PARAMETER (ONE = 1.0D0)

        IF (FORMB .EQ. 'F' .OR.  MAX(MAX(M,N),K) .LT. 40) THEN ! .or.MAX(MAX(M,N),K).lt.MACHINE(4)) THEN
            CALL DGEMM('N', TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
            RETURN
        ENDIF

        IF (TRANSB .EQ. 'N') THEN
            ! CALL DGEMM('N', TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
            ! RETURN

            IF (N.GE.K) THEN
                !   K    N
                !   AA BBbbb CCccc
                ! M AA  Bbbb CCccc
                !   AA       CCccc
                !   AA       CCccc
                IF (BETA.EQ.0D0) THEN
                    DO J=1,K
                        DO I=1,M
                            C(I,J) = A(I,J)
                        END DO
                    END DO
                    CALL DTRMM('R','U','N','N', M, K, ALPHA, B, LDB, C, LDC)
                    IF (FORMB .EQ. 'Q') THEN
                        CALL DLA_QTRMM2('R','N', M, K, ALPHA, B, LDB, A, LDA, C, LDC)
                    END IF
                ELSE
                    DO J=1,K
                        DO I=1,M
                            WORK(I+(J-1)*M) = A(I,J)
                        END DO
                    END DO
                    CALL DTRMM('R','U','N','N', M, K, ALPHA, B, LDB, WORK, M)
                    IF (FORMB .EQ.'Q') THEN
                        CALL DLA_QTRMM2('R','N', M, K, ALPHA, B, LDB, A, LDA, WORK, M)
                    END IF
                    DO J=1,K
                        DO I=1,M
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (N.GT.K) THEN
                    CALL DGEMM('N','N', M, N-K, K, ALPHA, A(1,1), LDA, B(1,K+1), LDB, BETA, C(1,K+1), LDC)
                END IF
            ELSE
                !   K    N
                !   aAAA bbb CCC
                ! M aAAA BBB CCC
                !   aAAA  BB CCC
                !   aAAA   B CCC
                !   aAAA     CCC
                IF (BETA.EQ.0D0) THEN
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = A(I,J+K-N)
                        END DO
                    END DO
                    CALL DTRMM('R','U','N','N', M, N, ALPHA, B(K-N+1,1), LDB, C, LDC)
                    IF (FORMB .EQ. 'Q') THEN
                        CALL DLA_QTRMM2('R','N', M, N, ALPHA, B(K-N+1,1), LDB, A(1,K-N+1), LDA, C, LDC)
                    END IF
                ELSE
                    DO J=1,N-1
                        DO I=1,M-1
                            WORK(I+(J-1)*M) = A(I,J+K-N)
                        END DO
                    END DO
                    CALL DTRMM('R','U','N','N', M, N, ALPHA, B(K-N+1,1), LDB, WORK, M)
                    IF (FORMB .EQ. 'Q') THEN
                        CALL DLA_QTRMM2('R','N', M, N, ALPHA, B(K-N+1,1), LDB, A(1,K-N+1), LDA, WORK, M)
                    END IF
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (K.GT.N) THEN
                    CALL DGEMM('N','N', M, N, K-N, ALPHA, A(1,1), LDA, B(1,1), LDB, ONE, C(1,1), LDC)
                END IF
            END IF
        ELSE ! TRANSB = B
            ! CALL DGEMM('N', TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)

            IF (K.GE.N) THEN
                !    K     N
                ! M AAaaa  B   CC
                !   AAaaa  BB  CC
                !   AAaaa  bb  CC
                !          bb
                !          bb
                IF (BETA.EQ.0D0) THEN
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = A(I,J)
                        END DO
                    END DO
                    CALL DTRMM('R','U','T', 'N', M, N, ALPHA, B, LDB, C, LDC)
                    IF (FORMB .EQ.'Q') THEN
                        CALL DLA_QTRMM2('R','T', M, N, ALPHA, B, LDB, A, LDA, C, LDC)
                    END IF
                ELSE
                    DO J=1,N
                        DO I=1,M
                            WORK(I+(J-1)*M) = A(I,J)
                        END DO
                    END DO
                    CALL DTRMM('R','U','T', 'N', M, N, ALPHA, B, LDB, WORK, M)
                    IF (FORMB .EQ. 'Q') THEN
                        CALL DLA_QTRMM2('R','T', M, N, ALPHA, B, LDB, A, LDA, WORK, M)
                    END IF
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (K.GT.N) THEN
                    CALL DGEMM('N', 'T', M, N, K-N, ALPHA, A(1,N+1), LDA, B(1,N+1), LDB, ONE, C(1,1), LDC)
                END IF
            ELSE
                !    K     N
                ! M AAA  bbB    ccCCC
                !   AAA  bbBB   ccCCC
                !   AAA  bbBBB  ccCCC
                !   AAA         ccCCC
                IF (BETA.EQ.0D0) THEN
                    DO J=1,K
                        DO I=1,M
                            C(I,J+N-K) = A(I,J)
                        END DO
                    END DO
                    CALL DTRMM('R','U','T', 'N', M, K, ALPHA, B(N-K+1,1), LDB, C(1,N-K+1), LDC)
                    IF (FORMB .EQ.'Q') THEN
                        CALL DLA_QTRMM2('R','T', M, K, ALPHA, B(N-K+1,1), LDB, A, LDA, C(1,N-K+1), LDC)
                    END IF
                ELSE
                    DO J=1,K
                        DO I=1,M
                            WORK(I+(J-1)*M) = A(I,J)
                        END DO
                    END DO
                    CALL DTRMM('R','U','T', 'N', M, K, ALPHA, B(N-K+1,1), LDB, WORK, M)
                    IF (FORMB .EQ.'Q') THEN
                        CALL DLA_QTRMM2('R','T', M, K, ALPHA, B(N-K+1,1), LDB, A, LDA, WORK, M)
                    END IF
                    DO J=1,K
                        DO I=1,M
                            C(I,J+N-K) = BETA * C(I,J+N-K) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (N.GT.K) THEN
                    CALL DGEMM('N', 'T', M, N-K, K, ALPHA, A(1,1), LDA, B(1,1), LDB, BETA, C(1,1), LDC)
                END IF
            END IF
        END IF
    END SUBROUTINE


END SUBROUTINE
