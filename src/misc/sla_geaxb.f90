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
SUBROUTINE SLA_GEAXB(TRANSA, TRANSB, FORMA, FORMB, M, N, K, L, ALPHA, A, LDA, X, LDX, B, LDB, BETA, C, LDC, WORK)
    IMPLICIT NONE
    ! Input Arguments

    CHARACTER(1) TRANSA, TRANSB, FORMA, FORMB
    INTEGER M,N,K,L, LDA, LDX, LDB, LDC
    REAL ALPHA, BETA
    REAL A(LDA,*), B(LDB,*), C(LDC,*), X(LDX,*)
    REAL WORK(*)


    ! Local Variables
    REAL MULT1A, MULT1B, MULT2A, MULT2B
    REAL ONE, ZERO, TWO
    PARAMETER (ONE = 1.0, ZERO = 0.0, TWO = 2.0)

    ! External CODE
    EXTERNAL SGEMM
    EXTERNAL SLA_QTRMM2
    EXTERNAL STRMM
    EXTERNAL LSAME
    LOGICAL LSAME

    ! Check Input
    ! Skipped for Performance Reasons


    ! Quick Return
    IF (M .EQ. 0 .OR. N.EQ.0) THEN
        RETURN
    ENDIF


    IF (FORMA.eq.'F') THEN
        MULT1A = (TWO*REAL(M))*REAL(N)*REAL(K)
        MULT2B = (TWO*REAL(M))*REAL(N)*REAL(L)
    ELSE
        MULT1A = (ONE*REAL(M))*REAL(M)*REAL(K)+(TWO*REAL(M))*REAL(N-M)*REAL(K)
        MULT2B = (ONE*REAL(M))*REAL(N)*REAL(L)+(TWO*REAL(M))*REAL(N-M)*REAL(L)
    ENDIF

    IF (FORMB.eq.'F') THEN
        MULT1B = (TWO*REAL(M))*REAL(K)*REAL(L)
        MULT2A = (TWO*REAL(N))*REAL(K)*REAL(L)
    ELSE
        MULT1B = (ONE*REAL(M))*REAL(L)*REAL(L)+(TWO*REAL(M))*REAL(L)*REAL(K-L)
        MULT2A = (ONE*REAL(N))*REAL(L)*REAL(L)+(TWO*REAL(N))*REAL(L)*REAL(K-L)
    ENDIF

    IF ((MULT1A+MULT1B) .LT. (MULT2A+MULT2B)) THEN

        CALL SGEMM_LEFT (FORMA, TRANSA, M, K, N, ALPHA, A, LDA, X, LDX, ZERO, WORK, M,  WORK(M*K+1))
        CALL SGEMM_RIGHT(FORMB, TRANSB, M, L, K, ONE, WORK, M, B, LDB, BETA, C, LDC, WORK(M*K+1))

    ELSE
        CALL SGEMM_RIGHT(FORMB, TRANSB, N, L, K, ALPHA, X, LDX, B, LDB, ZERO, WORK, N, WORK(N*L+1))
        CALL SGEMM_LEFT (FORMA, TRANSA, M, L, N, ONE, A, LDA, WORK, N, BETA, C, LDC, WORK(N*L+1))
    END IF
    RETURN

CONTAINS

    SUBROUTINE SGEMM_LEFT(FORMA, TRANSA, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, WORK)
        IMPLICIT NONE

        INTEGER M,N,K,LDA,LDB,LDC
        REAL A(LDA, *), B(LDB, *) , C(LDC, *)
        REAL WORK(*)
        REAL ALPHA, BETA
        CHARACTER(1) TRANSA,FORMA

        ! Local Variables
        REAL ONE
        PARAMETER (ONE = 1.0)
        INTEGER I,J

        ! External Functions
        EXTERNAL SGEMM, STRMM, SLA_QTRMM2

        IF (FORMA .EQ.'F' .OR. MAX(MAX(M,N),K) .LT. 40) THEN ! .or.MAX(MAX(M,N),K).lt.MACHINE(4)) THEN
            CALL SGEMM(TRANSA, 'N', M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
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
                    CALL STRMM('L','U','N','N', M, N, ALPHA, A, LDA, C, LDC)
                    IF (FORMA .EQ. 'Q') THEN
                        CALL SLA_QTRMM2('L','N', M, N, ALPHA, A, LDA, B, LDB, C, LDC)
                    END IF
                ELSE
                    DO J=1,N
                        DO I=1,M
                            WORK(I+(J-1)*M) = B(I,J)
                        END DO
                    END DO
                    CALL STRMM('L','U','N','N', M, N, ALPHA, A, LDA, WORK, M)
                    IF (FORMA .EQ.'Q') THEN
                        CALL SLA_QTRMM2('L','N', M, N, ALPHA, A, LDA, B, LDB, WORK, M)
                    END IF
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (K.GT.M) THEN
                    CALL SGEMM('N','N', M, N, K-M, ALPHA, A(1,M+1), LDA, B(M+1,1), LDB, ONE, C, LDC)
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
                    CALL STRMM('L','U','N','N', K, N, ALPHA, A(M-K+1,1), LDA, C(M-K+1,1), LDC)
                    IF (FORMA .EQ. 'Q') THEN
                        CALL SLA_QTRMM2('L','N', K, N, ALPHA, A(M-K+1,1), LDA, B, LDB, C(M-K+1,1), LDC)
                    END IF
                ELSE
                    DO J=1,N
                        DO I=1,K
                            WORK(I+(J-1)*K) = B(I,J)
                        END DO
                    END DO
                    CALL STRMM('L','U','N','N', K, N, ALPHA, A(M-K+1,1), LDA, WORK, K)
                    IF (FORMA.eq.'Q') THEN
                        CALL SLA_QTRMM2('L','N', K, N, ALPHA, A(M-K+1,1), LDA, B, LDB, WORK, K)
                    END IF
                    DO J=1,N
                        DO I=1,K
                            C(I+M-K-1,J) = BETA * C(I+M-K-1,J) + WORK(I+J*K)
                        END DO
                    END DO
                END IF
                IF (M.GT.K) THEN
                    CALL SGEMM('N','N', M-K, N, K, ALPHA, A(1,1), LDA, B(1,1), LDB, BETA, C(1,1), LDC )
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
                ! CALL SGEMM(TRANSA, 'N', M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
                IF (BETA.EQ.0D0) THEN
                    DO J=1,N
                        DO I=1,K
                            C(I,J) = B(I,J)
                        END DO
                    END DO
                    CALL STRMM('L','U','T', 'N', K, N, ALPHA, A, LDA, C, LDC)
                    IF (FORMA .EQ. 'Q') THEN
                        CALL SLA_QTRMM2('L','T', K, N, ALPHA, A, LDA, B, LDB, C, LDC)
                    END IF
                ELSE
                    DO J=1,N
                        DO I=1,K
                            WORK(I+(J-1)*K) = B(I,J)
                        END DO
                    END DO
                    CALL STRMM('L','U','T', 'N', K, N, ALPHA, A, LDA, WORK, K)
                    IF (FORMA .EQ. 'Q') THEN
                        CALL SLA_QTRMM2('L','T', K, N, ALPHA, A, LDA, B, LDB, WORK, K)
                    END IF
                    DO J=1,N
                        DO I=1,K
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*K)
                        END DO
                    END DO
                END IF
                IF (M.GT.K) THEN
                    CALL SGEMM('T','N', M-K, N, K, ALPHA, A(1,K+1), LDA, B(1,1), LDB, BETA, C(K+1,1), LDC)
                END IF
            ELSE ! M< K
                !   K      N
                !  aaA   bbbb  CCCC
                !M aaAA  bbbb  CCCC
                !  aaAAA BBBB  CCCC
                !        BBBB
                !        BBBB
                ! CALL SGEMM(TRANSA, 'N', M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)

                IF (BETA.EQ.0D0) THEN
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = B(I+K-M,J)
                        END DO
                    END DO
                    CALL STRMM('L','U','T', 'N', M, N, ALPHA, A(K-M+1,1), LDA, C, LDC)
                    IF (FORMA.eq.'Q') THEN
                        CALL SLA_QTRMM2('L','T', M, N, ALPHA, A(K-M+1,1), LDA, B(K-M+1,1), LDB, C, LDC)
                    END IF
                ELSE
                    DO J=1,N
                        DO I=1,M
                            WORK(I+(J-1)*M) = B(I+K-M,J)
                        END DO
                    END DO
                    CALL STRMM('L','U','T', 'N', M, N, ALPHA, A(K-M+1,1), LDA, WORK, M)
                    IF (FORMA.eq.'Q') THEN
                        CALL SLA_QTRMM2('L','T', M, N, ALPHA, A(K-M+1,1), LDA, B(K-M+1,1), LDB, WORK, M)
                    END IF
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (K.GT.M) THEN
                    CALL SGEMM('T','N', M, N, K-M, ALPHA, A(1,1), LDA, B(1,1), LDB, ONE, C(1,1), LDC)
                END IF
            END IF
        END IF
    END SUBROUTINE

    SUBROUTINE SGEMM_RIGHT(FORMB, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, WORK)
        IMPLICIT NONE

        INTEGER M,N,K,LDA,LDB,LDC,I,J
        REAL A(LDA,*), B(LDB, *), C(LDC,*)
        REAL WORK(*)
        REAL ALPHA, BETA
        CHARACTER(1) TRANSB,FORMB

        ! Local Variables
        REAL ONE
        PARAMETER (ONE = 1.0)

        IF (FORMB .EQ. 'F' .OR.  MAX(MAX(M,N),K) .LT. 40) THEN ! .or.MAX(MAX(M,N),K).lt.MACHINE(4)) THEN
            CALL SGEMM('N', TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
            RETURN
        ENDIF

        IF (TRANSB .EQ. 'N') THEN
            ! CALL SGEMM('N', TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
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
                    CALL STRMM('R','U','N','N', M, K, ALPHA, B, LDB, C, LDC)
                    IF (FORMB .EQ. 'Q') THEN
                        CALL SLA_QTRMM2('R','N', M, K, ALPHA, B, LDB, A, LDA, C, LDC)
                    END IF
                ELSE
                    DO J=1,K
                        DO I=1,M
                            WORK(I+(J-1)*M) = A(I,J)
                        END DO
                    END DO
                    CALL STRMM('R','U','N','N', M, K, ALPHA, B, LDB, WORK, M)
                    IF (FORMB .EQ.'Q') THEN
                        CALL SLA_QTRMM2('R','N', M, K, ALPHA, B, LDB, A, LDA, WORK, M)
                    END IF
                    DO J=1,K
                        DO I=1,M
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (N.GT.K) THEN
                    CALL SGEMM('N','N', M, N-K, K, ALPHA, A(1,1), LDA, B(1,K+1), LDB, BETA, C(1,K+1), LDC)
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
                    CALL STRMM('R','U','N','N', M, N, ALPHA, B(K-N+1,1), LDB, C, LDC)
                    IF (FORMB .EQ. 'Q') THEN
                        CALL SLA_QTRMM2('R','N', M, N, ALPHA, B(K-N+1,1), LDB, A(1,K-N+1), LDA, C, LDC)
                    END IF
                ELSE
                    DO J=1,N-1
                        DO I=1,M-1
                            WORK(I+(J-1)*M) = A(I,J+K-N)
                        END DO
                    END DO
                    CALL STRMM('R','U','N','N', M, N, ALPHA, B(K-N+1,1), LDB, WORK, M)
                    IF (FORMB .EQ. 'Q') THEN
                        CALL SLA_QTRMM2('R','N', M, N, ALPHA, B(K-N+1,1), LDB, A(1,K-N+1), LDA, WORK, M)
                    END IF
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (K.GT.N) THEN
                    CALL SGEMM('N','N', M, N, K-N, ALPHA, A(1,1), LDA, B(1,1), LDB, ONE, C(1,1), LDC)
                END IF
            END IF
        ELSE ! TRANSB = B
            ! CALL SGEMM('N', TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)

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
                    CALL STRMM('R','U','T', 'N', M, N, ALPHA, B, LDB, C, LDC)
                    IF (FORMB .EQ.'Q') THEN
                        CALL SLA_QTRMM2('R','T', M, N, ALPHA, B, LDB, A, LDA, C, LDC)
                    END IF
                ELSE
                    DO J=1,N
                        DO I=1,M
                            WORK(I+(J-1)*M) = A(I,J)
                        END DO
                    END DO
                    CALL STRMM('R','U','T', 'N', M, N, ALPHA, B, LDB, WORK, M)
                    IF (FORMB .EQ. 'Q') THEN
                        CALL SLA_QTRMM2('R','T', M, N, ALPHA, B, LDB, A, LDA, WORK, M)
                    END IF
                    DO J=1,N
                        DO I=1,M
                            C(I,J) = BETA * C(I,J) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (K.GT.N) THEN
                    CALL SGEMM('N', 'T', M, N, K-N, ALPHA, A(1,N+1), LDA, B(1,N+1), LDB, ONE, C(1,1), LDC)
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
                    CALL STRMM('R','U','T', 'N', M, K, ALPHA, B(N-K+1,1), LDB, C(1,N-K+1), LDC)
                    IF (FORMB .EQ.'Q') THEN
                        CALL SLA_QTRMM2('R','T', M, K, ALPHA, B(N-K+1,1), LDB, A, LDA, C(1,N-K+1), LDC)
                    END IF
                ELSE
                    DO J=1,K
                        DO I=1,M
                            WORK(I+(J-1)*M) = A(I,J)
                        END DO
                    END DO
                    CALL STRMM('R','U','T', 'N', M, K, ALPHA, B(N-K+1,1), LDB, WORK, M)
                    IF (FORMB .EQ.'Q') THEN
                        CALL SLA_QTRMM2('R','T', M, K, ALPHA, B(N-K+1,1), LDB, A, LDA, WORK, M)
                    END IF
                    DO J=1,K
                        DO I=1,M
                            C(I,J+N-K) = BETA * C(I,J+N-K) + WORK(I+(J-1)*M)
                        END DO
                    END DO
                END IF
                IF (N.GT.K) THEN
                    CALL SGEMM('N', 'T', M, N-K, K, ALPHA, A(1,1), LDA, B(1,1), LDB, BETA, C(1,1), LDC)
                END IF
            END IF
        END IF
    END SUBROUTINE


END SUBROUTINE
