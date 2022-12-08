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
! Arguments
CHARACTER TRANSA(1), TRANSB(1)
DOUBLE PRECISION SGN1, SGN2, SCALE
INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF,*)
DOUBLE PRECISION WORK(*)


! Local Variables
    INTEGER ININFO
INTEGER K, L, KB, LB, KH, LE, LH, KE, IT, J, JP, I, IX
INTEGER DIMMAT, INFO1, IPIV(8), JPIV(8)
DOUBLE PRECISION  TEMPX(16)
DOUBLE PRECISION MAT(8,16), MAT4(4,4), RHS(8)
DOUBLE PRECISION SCAL

LOGICAL BTRANSA, BTRANSB
DOUBLE PRECISION ZERO , ONE
PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
INTEGER IONE
PARAMETER(IONE = 1)


DOUBLE PRECISION AL(MAXNB,MAXNB), BL(MAXNB,MAXNB), CL(MAXNB,MAXNB), DL(MAXNB,MAXNB), EL(MAXNB,MAXNB), FL(MAXNB, MAXNB)
INTEGER LDAL, LDBL, LDCL, LDDL, LDEL, LDFL

!dir$ attributes align: 64:: AL
!dir$ attributes align: 64:: BL
!dir$ attributes align: 64:: CL
!dir$ attributes align: 64:: DL
!dir$ attributes align: 64:: EL
!dir$ attributes align: 64:: FL
!dir$ attributes align: 64:: MAT
!dir$ attributes align: 64:: TEMPX
!dir$ attributes align: 64:: MAT4
!dir$ attributes align: 64:: RHS
!dir$ attributes align: 64:: IPIV
!dir$ attributes align: 64:: JPIV
!IBM* ALIGN(64,AL,BL,CL,DL,EL,FL, WORKA,WORKB,MAT,RHS,IPIV,JPIV)


! EXTERNAL FUNCTIONS
EXTERNAL DGER
EXTERNAL DGESC2
EXTERNAL DGETC2X
EXTERNAL DLA_SMALL_SOLVE8
EXTERNAL DLA_SMALL_SOLVE4
EXTERNAL XERROR_HANDLER
EXTERNAL LSAME
LOGICAL LSAME

LDAL =  MAXNB
LDBL =  MAXNB
LDCL = MAXNB
LDDL = MAXNB
LDEL = MAXNB
LDFL = MAXNB


! Check Input
    ININFO = INFO
INFO = 0
BTRANSA = LSAME(TRANSA, 'T')
BTRANSB = LSAME(TRANSB, 'T')
IF ( .NOT. BTRANSA .AND. .NOT. LSAME(TRANSA, 'N')) THEN
    INFO = -1
ELSE IF ( .NOT. BTRANSB .AND. .NOT. LSAME(TRANSB, 'N')) THEN
    INFO = -2
ELSE IF ( SGN1 .NE. -ONE .AND. SGN1 .NE. ONE ) THEN
    INFO = -3
ELSE IF ( SGN2 .NE. -ONE .AND. SGN2 .NE. ONE ) THEN
    INFO = -4
ELSE IF ( M .LT. 0 .OR. M .GT. MAXNB) THEN
    INFO = -5
ELSE IF ( N .LT. 0 .OR. N .GT. MAXNB) THEN
    INFO = -6
ELSE IF ( LDA .LT. MAX(1, M)) THEN
    INFO = -8
ELSE IF ( LDB .LT. MAX(1, N)) THEN
    INFO = -10
ELSE IF ( LDC .LT. MAX(1, M)) THEN
    INFO = -12
ELSE IF ( LDD .LT. MAX(1, N)) THEN
    INFO = -14
ELSE IF ( LDE .LT. MAX(1, M)) THEN
    INFO = -16
ELSE IF ( LDF .LT. MAX(1, M)) THEN
    INFO = -18
END IF

IF ( INFO .LT. 0 ) THEN
    CALL XERROR_HANDLER('DLA_TGCSYLV_L2_LOCAL_COPY', -INFO)
    RETURN
END IF


IF ( ININFO .EQ. -1) THEN
    INFO = 1
    RETURN
END IF

! Quick Return
SCALE = ONE
IF ( M .EQ. 0 .OR. N .EQ. 0 ) THEN
    RETURN
END IF

! Prepare
INFO = 0
IF( ONE.EQ.ZERO) WORK(1) = WORK(1)


IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB ) THEN
    AL(1:M, 1:M) = A(1:M, 1:M)
    CL(1:M, 1:M) = C(1:M, 1:M)
    BL(1:N, 1:N) = B(1:N, 1:N)
    DL(1:N, 1:N) = D(1:N, 1:N)
    EL(1:M, 1:N) = E(1:M, 1:N)
    FL(1:M, 1:N) = F(1:M, 1:N)

    L = 1
    DO WHILE ( L .LE. N )
        LB = 1
        IF ( L .LT. N ) THEN
            IF ( BL(L+1,L) .NE. ZERO ) THEN
                LB = 2
            ELSE IF (DL(L+1,L) .NE. ZERO )  THEN
                LB = 2
            END IF
        END IF
        LH = L + LB - 1

        K = M
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
                ! Solve inner system
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2
                    MAT4(1,1) = AL(K,K)
                    MAT4(1,2) = SGN1 * BL(L,L)
                    MAT4(2,1) = CL(K,K)
                    MAT4(2,2) = SGN2 * DL(L,L)
                    RHS(1) = EL(K,L)
                    RHS(2) = FL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 4
                    MAT4(1:4,1:4) = ZERO
                    MAT4(1,1) = AL(K,K)
                    MAT4(1,3) = SGN1 * BL(L,L)
                    MAT4(1,4) = SGN1 * BL(LH, L)

                    MAT4(2,2) = AL(K,K)
                    MAT4(2,3) = SGN1 * BL(L,LH)
                    MAT4(2,4) = SGN1 * BL(LH, LH)

                    MAT4(3,1) = CL(K,K)
                    MAT4(3,3) = SGN2 * DL(L,L)
                    MAT4(3,4) = SGN2 * DL(LH, L)

                    MAT4(4,2) = CL(K,K)
                    MAT4(4,3) = SGN2 * DL(L,LH)
                    MAT4(4,4) = SGN2 * DL(LH, LH)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(K, LH)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(K, LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 4
                    MAT4(1:4,1:4) = ZERO

                    MAT4(1,1) = AL(K,K)
                    MAT4(1,2) = AL(K,KH)
                    MAT4(1,3) = SGN1 * BL(L,L)

                    MAT4(2,1) = AL(KH,K)
                    MAT4(2,2) = AL(KH, KH)
                    MAT4(2,4) = SGN1 * BL(L,L)

                    MAT4(3,1) = CL(K,K)
                    MAT4(3,2) = CL(K,KH)
                    MAT4(3,3) = SGN2 * DL(L,L)

                    MAT4(4,2) = CL(KH, KH)
                    MAT4(4,4) = SGN2 * DL(L,L)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(KH, L)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(KH, L)
                ELSE
                    DIMMAT = 8

                    MAT(1:8,1:8) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,5) = SGN1 * BL(L,L)
                    MAT(1,7) = SGN1 * BL(LH, L)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,6) = SGN1 * BL(L,L)
                    MAT(2,8) = SGN1 * BL(LH, L)

                    MAT(3,3) = AL(K,K)
                    MAT(3,4) = AL(K,KH)
                    MAT(3,5) = SGN1 * BL(L,LH)
                    MAT(3,7) = SGN1 * BL(LH,LH)

                    MAT(4,3) = AL(KH,K)
                    MAT(4,4) = AL(KH,KH)
                    MAT(4,6) = SGN1 * BL(L,LH)
                    MAT(4,8) = SGN1 * BL(LH,LH)

                    MAT(5,1) = CL(K,K)
                    MAT(5,2) = CL(K,KH)
                    MAT(5,5) = SGN2 * DL(L,L)
                    MAT(5,7) = SGN2 * DL(LH, L)

                    MAT(6,2) = CL(KH,KH)
                    MAT(6,6) = SGN2 * DL(L,L)
                    MAT(6,8) = SGN2 * DL(LH, L)

                    MAT(7,3) = CL(K,K)
                    MAT(7,4) = CL(K,KH)
                    MAT(7,5) = SGN2 * DL(L,LH)
                    MAT(7,7) = SGN2 * DL(LH,LH)

                    MAT(8,4) = CL(KH,KH)
                    MAT(8,6) = SGN2 * DL(L,LH)
                    MAT(8,8) = SGN2 * DL(LH,LH)


                    MAT(1,9) = EL(K,L)
                    MAT(2,9) = EL(KH,L)
                    MAT(3,9) = EL(K,LH)
                    MAT(4,9) = EL(KH,LH)

                    MAT(5,9) = FL(K,L)
                    MAT(6,9) = FL(KH,L)
                    MAT(7,9) = FL(K,LH)
                    MAT(8,9) = FL(KH,LH)
                END IF

                IF ( IT .EQ. 0 ) THEN
                    IF ( DIMMAT .EQ. 8 ) THEN
                        ! CALL DLA_SMALL_SOLVE8( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                        INCLUDE 'inc_solve8.f90'
                    ELSE
                        CALL DLA_SMALL_SOLVE4( DIMMAT, MAT4, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    END IF
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    IF (DIMMAT .EQ. 8) THEN
                        CALL DGETC2X(DIMMAT, MAT, 8, IPIV, JPIV, INFO1)
                        IF ( INFO1 .NE. 0 ) INFO = INFO1
                        CALL DGESC2(DIMMAT, MAT, 8, MAT(1,9), IPIV, JPIV, SCAL)
                    ELSE
                        CALL DGETC2X(DIMMAT, MAT4, 4, IPIV, JPIV, INFO1)
                        IF ( INFO1 .NE. 0 ) INFO = INFO1
                        CALL DGESC2(DIMMAT, MAT4, 4, RHS, IPIV, JPIV, SCAL)

                    END IF
                    IF ( SCAL .NE. ONE ) THEN
                        EL(1:M, 1:N) = SCAL*EL(1:M, 1:N)
                        FL(1:M, 1:N) = SCAL*FL(1:M, 1:N)

                        SCALE = SCALE * SCAL
                        INFO = 1
                    END IF
                    EXIT
                END IF

            END DO

            IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                EL(K,L) = RHS(1)
                FL(K,L) = RHS(2)

            ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN

                EL(K, L)  = RHS(1)
                EL(K, LH)  = RHS(2)
                FL(K, L)  = RHS(3)
                FL(K, LH)  = RHS(4)

            ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                EL(K, L)  = RHS(1)
                EL(KH, L)  = RHS(2)
                FL(K, L)  = RHS(3)
                FL(KH, L)  = RHS(4)
            ELSE
                EL(K,L)   = MAT(1,9)
                EL(KH,L)  = MAT(2,9)
                EL(K,LH)  = MAT(3,9)
                EL(KH,LH) = MAT(4,9)

                FL(K,L)   = MAT(5,9)
                FL(KH,L)  = MAT(6,9)
                FL(K,LH)  = MAT(7,9)
                FL(KH,LH) = MAT(8,9)
            END IF


            ! Update January 2021
            ! in current row
            IF ( K .GT. IONE ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    EL(1:K-1,L) = EL(1:K-1,L) - EL(K,L) * AL(1:K-1,K)
                    FL(1:K-1,L) = FL(1:K-1,L) - EL(K,L) * CL(1:K-1,K)

                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    EL(1:K-1,L) =  EL(1:K-1,L) - EL(K,L) * AL(1:K-1,K)
                    EL(1:K-1,LH) =  EL(1:K-1,LH) - EL(K,LH) * AL(1:K-1,K)
                    FL(1:K-1,L) =  FL(1:K-1,L) - EL(K,L) * CL(1:K-1,K)
                    FL(1:K-1,LH) =  FL(1:K-1,LH) - EL(K,LH) * CL(1:K-1,K)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    EL(1:K-1, L) = EL(1:K-1,L) - AL(1:K-1,K) * EL(K,L) - AL(1:K-1,KH) *EL(KH,L)
                    FL(1:K-1, L) = FL(1:K-1,L) - CL(1:K-1,K) * EL(K,L) - CL(1:K-1,KH) *EL(KH,L)

                ELSE
                    EL(1:K-1,L) = EL(1:K-1,L) - AL(1:K-1,K) * EL(K,L) - AL(1:K-1,KH) * EL(KH,L)
                    EL(1:K-1,LH) = EL(1:K-1,LH) - AL(1:K-1,K) * EL(K,LH) - AL(1:K-1,KH) * EL(KH,LH)

                    FL(1:K-1,L) = FL(1:K-1,L) - CL(1:K-1,K) * EL(K,L) - CL(1:K-1,KH) * EL(KH,L)
                    FL(1:K-1,LH) = FL(1:K-1,LH) - CL(1:K-1,K) * EL(K,LH) - CL(1:K-1,KH) * EL(KH,LH)

                END IF

            END IF

            K = K - 1

        END DO

        ! Update January 2021
        IF ( LH .LT. N ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = LH+1, N
                    EL(1:M,IT) =  EL(1:M,IT) - FL(1:M,L) * (BL(L,IT) * SGN1)
                END DO

                DO IT = LH+1, N
                    FL(1:M,IT) =  FL(1:M,IT) - FL(1:M,L) * (DL(L,IT) * SGN2)
                END DO

            ELSE
                DO IT = LH+1, N
                    EL(1:M,IT) =  EL(1:M,IT) - FL(1:M,L) * (BL(L,IT)*SGN1)  - FL(1:M, LH) * (BL(LH, IT)*SGN1)
                END DO
                DO IT = LH+1, N
                    FL(1:M,IT) =  FL(1:M,IT) - FL(1:M,L) * (DL(L,IT)*SGN2)  - FL(1:M, LH) * (DL(LH, IT)*SGN2)
                END DO
            ENDIF
        END IF

        L = LH + 1
    END DO

    E(1:M, 1:N) = EL(1:M, 1:N)
    F(1:M, 1:N) = FL(1:M, 1:N)

ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
    ! (N,T)
    AL(1:M, 1:M) = A(1:M, 1:M)
    CL(1:M, 1:M) = C(1:M, 1:M)
    BL(1:N, 1:N) = B(1:N, 1:N)
    DL(1:N, 1:N) = D(1:N, 1:N)
    EL(1:M, 1:N) = E(1:M, 1:N)
    FL(1:M, 1:N) = F(1:M, 1:N)

    L = N
    DO WHILE ( L .GT. 0 )
        LB = 1
        LH = L
        IF ( L .GT. 1 ) THEN
            IF (B(L,L-1) .NE. ZERO) THEN
                LB = 2
                L = L - 1
            ELSE IF (D(L,L-1) .NE. ZERO)  THEN
                LB = 2
                L = L - 1
            END IF
        END IF
        LE = L - 1

        K = M
        DO WHILE (K .GT. 0)
            KB = 1
            KH = K
            IF ( K .GT. 1 ) THEN
                IF ( A(K,K-1) .NE. ZERO ) THEN
                    KB = 2
                    K = K - 1
                END IF
            END IF
            KE = K - 1

            IT = 0
            DO

                ! Solve Inner System
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = SGN1 * BL(L,L)
                    MAT(2,1) = CL(K,K)
                    MAT(2,2) = SGN2 * DL(L,L)
                    RHS(1) = EL(K,L)
                    RHS(2) = FL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,3) = SGN1 * BL(L,L)
                    MAT(1,4) = SGN1 * BL(L, LH)

                    MAT(2,2) = AL(K,K)
                    MAT(2,3) = SGN1 * BL(LH,L)
                    MAT(2,4) = SGN1 * BL(LH, LH)

                    MAT(3,1) = CL(K,K)
                    MAT(3,3) = SGN2 * DL(L,L)
                    MAT(3,4) = SGN2 * DL(L, LH)

                    MAT(4,2) = CL(K,K)
                    MAT(4,3) = SGN2 * DL(LH, L)
                    MAT(4,4) = SGN2 * DL(LH, LH)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(K, LH)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(K, LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,3) = SGN1 * BL(L,L)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH, KH)
                    MAT(2,4) = SGN1 * BL(L,L)

                    MAT(3,1) = CL(K,K)
                    MAT(3,2) = CL(K,KH)
                    MAT(3,3) = SGN2 * DL(L,L)

                    MAT(4,2) = CL(KH, KH)
                    MAT(4,4) = SGN2 * DL(L,L)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(KH, L)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(KH, L)
                ELSE
                    DIMMAT = 8
                    MAT(1:8,1:8) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,5) = SGN1 * BL(L,L)
                    MAT(1,7) = SGN1 * BL(L, LH)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,6) = SGN1 * BL(L,L)
                    MAT(2,8) = SGN1 * BL(L, LH)

                    MAT(3,3) = AL(K,K)
                    MAT(3,4) = AL(K,KH)
                    MAT(3,5) = SGN1 * BL(LH,L)
                    MAT(3,7) = SGN1 * BL(LH,LH)

                    MAT(4,3) = AL(KH,K)
                    MAT(4,4) = AL(KH,KH)
                    MAT(4,6) = SGN1 * BL(LH,L)
                    MAT(4,8) = SGN1 * BL(LH,LH)

                    MAT(5,1) = CL(K,K)
                    MAT(5,2) = CL(K,KH)
                    MAT(5,5) = SGN2 * DL(L,L)
                    MAT(5,7) = SGN2 * DL(L, LH)

                    MAT(6,2) = CL(KH,KH)
                    MAT(6,6) = SGN2 * DL(L,L)
                    MAT(6,8) = SGN2 * DL(L, LH)

                    MAT(7,3) = CL(K,K)
                    MAT(7,4) = CL(K,KH)
                    MAT(7,5) = SGN2 * DL(LH,L)
                    MAT(7,7) = SGN2 * DL(LH,LH)

                    MAT(8,4) = CL(KH,KH)
                    MAT(8,6) = SGN2 * DL(LH,L)
                    MAT(8,8) = SGN2 * DL(LH,LH)


                    RHS(1) = EL(K,L)
                    RHS(2) = EL(KH,L)
                    RHS(3) = EL(K,LH)
                    RHS(4) = EL(KH,LH)

                    RHS(5) = FL(K,L)
                    RHS(6) = FL(KH,L)
                    RHS(7) = FL(K,LH)
                    RHS(8) = FL(KH,LH)

                END IF
                IF ( IT .EQ. 0 ) THEN
                    CALL DLA_SMALL_SOLVE8( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL DGETC2X(DIMMAT, MAT, 8, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL DGESC2(DIMMAT, MAT, 8, RHS, IPIV, JPIV, SCAL)
                    IF ( SCAL .NE. ONE ) THEN
                        EL(1:M, 1:N) = SCAL*EL(1:M, 1:N)
                        FL(1:M, 1:N) = SCAL*FL(1:M, 1:N)
                        SCALE = SCALE * SCAL
                        INFO = 1
                    END IF
                    EXIT
                END IF

            END DO


            IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                EL(K,L) = RHS(1)
                FL(K,L) = RHS(2)

            ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN

                EL(K, L)  = RHS(1)
                EL(K, LH)  = RHS(2)
                FL(K, L)  = RHS(3)
                FL(K, LH)  = RHS(4)

            ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                EL(K, L)  = RHS(1)
                EL(KH, L)  = RHS(2)
                FL(K, L)  = RHS(3)
                FL(KH, L)  = RHS(4)
            ELSE
                EL(K,L) = RHS(1)
                EL(KH,L)  = RHS(2)
                EL(K,LH)  = RHS(3)
                EL(KH,LH)  = RHS(4)

                FL(K,L) = RHS(5)
                FL(KH,L)  = RHS(6)
                FL(K,LH)  = RHS(7)
                FL(KH,LH)  = RHS(8)
            END IF

            IF ( K .GT. IONE ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    EL(1:K-1,L) = EL(1:K-1,L) - EL(K,L) * AL(1:K-1,K)
                    FL(1:K-1,L) = FL(1:K-1,L) - EL(K,L) * CL(1:K-1,K)

                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    EL(1:K-1,L) =  EL(1:K-1,L) - EL(K,L) * AL(1:K-1,K)
                    EL(1:K-1,LH) =  EL(1:K-1,LH) - EL(K,LH) * AL(1:K-1,K)
                    FL(1:K-1,L) =  FL(1:K-1,L) - EL(K,L) * CL(1:K-1,K)
                    FL(1:K-1,LH) =  FL(1:K-1,LH) - EL(K,LH) * CL(1:K-1,K)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    EL(1:K-1, L) = EL(1:K-1,L) - AL(1:K-1,K) * EL(K,L) - AL(1:K-1,KH) *EL(KH,L)
                    FL(1:K-1, L) = FL(1:K-1,L) - CL(1:K-1,K) * EL(K,L) - CL(1:K-1,KH) *EL(KH,L)

                ELSE
                    EL(1:K-1,L) = EL(1:K-1,L) - AL(1:K-1,K) * EL(K,L) - AL(1:K-1,KH) * EL(KH,L)
                    EL(1:K-1,LH) = EL(1:K-1,LH) - AL(1:K-1,K) * EL(K,LH) - AL(1:K-1,KH) * EL(KH,LH)

                    FL(1:K-1,L) = FL(1:K-1,L) - CL(1:K-1,K) * EL(K,L) - CL(1:K-1,KH) * EL(KH,L)
                    FL(1:K-1,LH) = FL(1:K-1,LH) - CL(1:K-1,K) * EL(K,LH) - CL(1:K-1,KH) * EL(KH,LH)

                END IF

            END IF
            K = K - 1
        END DO

        ! Update January 2021
        IF ( L .GT. 1 ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = 1, LE
                    EL(1:M, IT) = EL(1:M,IT) - SGN1*FL(1:M,L) *BL(IT,L)
                    FL(1:M, IT) = FL(1:M,IT) - SGN2*FL(1:M,L) *DL(IT,L)

                END DO

            ELSE
                DO IT = 1, LE
                    EL(1:M, IT) = EL(1:M,IT) - SGN1*FL(1:M,L) *BL(IT,L) -SGN1 *FL(1:M,LH) * BL(IT,LH)
                    FL(1:M, IT) = FL(1:M,IT) - SGN2*FL(1:M,L) *DL(IT,L) -SGN2 *FL(1:M,LH) * DL(IT,LH)
                END DO

            ENDIF
        END IF
        L = L - 1
    END DO

    E(1:M, 1:N) = EL(1:M, 1:N)
    F(1:M, 1:N) = FL(1:M, 1:N)

ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
    ! (T, N)
    AL(1:M, 1:M) = TRANSPOSE(A(1:M,1:M))
    CL(1:M, 1:M) = TRANSPOSE(C(1:M,1:M))
    BL(1:N, 1:N) = B(1:N, 1:N)
    DL(1:N, 1:N) = D(1:N, 1:N)
    EL(1:M, 1:N) = E(1:M, 1:N)
    FL(1:M, 1:N) = F(1:M, 1:N)

    L = 1
    DO WHILE ( L .LE. N )
        LB = 1
        IF ( L .LT. N ) THEN
            IF ( B(L+1,L) .NE. ZERO ) THEN
                LB =2
            ELSE IF ( D(L+1,L) .NE. ZERO )  THEN
                LB = 2
            END IF
        END IF
        LH = L + LB - 1

        K = 1
        DO WHILE ( K .LE. M )
            KB = 1
            IF ( K .LT. M ) THEN
                IF ( ( A(K+1,K) .NE. ZERO )) THEN
                    KB = 2
                END IF
            END IF
            KH = K + KB - 1

            IT = 0
            DO

                ! Solve Inner System
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = SGN1 * BL(L,L)
                    MAT(2,1) = CL(K,K)
                    MAT(2,2) = SGN2 * DL(L,L)
                    RHS(1) = EL(K,L)
                    RHS(2) = FL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,3) = SGN1 * BL(L,L)
                    MAT(1,4) = SGN1 * BL(LH, L)

                    MAT(2,2) = AL(K,K)
                    MAT(2,3) = SGN1 * BL(L,LH)
                    MAT(2,4) = SGN1 * BL(LH, LH)

                    MAT(3,1) = CL(K,K)
                    MAT(3,3) = SGN2 * DL(L,L)
                    MAT(3,4) = SGN2 * DL(LH, L)

                    MAT(4,2) = CL(K,K)
                    MAT(4,3) = SGN2 * DL(L,LH)
                    MAT(4,4) = SGN2 * DL(LH, LH)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(K, LH)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(K, LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,3) = SGN1 * BL(L,L)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = A(KH, KH)
                    MAT(2,4) = SGN1 * BL(L,L)

                    MAT(3,1) = CL(K,K)
                    MAT(3,3) = SGN2 * DL(L,L)

                    MAT(4,1) = CL(K, KH)
                    MAT(4,2) = CL(KH, KH)
                    MAT(4,4) = SGN2 * DL(L,L)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(KH, L)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(KH, L)

                ELSE
                    DIMMAT = 8
                    MAT(1:8,1:8) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,5) = SGN1 * BL(L,L)
                    MAT(1,7) = SGN1 * BL(LH, L)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,6) = SGN1 * BL(L,L)
                    MAT(2,8) = SGN1 * BL(LH, L)

                    MAT(3,3) = AL(K,K)
                    MAT(3,4) = AL(K,KH)
                    MAT(3,5) = SGN1 * BL(L,LH)
                    MAT(3,7) = SGN1 * BL(LH,LH)

                    MAT(4,3) = AL(KH,K)
                    MAT(4,4) = AL(KH,KH)
                    MAT(4,6) = SGN1 * BL(L,LH)
                    MAT(4,8) = SGN1 * BL(LH,LH)

                    MAT(5,1) = CL(K,K)
                    MAT(5,5) = SGN2 * DL(L,L)
                    MAT(5,7) = SGN2 * DL(LH, L)

                    MAT(6,1) = CL(KH,K)
                    MAT(6,2) = CL(KH,KH)
                    MAT(6,6) = SGN2 * DL(L,L)
                    MAT(6,8) = SGN2 * DL(LH, L)

                    MAT(7,3) = CL(K,K)
                    MAT(7,5) = SGN2 * DL(L,LH)
                    MAT(7,7) = SGN2 * DL(LH,LH)

                    MAT(8,3) = CL(KH,K)
                    MAT(8,4) = CL(KH,KH)
                    MAT(8,6) = SGN2 * DL(L,LH)
                    MAT(8,8) = SGN2 * DL(LH,LH)


                    RHS(1) = EL(K,L)
                    RHS(2) = EL(KH,L)
                    RHS(3) = EL(K,LH)
                    RHS(4) = EL(KH,LH)

                    RHS(5) = FL(K,L)
                    RHS(6) = FL(KH,L)
                    RHS(7) = FL(K,LH)
                    RHS(8) = FL(KH,LH)

                END IF
                IF ( IT .EQ. 0 ) THEN
                    CALL DLA_SMALL_SOLVE8( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL DGETC2X(DIMMAT, MAT, 8, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL DGESC2(DIMMAT, MAT, 8, RHS, IPIV, JPIV, SCAL)
                    IF ( SCAL .NE. ONE ) THEN
                        EL(1:M, 1:N) = SCAL*EL(1:M, 1:N)
                        FL(1:M, 1:N) = SCAL*FL(1:M, 1:N)
                        SCALE = SCALE * SCAL
                        INFO = 1
                    END IF
                    EXIT
                END IF

            END DO

            IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                EL(K,L) = RHS(1)
                FL(K,L) = RHS(2)

            ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN

                EL(K, L)  = RHS(1)
                EL(K, LH)  = RHS(2)
                FL(K, L)  = RHS(3)
                FL(K, LH)  = RHS(4)

            ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                EL(K, L)  = RHS(1)
                EL(KH, L)  = RHS(2)
                FL(K, L)  = RHS(3)
                FL(KH, L)  = RHS(4)
            ELSE
                EL(K,L) = RHS(1)
                EL(KH,L)  = RHS(2)
                EL(K,LH)  = RHS(3)
                EL(KH,LH)  = RHS(4)

                FL(K,L) = RHS(5)
                FL(KH,L)  = RHS(6)
                FL(K,LH)  = RHS(7)
                FL(KH,LH)  = RHS(8)
            END IF


            ! in current row

            IF ( KH .LT. M ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    EL(KH+1:M, L) = EL(KH+1:M, L) - AL(KH+1:M,K) * EL(K,L)
                    FL(KH+1:M, L) = FL(KH+1:M, L) - CL(KH+1:M,K) * EL(K,L)

                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    EL(KH+1:M, L) = EL(KH+1:M, L) - AL(KH+1:M,K) * EL(K,L)
                    EL(KH+1:M, LH) = EL(KH+1:M, LH) - AL(KH+1:M,K) * EL(K,LH)

                    FL(KH+1:M, L) = FL(KH+1:M, L) - CL(KH+1:M,K) * EL(K,L)
                    FL(KH+1:M, LH) = FL(KH+1:M, LH) - CL(KH+1:M,K) * EL(K,LH)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    EL(KH+1:M,L) = EL(KH+1:M,L) - AL(KH+1:M,K)*EL(K,L) - AL(KH+1:M,KH)*EL(KH,L)
                    FL(KH+1:M,L) = FL(KH+1:M,L) - CL(KH+1:M,K)*EL(K,L) - CL(KH+1:M,KH)*EL(KH,L)
                ELSE
                    EL(KH+1:M,L) = EL(KH+1:M,L) - AL(KH+1:M,K) * EL(K,L) - AL(KH+1:M,KH) * EL(KH,L)
                    EL(KH+1:M,LH) = EL(KH+1:M,LH) - AL(KH+1:M,K) * EL(K,LH) - AL(KH+1:M,KH) * EL(KH,LH)
                    FL(KH+1:M,L) = FL(KH+1:M,L) - CL(KH+1:M,K) * EL(K,L) - CL(KH+1:M,KH) * EL(KH,L)
                    FL(KH+1:M,LH) = FL(KH+1:M,LH) - CL(KH+1:M,K) * EL(K,LH) - CL(KH+1:M,KH) * EL(KH,LH)
                END IF
            END IF

            K = KH + 1
        END DO

        ! Update January 2021
        IF ( LH .LT. N ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = LH+1, N
                    EL(1:M,IT) =  EL(1:M,IT) - SGN1*FL(1:M,L) * BL(L,IT)
                    FL(1:M,IT) =  FL(1:M,IT) - SGN2*FL(1:M,L) * DL(L,IT)

                END DO
            ELSE
                DO IT = LH+1, N
                    EL(1:M,IT) =  EL(1:M,IT) - SGN1*FL(1:M,L) * BL(L,IT)  - SGN1*FL(1:M, LH) * B(LH, IT)
                    FL(1:M,IT) =  FL(1:M,IT) - SGN2*FL(1:M,L) * DL(L,IT)  - SGN2*FL(1:M, LH) * D(LH, IT)
                END DO
            ENDIF
        END IF

        L = LH + 1
    END DO


    E(1:M, 1:N) = EL(1:M, 1:N)
    F(1:M, 1:N) = FL(1:M, 1:N)


ELSE IF ( BTRANSA .AND. BTRANSB) THEN
    ! (T, T)
    AL(1:M, 1:M) = TRANSPOSE(A(1:M,1:M))
    CL(1:M, 1:M) = TRANSPOSE(C(1:M,1:M))
    BL(1:N, 1:N) = B(1:N, 1:N)
    DL(1:N, 1:N) = D(1:N, 1:N)
    EL(1:M, 1:N) = E(1:M, 1:N)
    FL(1:M, 1:N) = F(1:M, 1:N)

    L = N
    DO WHILE ( L .GT. 0 )
        LB = 1
        LH = L
        IF ( L .GT. 1 ) THEN
            IF (B(L,L-1) .NE. ZERO) THEN
                LB = 2
                L = L - 1
            ELSE IF (D(L,L-1) .NE. ZERO)  THEN
                LB = 2
                L = L - 1
            END IF
        END IF
        LE = L - 1


        K = 1
        DO WHILE ( K .LE. M )
            KB = 1
            IF ( K .LT. M ) THEN
                IF ( ( A(K+1,K) .NE. ZERO )) THEN
                    KB = 2
                END IF
            END IF
            KH = K + KB - 1

            IT = 0
            DO
                ! Solve Inner System
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = SGN1 * BL(L,L)
                    MAT(2,1) = CL(K,K)
                    MAT(2,2) = SGN2 * DL(L,L)
                    RHS(1) = EL(K,L)
                    RHS(2) = FL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,3) = SGN1 * BL(L,L)
                    MAT(1,4) = SGN1 * BL(L, LH)

                    MAT(2,2) = AL(K,K)
                    MAT(2,3) = SGN1 * BL(LH,L)
                    MAT(2,4) = SGN1 * BL(LH, LH)

                    MAT(3,1) = CL(K,K)
                    MAT(3,3) = SGN2 * DL(L,L)
                    MAT(3,4) = SGN2 * DL(L, LH)

                    MAT(4,2) = CL(K,K)
                    MAT(4,3) = SGN2 * DL(LH,L)
                    MAT(4,4) = SGN2 * DL(LH, LH)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(K, LH)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(K, LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,3) = SGN1 * BL(L,L)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = A(KH, KH)
                    MAT(2,4) = SGN1 * BL(L,L)

                    MAT(3,1) = CL(K,K)
                    MAT(3,3) = SGN2 * DL(L,L)

                    MAT(4,1) = C(K, KH)
                    MAT(4,2) = C(KH, KH)
                    MAT(4,4) = SGN2 * DL(L,L)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(KH, L)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(KH, L)

                ELSE
                    DIMMAT = 8
                    MAT(1:8,1:8) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,5) = SGN1 * BL(L,L)
                    MAT(1,7) = SGN1 * BL(L, LH)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,6) = SGN1 * BL(L,L)
                    MAT(2,8) = SGN1 * BL(L, LH)

                    MAT(3,3) = AL(K,K)
                    MAT(3,4) = AL(K,KH)
                    MAT(3,5) = SGN1 * BL(LH,L)
                    MAT(3,7) = SGN1 * BL(LH,LH)

                    MAT(4,3) = AL(KH,K)
                    MAT(4,4) = AL(KH,KH)
                    MAT(4,6) = SGN1 * BL(LH,L)
                    MAT(4,8) = SGN1 * BL(LH,LH)

                    MAT(5,1) = CL(K,K)
                    MAT(5,5) = SGN2 * DL(L,L)
                    MAT(5,7) = SGN2 * DL(L, LH)

                    MAT(6,1) = CL(KH,K)
                    MAT(6,2) = CL(KH,KH)
                    MAT(6,6) = SGN2 * DL(L,L)
                    MAT(6,8) = SGN2 * DL(L, LH)

                    MAT(7,3) = CL(K,K)
                    MAT(7,5) = SGN2 * DL(LH,L)
                    MAT(7,7) = SGN2 * DL(LH,LH)

                    MAT(8,3) = CL(KH,K)
                    MAT(8,4) = CL(KH,KH)
                    MAT(8,6) = SGN2 * DL(LH,L)
                    MAT(8,8) = SGN2 * DL(LH,LH)


                    RHS(1) = EL(K,L)
                    RHS(2) = EL(KH,L)
                    RHS(3) = EL(K,LH)
                    RHS(4) = EL(KH,LH)

                    RHS(5) = FL(K,L)
                    RHS(6) = FL(KH,L)
                    RHS(7) = FL(K,LH)
                    RHS(8) = FL(KH,LH)

                END IF

                IF ( IT .EQ. 0 ) THEN
                    CALL DLA_SMALL_SOLVE8( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL DGETC2X(DIMMAT, MAT, 8, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL DGESC2(DIMMAT, MAT, 8, RHS, IPIV, JPIV, SCAL)
                    IF ( SCAL .NE. ONE ) THEN
                        EL(1:M, 1:N) = SCAL*EL(1:M, 1:N)
                        FL(1:M, 1:N) = SCAL*FL(1:M, 1:N)

                        SCALE = SCALE * SCAL
                        INFO = 1
                    END IF
                    EXIT
                END IF

            END DO

            IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                EL(K,L) = RHS(1)
                FL(K,L) = RHS(2)

            ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN

                EL(K, L)  = RHS(1)
                EL(K, LH)  = RHS(2)
                FL(K, L)  = RHS(3)
                FL(K, LH)  = RHS(4)

            ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                EL(K, L)  = RHS(1)
                EL(KH, L)  = RHS(2)
                FL(K, L)  = RHS(3)
                FL(KH, L)  = RHS(4)
            ELSE
                EL(K,L) = RHS(1)
                EL(KH,L)  = RHS(2)
                EL(K,LH)  = RHS(3)
                EL(KH,LH)  = RHS(4)

                FL(K,L) = RHS(5)
                FL(KH,L)  = RHS(6)
                FL(K,LH)  = RHS(7)
                FL(KH,LH)  = RHS(8)
            END IF

            IF ( KH .LT. M ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    EL(KH+1:M, L) = EL(KH+1:M, L) - AL(KH+1:M,K) * EL(K,L)
                    FL(KH+1:M, L) = FL(KH+1:M, L) - CL(KH+1:M,K) * EL(K,L)

                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    EL(KH+1:M, L) = EL(KH+1:M, L) - AL(KH+1:M,K) * EL(K,L)
                    EL(KH+1:M, LH) = EL(KH+1:M, LH) - AL(KH+1:M,K) * EL(K,LH)

                    FL(KH+1:M, L) = FL(KH+1:M, L) - CL(KH+1:M,K) * EL(K,L)
                    FL(KH+1:M, LH) = FL(KH+1:M, LH) - CL(KH+1:M,K) * EL(K,LH)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    EL(KH+1:M,L) = EL(KH+1:M,L) - AL(KH+1:M,K)*EL(K,L) - AL(KH+1:M,KH)*EL(KH,L)
                    FL(KH+1:M,L) = FL(KH+1:M,L) - CL(KH+1:M,K)*EL(K,L) - CL(KH+1:M,KH)*EL(KH,L)
                ELSE
                    EL(KH+1:M,L) = EL(KH+1:M,L) - AL(KH+1:M,K) * EL(K,L) - AL(KH+1:M,KH) * EL(KH,L)
                    EL(KH+1:M,LH) = EL(KH+1:M,LH) - AL(KH+1:M,K) * EL(K,LH) - AL(KH+1:M,KH) * EL(KH,LH)
                    FL(KH+1:M,L) = FL(KH+1:M,L) - CL(KH+1:M,K) * EL(K,L) - CL(KH+1:M,KH) * EL(KH,L)
                    FL(KH+1:M,LH) = FL(KH+1:M,LH) - CL(KH+1:M,K) * EL(K,LH) - CL(KH+1:M,KH) * EL(KH,LH)
                END IF
            END IF

            K = KH + 1
        END DO

        ! Update January 2021
        IF ( L .GT. 1 ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = 1, LE
                    EL(1:M, IT) = EL(1:M,IT) - SGN1*FL(1:M,L) *BL(IT,L)
                    FL(1:M, IT) = FL(1:M,IT) - SGN2*FL(1:M,L) *DL(IT,L)

                END DO

            ELSE
                DO IT = 1, LE
                    EL(1:M, IT) = EL(1:M,IT) - SGN1*FL(1:M,L) *BL(IT,L) -SGN1 *FL(1:M,LH) * BL(IT,LH)
                    FL(1:M, IT) = FL(1:M,IT) - SGN2*FL(1:M,L) *DL(IT,L) -SGN2 *FL(1:M,LH) * DL(IT,LH)
                END DO

            ENDIF
        END IF

        L = L - 1
    END DO

    E(1:M, 1:N) = EL(1:M, 1:N)
    F(1:M, 1:N) = FL(1:M, 1:N)

END IF



