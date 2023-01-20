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

! Arguments
CHARACTER TRANSA(1), TRANSB(1)
REAL SGN1, SGN2, SCALE
INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
REAL A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF,*)
REAL WORK(*)


! Local Variables
INTEGER ININFO
INTEGER K, L, KB, LB, KH, LE, LH, KE, IT
INTEGER DIMMAT, INFO1, IPIV(8), JPIV(8)
REAL MAT(8,8), RHS(8)
REAL SCAL

LOGICAL BTRANSA, BTRANSB
REAL ZERO , ONE
PARAMETER (ZERO = 0.0, ONE = 1.0)
INTEGER IONE
PARAMETER(IONE = 1)

REAL AL(MAXNB,MAXNB), BL(MAXNB,MAXNB), CL(MAXNB,MAXNB), DL(MAXNB,MAXNB), EL(MAXNB,MAXNB), FL(MAXNB, MAXNB)
INTEGER LDAL, LDBL, LDCL, LDDL, LDEL, LDFL

!dir$ attributes align: 64:: AL
!dir$ attributes align: 64:: BL
!dir$ attributes align: 64:: CL
!dir$ attributes align: 64:: DL
!dir$ attributes align: 64:: EL
!dir$ attributes align: 64:: FL
!dir$ attributes align: 64:: MAT
!dir$ attributes align: 64:: RHS
!dir$ attributes align: 64:: IPIV
!dir$ attributes align: 64:: JPIV
!IBM* ALIGN(64,AL,BL,CL,DL,EL,FL,MAT,RHS,IPIV,JPIV)



! EXTERNAL FUNCTIONS
EXTERNAL SGER
EXTERNAL SGESC2
EXTERNAL SGETC2X
EXTERNAL SLA_SMALL_SOLVE8
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
    CALL XERROR_HANDLER('SLA_TGCSYLV_DUAL_L2_LOCAL_COPY', -INFO)
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
    ! ( N, N ) - in the standard equation. This means
    ! ( T, T ) for the dual one.

    ! Transpose is costly, skip the lower triangular part.
    DO L = 1, M
        DO K = MAX(1,L-1), M
            AL(K,L) = A(L,K)
            CL(K,L) = C(L,K)
        END DO
    END DO

    ! Due to the "good" memory access here, it is not worth to skip the lower triangular part.
    ! DO L =1, N
    !     DO K = 1, MIN(N,L+1)
    !         BL(K,L) = SGN1*B(K,L)
    !         DL(K,L) = SGN2*D(K,L)
    !     END DO
    ! END DO
    !
    BL(1:N,1:N) = SGN1*B(1:N,1:N)
    DL(1:N,1:N) = SGN2*D(1:N,1:N)
    EL(1:M, 1:N) = E(1:M,1:N)
    FL(1:M, 1:N) = F(1:M,1:N)

    L = N
    DO WHILE ( L .GT. 0 )
        LB = 1
        LH = L
        IF ( L .GT. 1 ) THEN
            IF (BL(L,L-1) .NE. ZERO) THEN
                LB = 2
                L = L - 1
            ELSE IF (DL(L,L-1) .NE. ZERO)  THEN
                LB = 2
                L = L - 1
            END IF
        END IF
        LE = L - 1


        K = 1
        DO WHILE ( K .LE. M )
            KB = 1
            IF ( K .LT. M ) THEN
                IF ( ( AL(K,K+1) .NE. ZERO )) THEN
                    KB = 2
                END IF
            END IF
            KH = K + KB - 1

            IT = 0
            DO
                ! Solve inner system
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = CL(K,K)
                    MAT(2,1) = BL(L,L)
                    MAT(2,2) = DL(L,L)
                    RHS(1) = EL(K,L)
                    RHS(2) = FL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,3) = CL(K,K)

                    MAT(2,2) = AL(K,K)
                    MAT(2,4) = CL(K,K)

                    MAT(3,1) = BL(L,L)
                    MAT(3,2) = BL(L, LH)
                    MAT(3,3) = DL(L,L)
                    MAT(3,4) = DL(L, LH)

                    MAT(4,1) = BL(LH,L)
                    MAT(4,2) = BL(LH, LH)
                    MAT(4,3) = DL(LH,L)
                    MAT(4,4) = DL(LH, LH)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(K, LH)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(K, LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,3) = CL(K,K)
                    MAT(1,4) = CL(K,KH)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,3) = CL(KH,K)
                    MAT(2,4) = CL(KH,KH)

                    MAT(3,1) = BL(L,L)
                    MAT(3,3) = DL(L,L)
                    MAT(4,2) = BL(L,L)
                    MAT(4,4) = DL(L,L)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(KH, L)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(KH, L)
                ELSE
                    DIMMAT = 8

                    MAT(1:8,1:8) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,5) = CL(K,K)
                    MAT(1,6) = CL(K,KH)
                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,5) = CL(KH,K)
                    MAT(2,6) = CL(KH,KH)

                    MAT(3,3) = AL(K,K)
                    MAT(3,4) = AL(K,KH)
                    MAT(3,7) = CL(K,K)
                    MAT(3,8) = CL(K,KH)
                    MAT(4,3) = AL(KH,K)
                    MAT(4,4) = AL(KH,KH)
                    MAT(4,7) = CL(KH,K)
                    MAT(4,8) = CL(KH,KH)

                    MAT(5,1) = BL(L,L)
                    MAT(5,3) = BL(L,LH)
                    MAT(5,5) = DL(L,L)
                    MAT(5,7) = DL(L,LH)

                    MAT(6,2) = BL(L,L)
                    MAT(6,4) = BL(L,LH)
                    MAT(6,6) = DL(L,L)
                    MAT(6,8) = DL(L,LH)

                    MAT(7,1) = BL(LH,L)
                    MAT(7,3) = BL(LH,LH)
                    MAT(7,5) = DL(LH,L)
                    MAT(7,7) = DL(LH,LH)

                    MAT(8,2) = BL(LH,L)
                    MAT(8,4) = BL(LH,LH)
                    MAT(8,6) = DL(LH,L)
                    MAT(8,8) = DL(LH,LH)

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
                    CALL SLA_SMALL_SOLVE8( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL SGETC2X(DIMMAT, MAT, 8, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL SGESC2(DIMMAT, MAT, 8, RHS, IPIV, JPIV, SCAL)
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
                    EL(KH+1:M, L) = EL(KH+1:M, L) - CL(KH+1:M,K) * FL(K,L)

                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    EL(KH+1:M, L)  = EL(KH+1:M, L) -  AL(KH+1:M,K) * EL(K,L)
                    EL(KH+1:M, LH) = EL(KH+1:M, LH) - AL(KH+1:M,K) * EL(K,LH)

                    EL(KH+1:M, L)  = EL(KH+1:M, L)  - CL(KH+1:M,K) * FL(K,L)
                    EL(KH+1:M, LH) = EL(KH+1:M, LH) - CL(KH+1:M,K) * FL(K,LH)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    EL(KH+1:M,L) = EL(KH+1:M,L) - AL(KH+1:M,K)*EL(K,L) - AL(KH+1:M,KH)*EL(KH,L)
                    EL(KH+1:M,L) = EL(KH+1:M,L) - CL(KH+1:M,K)*FL(K,L) - CL(KH+1:M,KH)*FL(KH,L)
                ELSE
                    EL(KH+1:M,L)  = EL(KH+1:M,L)  - AL(KH+1:M,K) * EL(K,L)  - AL(KH+1:M,KH) * EL(KH,L)
                    EL(KH+1:M,LH) = EL(KH+1:M,LH) - AL(KH+1:M,K) * EL(K,LH) - AL(KH+1:M,KH) * EL(KH,LH)
                    EL(KH+1:M,L)  = EL(KH+1:M,L)  - CL(KH+1:M,K) * FL(K,L)  - CL(KH+1:M,KH) * FL(KH,L)
                    EL(KH+1:M,LH) = EL(KH+1:M,LH) - CL(KH+1:M,K) * FL(K,LH) - CL(KH+1:M,KH) * FL(KH,LH)
                END IF
            END IF

            K = KH + 1
        END DO

        ! Update January 2021
        IF ( L .GT. 1 ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = 1,LE
                    FL(1:M, IT) = FL(1:M,IT) - EL(1:M,L) *BL(IT,L)
                    FL(1:M, IT) = FL(1:M,IT) - FL(1:M,L) *DL(IT,L)
                END DO
                ! CALL SGER(M, LE, -ONE, EL(1,L), IONE, BL(1,L), IONE, FL(1,1), LDFL)
                ! CALL SGER(M, LE, -ONE, FL(1,L), IONE, DL(1,L), IONE, FL(1,1), LDFL)
            ELSE
                DO IT = 1, LE
                    FL(1:M, IT) = FL(1:M,IT) - EL(1:M,L) *BL(IT,L) - EL(1:M,LH) * BL(IT,LH)
                    FL(1:M, IT) = FL(1:M,IT) - FL(1:M,L) *DL(IT,L) - FL(1:M,LH) * DL(IT,LH)
                END DO


                ! CALL SGER(M, LE, -ONE, EL(1,L),  IONE, BL(1,L),  IONE, FL(1,1), LDFL)
                ! CALL SGER(M, LE, -ONE, EL(1,LH), IONE, BL(1,LH), IONE, FL(1,1), LDFL)
                ! CALL SGER(M, LE, -ONE, FL(1,L),  IONE, DL(1,L),  IONE, FL(1,1), LDFL)
                ! CALL SGER(M, LE, -ONE, FL(1,LH), IONE, DL(1,LH), IONE, FL(1,1), LDFL)
            ENDIF
        END IF

        L = L - 1
    END DO

    E(1:M, 1:N) = EL(1:M,1:N)
    F(1:M, 1:N) = FL(1:M,1:N)





ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
    ! ( N, T ) - in the standard equation. This means
    ! ( T, N ) for the dual one.

    ! Transpose is costly, skip the lower triangular part.
    DO L = 1, M
        DO K = MAX(1,L-1), M
            AL(K,L) = A(L,K)
            CL(K,L) = C(L,K)
        END DO
    END DO

    ! Due to the "good" memory access here, it is not worth to skip the lower triangular part.
    ! DO L =1, N
    !     DO K = 1, MIN(N,L+1)
    !         BL(K,L) = SGN1*B(K,L)
    !         DL(K,L) = SGN2*D(K,L)
    !     END DO
    ! END DO
    !
    BL(1:N,1:N) = SGN1*B(1:N,1:N)
    DL(1:N,1:N) = SGN2*D(1:N,1:N)
    EL(1:M, 1:N) = E(1:M,1:N)
    FL(1:M, 1:N) = F(1:M,1:N)


    L = 1
    DO WHILE ( L .LE. N )
        LB = 1
        IF ( L .LT. N ) THEN
            IF ( BL(L+1,L) .NE. ZERO ) THEN
                LB =2
            ELSE IF ( DL(L+1,L) .NE. ZERO )  THEN
                LB = 2
            END IF
        END IF
        LH = L + LB - 1

        K = 1
        DO WHILE ( K .LE. M )
            KB = 1
            IF ( K .LT. M ) THEN
                IF ( ( AL(K,K+1) .NE. ZERO )) THEN
                    KB = 2
                END IF
            END IF
            KH = K + KB - 1

            IT = 0
            DO

                ! Solve inner system
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = CL(K,K)
                    MAT(2,1) = BL(L,L)
                    MAT(2,2) = DL(L,L)
                    RHS(1) = EL(K,L)
                    RHS(2) = FL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,3) = CL(K,K)

                    MAT(2,2) = AL(K,K)
                    MAT(2,4) = CL(K,K)

                    MAT(3,1) = BL(L,L)
                    MAT(3,2) = BL(LH, L)
                    MAT(3,3) = DL(L,L)
                    MAT(3,4) = DL(LH, L)

                    MAT(4,1) = BL(L,LH)
                    MAT(4,2) = BL(LH, LH)
                    MAT(4,3) = DL(L,LH)
                    MAT(4,4) = DL(LH, LH)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(K, LH)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(K, LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,3) = CL(K,K)
                    MAT(1,4) = CL(K, KH)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,3) = CL(KH,K)
                    MAT(2,4) = CL(KH, KH)

                    MAT(3,1) = BL(L,L)
                    MAT(3,3) = DL(L,L)
                    MAT(4,2) = BL(L,L)
                    MAT(4,4) = DL(L,L)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(KH, L)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(KH, L)
                ELSE
                    DIMMAT = 8

                    MAT(1:8,1:8) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,5) = CL(K,K)
                    MAT(1,6) = CL(K,KH)
                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,5) = CL(KH,K)
                    MAT(2,6) = CL(KH,KH)

                    MAT(3,3) = AL(K,K)
                    MAT(3,4) = AL(K,KH)
                    MAT(3,7) = CL(K,K)
                    MAT(3,8) = CL(K,KH)
                    MAT(4,3) = AL(KH,K)
                    MAT(4,4) = AL(KH,KH)
                    MAT(4,7) = CL(KH,K)
                    MAT(4,8) = CL(KH,KH)

                    MAT(5,1) =  BL(L,L)
                    MAT(5,3) =  BL(LH,L)
                    MAT(5,5) =  DL(L,L)
                    MAT(5,7) =  DL(LH,L)

                    MAT(6,2) =  BL(L,L)
                    MAT(6,4) =  BL(LH,L)
                    MAT(6,6) =  DL(L,L)
                    MAT(6,8) =  DL(LH,L)

                    MAT(7,1) =  BL(L,LH)
                    MAT(7,3) =  BL(LH,LH)
                    MAT(7,5) =  DL(L,LH)
                    MAT(7,7) =  DL(LH,LH)

                    MAT(8,2) =  BL(L,LH)
                    MAT(8,4) =  BL(LH,LH)
                    MAT(8,6) =  DL(L,LH)
                    MAT(8,8) =  DL(LH,LH)

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
                    CALL SLA_SMALL_SOLVE8( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL SGETC2X(DIMMAT, MAT, 8, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL SGESC2(DIMMAT, MAT, 8, RHS, IPIV, JPIV, SCAL)
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
                    EL(KH+1:M, L) = EL(KH+1:M, L) - CL(KH+1:M,K) * FL(K,L)

                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    EL(KH+1:M, L)  = EL(KH+1:M, L)  - AL(KH+1:M,K) * EL(K,L)
                    EL(KH+1:M, LH) = EL(KH+1:M, LH) - AL(KH+1:M,K) * EL(K,LH)

                    EL(KH+1:M, L)  = EL(KH+1:M, L)  - CL(KH+1:M,K) * FL(K,L)
                    EL(KH+1:M, LH) = EL(KH+1:M, LH) - CL(KH+1:M,K) * FL(K,LH)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    EL(KH+1:M,L) = EL(KH+1:M,L) - AL(KH+1:M,K)*EL(K,L) - AL(KH+1:M,KH)*EL(KH,L)
                    EL(KH+1:M,L) = EL(KH+1:M,L) - CL(KH+1:M,K)*FL(K,L) - CL(KH+1:M,KH)*FL(KH,L)
                ELSE
                    EL(KH+1:M,L)  = EL(KH+1:M,L)  - AL(KH+1:M,K) * EL(K,L)  - AL(KH+1:M,KH) * EL(KH,L)
                    EL(KH+1:M,LH) = EL(KH+1:M,LH) - AL(KH+1:M,K) * EL(K,LH) - AL(KH+1:M,KH) * EL(KH,LH)
                    EL(KH+1:M,L)  = EL(KH+1:M,L)  - CL(KH+1:M,K) * FL(K,L)  - CL(KH+1:M,KH) * FL(KH,L)
                    EL(KH+1:M,LH) = EL(KH+1:M,LH) - CL(KH+1:M,K) * FL(K,LH) - CL(KH+1:M,KH) * FL(KH,LH)
                END IF
            END IF

            K = KH + 1
        END DO

        IF ( LH .LT. N ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = LH+1, N
                    FL(1:M,IT) =  FL(1:M,IT) - EL(1:M,L) * BL(L,IT)
                    FL(1:M,IT) =  FL(1:M,IT) - FL(1:M,L) * DL(L,IT)

                END DO
            ELSE
                DO IT = LH+1, N
                    FL(1:M,IT) =  FL(1:M,IT) - EL(1:M,L) * BL(L,IT)  - EL(1:M, LH) * BL(LH, IT)
                    FL(1:M,IT) =  FL(1:M,IT) - FL(1:M,L) * DL(L,IT)  - FL(1:M, LH) * DL(LH, IT)
                END DO
            ENDIF

            ! IF ( LB .EQ. IONE ) THEN
            !     CALL SGER(M, N-LH, -SGN1*ONE, EL(1,L), IONE, B(L,LH+1), LDB, FL(1,LH+1), LDFL)
            !     CALL SGER(M, N-LH, -SGN2*ONE, FL(1,L), IONE, D(L,LH+1), LDD, FL(1,LH+1), LDFL)
            !
            ! ELSE
            !     CALL SGER(M, N-LH, -SGN1*ONE, EL(1,L),  IONE, B(L,LH+1),  LDB, FL(1,LH+1), LDFL)
            !     CALL SGER(M, N-LH, -SGN1*ONE, EL(1,LH), IONE, B(LH,LH+1), LDB, FL(1,LH+1), LDFL)
            !
            !     CALL SGER(M, N-LH, -SGN2*ONE, FL(1,L),  IONE, D(L,LH+1),  LDD, FL(1,LH+1), LDFL)
            !     CALL SGER(M, N-LH, -SGN2*ONE, FL(1,LH), IONE, D(LH,LH+1), LDD, FL(1,LH+1), LDFL)
            ! ENDIF
        END IF

        L = LH + 1
    END DO

    E(1:M, 1:N) = EL(1:M,1:N)
    F(1:M, 1:N) = FL(1:M,1:N)

ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
    ! ( T, N ) - in the standard equation. This means
    ! ( N, T ) for the dual one.

    AL(1:M,1:M) = A(1:M,1:M)
    CL(1:M,1:M) = C(1:M,1:M)

    BL(1:N,1:N) = SGN1*B(1:N,1:N)
    DL(1:N,1:N) = SGN2*D(1:N,1:N)
    EL(1:M, 1:N) = E(1:M,1:N)
    FL(1:M, 1:N) = F(1:M,1:N)

    L = N
    DO WHILE ( L .GT. 0 )
        LB = 1
        LH = L
        IF ( L .GT. 1 ) THEN
            IF (BL(L,L-1) .NE. ZERO) THEN
                LB = 2
                L = L - 1
            ELSE IF (DL(L,L-1) .NE. ZERO)  THEN
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
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = CL(K,K)
                    MAT(2,1) = BL(L,L)
                    MAT(2,2) = DL(L,L)
                    RHS(1) = EL(K,L)
                    RHS(2) = FL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,3) = CL(K,K)

                    MAT(2,2) = AL(K,K)
                    MAT(2,4) = CL(K,K)

                    MAT(3,1) = BL(L,L)
                    MAT(3,2) = BL(L, LH)
                    MAT(3,3) = DL(L,L)
                    MAT(3,4) = DL(L, LH)

                    MAT(4,1) = BL(LH,L)
                    MAT(4,2) = BL(LH, LH)
                    MAT(4,3) = DL(LH,L)
                    MAT(4,4) = DL(LH, LH)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(K, LH)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(K, LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,3) = CL(K,K)
                    MAT(1,4) = CL(K, KH)
                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,3) = CL(KH,K)
                    MAT(2,4) = CL(KH, KH)

                    MAT(3,1) = BL(L,L)
                    MAT(3,3) = DL(L,L)
                    MAT(4,2) = BL(L,L)
                    MAT(4,4) = DL(L,L)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(KH, L)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(KH, L)
                ELSE
                    DIMMAT = 8

                    MAT(1:8,1:8) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,5) = CL(K,K)
                    MAT(1,6) = CL(K,KH)
                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,5) = CL(KH,K)
                    MAT(2,6) = CL(KH,KH)

                    MAT(3,3) = AL(K,K)
                    MAT(3,4) = AL(K,KH)
                    MAT(3,7) = CL(K,K)
                    MAT(3,8) = CL(K,KH)
                    MAT(4,3) = AL(KH,K)
                    MAT(4,4) = AL(KH,KH)
                    MAT(4,7) = CL(KH,K)
                    MAT(4,8) = CL(KH,KH)

                    MAT(5,1) =  BL(L,L)
                    MAT(5,3) =  BL(L,LH)
                    MAT(5,5) =  DL(L,L)
                    MAT(5,7) =  DL(L,LH)

                    MAT(6,2) =  BL(L,L)
                    MAT(6,4) =  BL(L,LH)
                    MAT(6,6) =  DL(L,L)
                    MAT(6,8) =  DL(L,LH)

                    MAT(7,1) =  BL(LH,L)
                    MAT(7,3) =  BL(LH,LH)
                    MAT(7,5) =  DL(LH,L)
                    MAT(7,7) =  DL(LH,LH)

                    MAT(8,2) =  BL(LH,L)
                    MAT(8,4) =  BL(LH,LH)
                    MAT(8,6) =  DL(LH,L)
                    MAT(8,8) =  DL(LH,LH)

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
                    CALL SLA_SMALL_SOLVE8( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL SGETC2X(DIMMAT, MAT, 8, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL SGESC2(DIMMAT, MAT, 8, RHS, IPIV, JPIV, SCAL)
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
                    EL(1:K-1,L) = EL(1:K-1,L) - FL(K,L) * CL(1:K-1,K)

                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    EL(1:K-1,L)  =  EL(1:K-1,L)  - EL(K,L)  * AL(1:K-1,K)
                    EL(1:K-1,LH) =  EL(1:K-1,LH) - EL(K,LH) * AL(1:K-1,K)
                    EL(1:K-1,L)  =  EL(1:K-1,L)  - FL(K,L)  * CL(1:K-1,K)
                    EL(1:K-1,LH) =  EL(1:K-1,LH) - FL(K,LH) * CL(1:K-1,K)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    EL(1:K-1, L) = EL(1:K-1,L) - AL(1:K-1,K) * EL(K,L) - AL(1:K-1,KH) *EL(KH,L)
                    EL(1:K-1, L) = EL(1:K-1,L) - CL(1:K-1,K) * FL(K,L) - CL(1:K-1,KH) *FL(KH,L)

                ELSE
                    EL(1:K-1,L)  = EL(1:K-1,L)  - AL(1:K-1,K) * EL(K,L)  - AL(1:K-1,KH) * EL(KH,L)
                    EL(1:K-1,LH) = EL(1:K-1,LH) - AL(1:K-1,K) * EL(K,LH) - AL(1:K-1,KH) * EL(KH,LH)

                    EL(1:K-1,L)  = EL(1:K-1,L)  - CL(1:K-1,K) * FL(K,L)  - CL(1:K-1,KH) * FL(KH,L)
                    EL(1:K-1,LH) = EL(1:K-1,LH) - CL(1:K-1,K) * FL(K,LH) - CL(1:K-1,KH) * FL(KH,LH)

                END IF

            END IF
            K = K - 1
        END DO

        ! Update January 2021
        IF ( L .GT. 1 ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = 1, LE
                    FL(1:M, IT) = FL(1:M,IT) - EL(1:M,L) *BL(IT,L)
                    FL(1:M, IT) = FL(1:M,IT) - FL(1:M,L) *DL(IT,L)

                END DO

            ELSE
                DO IT = 1, LE
                    FL(1:M, IT) = FL(1:M,IT) - EL(1:M,L) *BL(IT,L) - EL(1:M,LH) * BL(IT,LH)
                    FL(1:M, IT) = FL(1:M,IT) - FL(1:M,L) *DL(IT,L) - FL(1:M,LH) * DL(IT,LH)
                END DO

            ENDIF

            ! IF ( LB .EQ. IONE ) THEN
            !     CALL SGER(M, LE, -SGN1*ONE, EL(1,L), IONE, B(1,L), IONE, FL(1,1), LDFL)
            !     CALL SGER(M, LE, -SGN2*ONE, FL(1,L), IONE, D(1,L), IONE, FL(1,1), LDFL)
            ! ELSE
            !     CALL SGER(M, LE, -SGN1*ONE, EL(1,L),  IONE, B(1,L),  IONE, FL(1,1), LDFL)
            !     CALL SGER(M, LE, -SGN1*ONE, EL(1,LH), IONE, B(1,LH), IONE, FL(1,1), LDFL)
            !     CALL SGER(M, LE, -SGN2*ONE, FL(1,L),  IONE, D(1,L),  IONE, FL(1,1), LDFL)
            !     CALL SGER(M, LE, -SGN2*ONE, FL(1,LH), IONE, D(1,LH), IONE, FL(1,1), LDFL)
            ! ENDIF
        END IF
        L = L - 1
    END DO

    E(1:M, 1:N) = EL(1:M,1:N)
    F(1:M, 1:N) = FL(1:M,1:N)

ELSE IF ( BTRANSA .AND. BTRANSB) THEN
    ! ( T, T ) - in the standard equation. This means
    ! ( N, N ) for the dual one.

    AL(1:M,1:M) = A(1:M,1:M)
    CL(1:M,1:M) = C(1:M,1:M)

    BL(1:N,1:N) = SGN1*B(1:N,1:N)
    DL(1:N,1:N) = SGN2*D(1:N,1:N)
    EL(1:M, 1:N) = E(1:M,1:N)
    FL(1:M, 1:N) = F(1:M,1:N)


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
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = CL(K,K)
                    MAT(2,1) = BL(L,L)
                    MAT(2,2) = DL(L,L)
                    RHS(1) = EL(K,L)
                    RHS(2) = FL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,3) = CL(K,K)

                    MAT(2,2) = AL(K,K)
                    MAT(2,4) = CL(K,K)

                    MAT(3,1) = BL(L,L)
                    MAT(3,2) = BL(LH, L)
                    MAT(3,3) = DL(L,L)
                    MAT(3,4) = DL(LH, L)

                    MAT(4,1) = BL(L,LH)
                    MAT(4,2) = BL(LH, LH)
                    MAT(4,3) = DL(L,LH)
                    MAT(4,4) = DL(LH, LH)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(K, LH)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(K, LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,3) = CL(K,K)
                    MAT(1,4) = CL(K, KH)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,3) = CL(KH,K)
                    MAT(2,4) = CL(KH, KH)

                    MAT(3,1) = BL(L,L)
                    MAT(3,3) = DL(L,L)
                    MAT(4,2) = BL(L,L)
                    MAT(4,4) = DL(L,L)

                    RHS(1) = EL(K, L)
                    RHS(2) = EL(KH, L)
                    RHS(3) = FL(K, L)
                    RHS(4) = FL(KH, L)
                ELSE
                    DIMMAT = 8

                    MAT(1:8,1:8) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,5) = CL(K,K)
                    MAT(1,6) = CL(K,KH)
                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,5) = CL(KH,K)
                    MAT(2,6) = CL(KH,KH)

                    MAT(3,3) = AL(K,K)
                    MAT(3,4) = AL(K,KH)
                    MAT(3,7) = CL(K,K)
                    MAT(3,8) = CL(K,KH)
                    MAT(4,3) = AL(KH,K)
                    MAT(4,4) = AL(KH,KH)
                    MAT(4,7) = CL(KH,K)
                    MAT(4,8) = CL(KH,KH)

                    MAT(5,1) = BL(L,L)
                    MAT(5,3) = BL(LH,L)
                    MAT(5,5) = DL(L,L)
                    MAT(5,7) = DL(LH,L)

                    MAT(6,2) = BL(L,L)
                    MAT(6,4) = BL(LH,L)
                    MAT(6,6) = DL(L,L)
                    MAT(6,8) = DL(LH,L)

                    MAT(7,1) = BL(L,LH)
                    MAT(7,3) = BL(LH,LH)
                    MAT(7,5) = DL(L,LH)
                    MAT(7,7) = DL(LH,LH)

                    MAT(8,2) = BL(L,LH)
                    MAT(8,4) = BL(LH,LH)
                    MAT(8,6) = DL(L,LH)
                    MAT(8,8) = DL(LH,LH)

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
                    CALL SLA_SMALL_SOLVE8( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL SGETC2X(DIMMAT, MAT, 8, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL SGESC2(DIMMAT, MAT, 8, RHS, IPIV, JPIV, SCAL)
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


            ! Update January 2021
            ! in current column
            IF ( K .GT. IONE ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    EL(1:K-1,L) = EL(1:K-1,L) - EL(K,L) * AL(1:K-1,K)
                    EL(1:K-1,L) = EL(1:K-1,L) - FL(K,L) * CL(1:K-1,K)

                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    EL(1:K-1,L)  =  EL(1:K-1,L)  - EL(K,L)  * AL(1:K-1,K)
                    EL(1:K-1,LH) =  EL(1:K-1,LH) - EL(K,LH) * AL(1:K-1,K)
                    EL(1:K-1,L)  =  EL(1:K-1,L)  - FL(K,L)  * CL(1:K-1,K)
                    EL(1:K-1,LH) =  EL(1:K-1,LH) - FL(K,LH) * CL(1:K-1,K)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    EL(1:K-1, L) = EL(1:K-1,L) - AL(1:K-1,K) * EL(K,L) - AL(1:K-1,KH) * EL(KH,L)
                    EL(1:K-1, L) = EL(1:K-1,L) - CL(1:K-1,K) * FL(K,L) - CL(1:K-1,KH) * FL(KH,L)

                ELSE
                    EL(1:K-1,L)  = EL(1:K-1,L)  - AL(1:K-1,K) * EL(K,L)  - AL(1:K-1,KH) * EL(KH,L)
                    EL(1:K-1,LH) = EL(1:K-1,LH) - AL(1:K-1,K) * EL(K,LH) - AL(1:K-1,KH) * EL(KH,LH)

                    EL(1:K-1,L)  = EL(1:K-1,L)  - CL(1:K-1,K) * FL(K,L)  - CL(1:K-1,KH) * FL(KH,L)
                    EL(1:K-1,LH) = EL(1:K-1,LH) - CL(1:K-1,K) * FL(K,LH) - CL(1:K-1,KH) * FL(KH,LH)

                END IF

            END IF

            K = K - 1

        END DO

        ! Update January 2021
        IF ( LH .LT. N ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = LH+1, N
                    FL(1:M,IT) =  FL(1:M,IT) - EL(1:M,L) * BL(L,IT)
                    FL(1:M,IT) =  FL(1:M,IT) - FL(1:M,L) * DL(L,IT)
                END DO

            ELSE
                DO IT = LH+1, N
                    FL(1:M,IT) =  FL(1:M,IT) - EL(1:M,L) * BL(L,IT)  - EL(1:M, LH) * BL(LH, IT)
                    FL(1:M,IT) =  FL(1:M,IT) - FL(1:M,L) * DL(L,IT)  - FL(1:M, LH) * DL(LH, IT)
                END DO
            ENDIF


            ! IF ( LB .EQ. IONE ) THEN
            !     CALL SGER(M, N-LH, -SGN1*ONE, EL(1,L), IONE, B(L,LH+1), LDB, FL(1,LH+1), LDFL)
            !     CALL SGER(M, N-LH, -SGN2*ONE, FL(1,L), IONE, D(L,LH+1), LDD, FL(1,LH+1), LDFL)
            !
            ! ELSE
            !     CALL SGER(M, N-LH, -SGN1*ONE, EL(1,L),  IONE, B(L,LH+1),  LDB, FL(1,LH+1), LDFL)
            !     CALL SGER(M, N-LH, -SGN1*ONE, FL(1,LH), IONE, B(LH,LH+1), LDB, FL(1,LH+1), LDFL)
            !
            !     CALL SGER(M, N-LH, -SGN2*ONE, FL(1,L),  IONE, D(L,LH+1),  LDD, FL(1,LH+1), LDFL)
            !     CALL SGER(M, N-LH, -SGN2*ONE, FL(1,LH), IONE, D(LH,LH+1), LDD, FL(1,LH+1), LDFL)
            ! ENDIF
        END IF

        L = LH + 1
    END DO

    E(1:M, 1:N) = EL(1:M,1:N)
    F(1:M, 1:N) = FL(1:M,1:N)


END IF

