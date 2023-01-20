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

CHARACTER TRANSA(1), TRANSB(1)
DOUBLE PRECISION SGN, SCALE
INTEGER M, N, LDA, LDB, LDX, INFO
DOUBLE PRECISION A(LDA, *), B(LDB, *) ,  X(LDX, *)
DOUBLE PRECISION WORK(*)


! Local Variables
    INTEGER ININFO
INTEGER K, L, KB, LB, KH, LE, LH, KE, IT
INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
DOUBLE PRECISION MAT(4,4), RHS(4)
DOUBLE PRECISION A11, A12, A21, A22
DOUBLE PRECISION B11, B12, B21, B22
DOUBLE PRECISION SCAL

LOGICAL BTRANSA, BTRANSB
DOUBLE PRECISION ZERO , ONE
PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

! EXTERNAL FUNCTIONS
EXTERNAL DLA_SMALL_SOLVE4
EXTERNAL LSAME
EXTERNAL XERROR_HANDLER
EXTERNAL DGETC2X
EXTERNAL DGESC2
LOGICAL LSAME

INTEGER IONE
PARAMETER ( IONE = 1)

DOUBLE PRECISION AL(MAXNB,MAXNB), BL(MAXNB,MAXNB), XL(MAXNB,MAXNB)
INTEGER LDAL, LDBL, LDXL

!dir$ attributes align: 64:: AL
!dir$ attributes align: 64:: BL
!dir$ attributes align: 64:: XL
!dir$ attributes align: 64:: MAT
!dir$ attributes align: 64:: RHS
!dir$ attributes align: 64:: IPIV
!dir$ attributes align: 64:: JPIV
!IBM* ALIGN(64,AL,BL,XL,MAT,RHS,IPIV,JPIV)

LDAL =  MAXNB
LDBL =  MAXNB
LDXL = MAXNB


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
ELSE IF ( M .LT. 0 .OR. M .GT. MAXNB ) THEN
    INFO = -4
ELSE IF ( N .LT. 0 .OR. N .GT. MAXNB) THEN
    INFO = -5
ELSE IF ( LDA .LT. MAX(1, M)) THEN
    INFO = -7
ELSE IF ( LDB .LT. MAX(1, N)) THEN
    INFO = -9
ELSE IF ( LDX .LT. MAX(1, M)) THEN
    INFO = -11
END IF

IF ( INFO .LT. 0 ) THEN
    CALL XERROR_HANDLER('DLA_TRSYLV_L2_LOCAL_COPY_XX', -INFO)
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
WORK(1)=WORK(1)

IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB ) THEN
    AL(1:M,1:M) = A(1:M, 1:M)
    BL(1:N,1:N) = B(1:N, 1:N)
    XL(1:M, 1:N) = X(1:M, 1:N)

    L = 1
    DO WHILE ( L .LE. N )
        LB = 1
        IF ( L .LT. N ) THEN
            IF ( BL(L+1,L) .NE. ZERO ) THEN
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
                    DIMMAT = 1

                    MAT(1,1) = AL(K,K) + SGN * BL(L,L)
                    RHS(1)   = XL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 2

                    A11 = AL(K,K)

                    B11 = SGN *BL(L,L)
                    B12 = SGN *BL(L,LH)
                    B21 = SGN *BL(LH,L)
                    B22 = SGN *BL(LH,LH)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) =       B21
                    MAT(2,1) =       B12
                    MAT(2,2) = A11 + B22

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(K,LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2

                    A11 = AL(K,K)
                    A12 = AL(K,KH)
                    A21 = AL(KH,K)
                    A22 = AL(KH,KH)

                    B11 = SGN*BL(L,L)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A12
                    MAT(2,1) = A21
                    MAT(2,2) = A22 + B11

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(KH,L)

                ELSE
                    DIMMAT = 4

                    A11 = AL(K,K)
                    A12 = AL(K,KH)
                    A21 = AL(KH,K)
                    A22 = AL(KH,KH)

                    B11 = SGN* BL(L,L)
                    B12 = SGN* BL(L,LH)
                    B21 = SGN* BL(LH, L)
                    B22 = SGN* BL(LH,LH)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A12
                    MAT(1,3) =       B21
                    MAT(1,4) = ZERO

                    MAT(2,1) = A21
                    MAT(2,2) = A22 + B11
                    MAT(2,3) = ZERO
                    MAT(2,4) =       B21

                    MAT(3,1) =       B12
                    MAT(3,2) = ZERO
                    MAT(3,3) = A11 + B22
                    MAT(3,4) = A12

                    MAT(4,1) = ZERO
                    MAT(4,2) =       B12
                    MAT(4,3) = A21
                    MAT(4,4) = A22 + B22

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(KH,L)
                    RHS(3) = XL(K, LH)
                    RHS(4) = XL(KH,LH)

                END IF

                IF ( IT .EQ. 0 ) THEN
                    CALL DLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL DGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL DGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
                    IF ( SCAL .NE. ONE ) THEN
                        XL(1:M, 1:N) = SCAL*XL(1:M, 1:N)
                        SCALE = SCALE * SCAL
                        INFO = 1
                    END IF
                    EXIT
                END IF

            END DO

            IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                XL(K,L) = RHS(1)
            ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                XL(K,L) = RHS(1)
                XL(K,LH) = RHS(2)
            ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                XL(K,L) = RHS(1)
                XL(KH,L) = RHS(2)
            ELSE
                XL(K,L) = RHS(1)
                XL(KH,L) = RHS(2)
                XL(K, LH) = RHS(3)
                XL(KH,LH) = RHS(4)
            END IF

            ! Update January 2021
            ! in current row
            IF ( K .GT. IONE ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    XL(1:K-1,L) = XL(1:K-1,L) - XL(K,L) *AL(1:K-1,K)
                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    XL(1:K-1,L) = XL(1:K-1,L) - XL(K,L) *AL(1:K-1,K)
                    XL(1:K-1,LH) = XL(1:K-1,LH) - XL(K,LH) * AL(1:K-1,K)
                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    XL(1:K-1, L) = XL(1:K-1,L) - AL(1:K-1,K) *XL(K,L) - AL(1:K-1,KH) *XL(KH,L)
                ELSE
                    XL(1:K-1,L) = XL(1:K-1,L) -   XL(K,L) * AL(1:K-1,K) -  XL(KH,L) * AL(1:K-1,KH)
                    XL(1:K-1,LH) = XL(1:K-1,LH) - XL(K,LH) * AL(1:K-1,K) - XL(KH,LH) * AL(1:K-1,KH)
                END IF

            END IF

            K = K - 1

        END DO

        ! Update January 2021
        IF ( LH .LT. N ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = LH+1, N
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, L)*BL(L,IT)
                END DO
            ELSE
                DO IT = LH+1, N
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, L)*BL(L,IT)
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, LH)*BL(LH,IT)
                END DO
            ENDIF
        END IF

        L = LH + 1
    END DO

    X(1:M, 1:N) = XL(1:M, 1:N)

ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
    ! (N,T)
    AL(1:M,1:M) = A(1:M, 1:M)
    BL(1:N,1:N) = B(1:N, 1:N)
    XL(1:M, 1:N) = X(1:M, 1:N)

    L = N
    DO WHILE ( L .GT. 0 )
        LB = 1
        LH = L
        IF ( L .GT. 1 ) THEN
            IF (BL(L,L-1) .NE. ZERO) THEN
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

                ! Solve Inner System
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 1

                    MAT(1,1) = AL(K,K) + SGN * BL(L,L)
                    RHS(1)   = XL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 2

                    A11 = AL(K,K)

                    B11 = SGN*BL(L,L)
                    B21 = SGN*BL(LH,L)
                    B12 = SGN*BL(L,LH)
                    B22 = SGN*BL(LH,LH)


                    MAT(1,1) = A11 + B11
                    MAT(1,2) =       B12
                    MAT(2,1) =       B21
                    MAT(2,2) = A11 + B22

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(K,LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2

                    A11 = AL(K,K)
                    A12 = AL(K,KH)
                    A21 = AL(KH,K)
                    A22 = AL(KH,KH)

                    B11 = SGN * BL(L,L)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A12
                    MAT(2,1) = A21
                    MAT(2,2) = A22 + B11

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(KH,L)

                ELSE
                    DIMMAT = 4

                    A11 = AL(K,K)
                    A12 = AL(K,KH)
                    A21 = AL(KH,K)
                    A22 = AL(KH,KH)

                    B11 = SGN*BL(L,L)
                    B12 = SGN*BL(L,LH)
                    B21 = SGN*BL(LH, L)
                    B22 = SGN*BL(LH,LH)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A12
                    MAT(1,3) =       B12
                    MAT(1,4) = ZERO

                    MAT(2,1) = A21
                    MAT(2,2) = A22 + B11
                    MAT(2,3) = ZERO
                    MAT(2,4) =     + B12

                    MAT(3,1) = B21
                    MAT(3,2) = ZERO
                    MAT(3,3) = A11 + B22
                    MAT(3,4) = A12

                    MAT(4,1) = ZERO
                    MAT(4,2) = B21
                    MAT(4,3) = A21
                    MAT(4,4) = A22 + B22

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(KH,L)
                    RHS(3) = XL(K, LH)
                    RHS(4) = XL(KH,LH)

                END IF
                IF ( IT .EQ. 0 ) THEN
                    CALL DLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL DGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL DGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
                    IF ( SCAL .NE. ONE ) THEN
                        XL(1:M, 1:N) = SCAL*XL(1:M, 1:N)
                        SCALE = SCALE * SCAL
                        INFO = 1
                    END IF
                    EXIT
                END IF

            END DO

            IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                XL(K,L) = RHS(1)
            ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                XL(K,L) = RHS(1)
                XL(K,LH) = RHS(2)
            ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                XL(K,L) = RHS(1)
                XL(KH,L) = RHS(2)
            ELSE
                XL(K,L) = RHS(1)
                XL(KH,L) = RHS(2)
                XL(K, LH) = RHS(3)
                XL(KH,LH) = RHS(4)
            END IF

            IF ( K .GT. 1 ) THEN
                IF ( KB .EQ. IONE .AND. LB .EQ. IONE ) THEN
                    XL(1:KE, L ) = XL(1:KE, L) - XL(K,L)* AL(1:KE, K)
                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    XL(1:KE,L) = XL(1:KE,L) - XL(K,L) * AL(1:KE,K)
                    XL(1:KE,LH) = XL(1:KE,LH) - XL(K,LH) *AL(1:KE,K)
                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    XL(1:KE, L) = XL(1:KE,L) - AL(1:KE,K) *XL(K,L) - AL(1:KE,KH) *XL(KH,L)
                ELSE
                    XL(1:KE,L) = XL(1:KE,L) -   XL(K,L)  * AL(1:KE,K) - XL(KH,L) * AL(1:KE,KH)
                    XL(1:KE,LH) = XL(1:KE,LH) - XL(K,LH) * AL(1:KE,K) - XL(KH,LH) * AL(1:KE,KH)
                END IF

            END IF

            K = K - 1
        END DO

        ! Update January 2021
        IF ( L .GT. 1 ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = 1,LE
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, L)*BL(IT,L)
                END DO
            ELSE
                DO IT = 1,LE
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, L)*BL(IT,L)
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, LH)*BL(IT, LH)
                END DO
            ENDIF

        END IF
        L = L - 1
    END DO

    X(1:M, 1:N) = XL(1:M, 1:N)

ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
    ! (T, N)
    AL(1:M,1:M) = TRANSPOSE(A(1:M, 1:M))
    BL(1:N,1:N) = B(1:N, 1:N)
    XL(1:M, 1:N) = X(1:M, 1:N)

    L = 1
    DO WHILE ( L .LE. N )
        LB = 1
        IF ( L .LT. N ) THEN
            IF ( BL(L+1,L) .NE. ZERO ) THEN
                LB =2
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

                ! Solve Inner System
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 1

                    MAT(1,1) = AL(K,K) + SGN * BL(L,L)
                    RHS(1)   = XL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 2

                    A11 = AL(K,K)


                    B11 = SGN * BL(L,L)
                    B21 = SGN * BL(LH,L)
                    B12 = SGN * BL(L,LH)
                    B22 = SGN * BL(LH,LH)


                    MAT(1,1) = A11 + B11
                    MAT(1,2) =       B21
                    MAT(2,1) =       B12
                    MAT(2,2) = A11 + B22

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(K,LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2

                    A11 = AL(K,K)
                    A12 = AL(KH,K)
                    A21 = AL(K,KH)
                    A22 = AL(KH,KH)

                    B11 = SGN * BL(L,L)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A21
                    MAT(2,1) = A12
                    MAT(2,2) = A22 + B11

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(KH,L)

                ELSE
                    DIMMAT = 4

                    A11 = AL(K,K)
                    A12 = AL(KH,K)
                    A21 = AL(K,KH)
                    A22 = AL(KH,KH)


                    B11 = SGN * BL(L,L)
                    B12 = SGN * BL(L,LH)
                    B21 = SGN * BL(LH, L)
                    B22 = SGN * BL(LH,LH)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A21
                    MAT(1,3) =       B21
                    MAT(1,4) = ZERO

                    MAT(2,1) = A12
                    MAT(2,2) = A22 + B11
                    MAT(2,3) = ZERO
                    MAT(2,4) =       B21

                    MAT(3,1) =        B12
                    MAT(3,2) = ZERO
                    MAT(3,3) = A11 +  B22
                    MAT(3,4) = A21

                    MAT(4,1) = ZERO
                    MAT(4,2) =        B12
                    MAT(4,3) = A12
                    MAT(4,4) = A22 + B22

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(KH,L)
                    RHS(3) = XL(K, LH)
                    RHS(4) = XL(KH,LH)

                END IF
                IF ( IT .EQ. 0 ) THEN
                    CALL DLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL DGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL DGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
                    IF ( SCAL .NE. ONE ) THEN
                        XL(1:M, 1:N) = SCAL*XL(1:M, 1:N)
                        SCALE = SCALE * SCAL
                        INFO = 1
                    END IF
                    EXIT
                END IF

            END DO

            IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                XL(K,L) = RHS(1)
            ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                XL(K,L) = RHS(1)
                XL(K,LH) = RHS(2)
            ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                XL(K,L) = RHS(1)
                XL(KH,L) = RHS(2)
            ELSE
                XL(K,L) = RHS(1)
                XL(KH,L) = RHS(2)
                XL(K, LH) = RHS(3)
                XL(KH,LH) = RHS(4)
            END IF

            ! in current row

            IF ( KH .LT. M ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    XL(KH+1:M, L) = XL(KH+1:M, L) - AL(KH+1:M,K) * XL(K,L)
                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    XL(KH+1:M, L) = XL(KH+1:M, L) - AL(KH+1:M,K) * XL(K,L)
                    XL(KH+1:M, LH) = XL(KH+1:M, LH) - AL(KH+1:M,K) * XL(K,LH)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    XL(KH+1:M,L) = XL(KH+1:M,L) - AL(KH+1:M,K)*XL(K,L) - AL(KH+1:M,KH)*XL(KH,L)
                ELSE
                    XL(KH+1:M,L) = XL(KH+1:M,L) -  XL(K,L) *AL(KH+1:M,K) - XL(KH,L) * AL(KH+1:M,KH)
                    XL(KH+1:M,LH) = XL(KH+1:M,LH) - XL(K,LH)*AL(KH+1:M,K) - XL(KH,LH) * AL(KH+1:M,KH)

                END IF
            END IF

            K = KH + 1
        END DO

        IF ( LH .LT. N ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = LH+1,N
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, L)*BL(L,IT)
                END DO
            ELSE
                DO IT = LH+1, N
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, L)*BL(L,IT)
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, LH)*BL(LH,IT)
                END DO
            ENDIF
        END IF

        L = LH + 1
    END DO

    X(1:M, 1:N) = XL(1:M, 1:N)

ELSE IF ( BTRANSA .AND. BTRANSB) THEN
    ! (T, T)

    AL(1:M,1:M) = TRANSPOSE(A(1:M, 1:M))
    BL(1:N,1:N) = B(1:N, 1:N)
    XL(1:M, 1:N) = X(1:M, 1:N)

    L = N
    DO WHILE ( L .GT. 0 )
        LB = 1
        LH = L
        IF ( L .GT. 1 ) THEN
            IF (BL(L,L-1) .NE. ZERO) THEN
                LB = 2
                L = L - 1
            END IF
        END IF
        LE = L - 1


        K = 1
        DO WHILE ( K .LE. M )
            KB = 1
            IF ( K .LT. M ) THEN
                IF ( ( AL(K, K+1) .NE. ZERO )) THEN
                    KB = 2
                END IF
            END IF
            KH = K + KB - 1

            IT = 0
            DO

                ! Solve Inner System
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 1

                    MAT(1,1) = AL(K,K) + SGN * BL(L,L)
                    RHS(1)   = XL(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 2

                    A11 = AL(K,K)

                    B11 = SGN * BL(L,L)
                    B21 = SGN * BL(LH,L)
                    B12 = SGN * BL(L,LH)
                    B22 = SGN * BL(LH,LH)


                    MAT(1,1) = A11 + B11
                    MAT(1,2) =       B12
                    MAT(2,1) =       B21
                    MAT(2,2) = A11 + B22

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(K,LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2

                    A11 = AL(K,K)
                    A12 = AL(KH,K)
                    A21 = AL(K,KH)
                    A22 = AL(KH,KH)

                    B11 = SGN * BL(L,L)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A21
                    MAT(2,1) = A12
                    MAT(2,2) = A22 + B11

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(KH,L)

                ELSE
                    DIMMAT = 4

                    A11 = AL(K,K)
                    A12 = AL(KH,K)
                    A21 = AL(K,KH)
                    A22 = AL(KH,KH)

                    B11 = SGN * BL(L,L)
                    B12 = SGN * BL(L,LH)
                    B21 = SGN * BL(LH, L)
                    B22 = SGN * BL(LH,LH)

                    MAT(1,1) = A11 + B11
                    MAT(1,2) = A21
                    MAT(1,3) =         + B12
                    MAT(1,4) = ZERO

                    MAT(2,1) = A12
                    MAT(2,2) = A22 + B11
                    MAT(2,3) =  ZERO
                    MAT(2,4) =  B12

                    MAT(3,1) =       B21
                    MAT(3,2) = ZERO
                    MAT(3,3) = A11 + B22
                    MAT(3,4) = A21

                    MAT(4,1) = ZERO
                    MAT(4,2) =       B21
                    MAT(4,3) = A12
                    MAT(4,4) = A22 + B22

                    RHS(1) = XL(K,L)
                    RHS(2) = XL(KH,L)
                    RHS(3) = XL(K, LH)
                    RHS(4) = XL(KH,LH)

                END IF
                IF ( IT .EQ. 0 ) THEN
                    CALL DLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL DGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL DGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
                    IF ( SCAL .NE. ONE ) THEN
                        XL(1:M, 1:N) = SCAL*XL(1:M, 1:N)
                        SCALE = SCALE * SCAL
                        INFO = 1
                    END IF
                    EXIT
                END IF

            END DO


            IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                XL(K,L) = RHS(1)
            ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                XL(K,L) = RHS(1)
                XL(K,LH) = RHS(2)
            ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                XL(K,L) = RHS(1)
                XL(KH,L) = RHS(2)
            ELSE
                XL(K,L) = RHS(1)
                XL(KH,L) = RHS(2)
                XL(K, LH) = RHS(3)
                XL(KH,LH) = RHS(4)
            END IF

            ! in current row
            IF ( KH .LT. M  ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    XL(KH+1:M, L) = XL(KH+1:M, L) - AL(KH+1:M,K) * XL(K,L)

                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    XL(KH+1:M, L) = XL(KH+1:M, L) - AL(KH+1:M,K) * XL(K,L)
                    XL(KH+1:M, LH) = XL(KH+1:M, LH) - AL(KH+1:M,K) * XL(K,LH)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    XL(KH+1:M,L) = XL(KH+1:M,L) - AL(KH+1:M,K)*XL(K,L) - AL(KH+1:M,KH)*XL(KH,L)

                ELSE

                    XL(KH+1:M,L) = XL(KH+1:M,L) -   XL(K,L) *AL(KH+1:M,K) - XL(KH,L) * AL(KH+1:M, KH)
                    XL(KH+1:M,LH) = XL(KH+1:M,LH) - XL(K,LH)*AL(KH+1:M,K) - XL(KH,LH) * AL(KH+1:M, KH)
                END IF
            END IF

            K = KH + 1
        END DO

        IF ( L .GT. IONE ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = 1,LE
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, L)*BL(IT,L)
                END DO
            ELSE
                DO IT = 1,LE
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, L)*BL(IT,L)
                    XL(1:M,IT) =  XL(1:M,IT) - SGN*XL(1:M, LH)*BL(IT, LH)
                END DO
            ENDIF

        END IF
        L = L - 1
    END DO

    X(1:M, 1:N) = XL(1:M, 1:N)

END IF


