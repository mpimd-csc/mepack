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


!> \brief Level-2 Bartels-Stewart Algorithm for the discrete time  Sylvester equation (Optimized)
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TRSYLV2_L2_LOCAL_COPY ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANSA(1), TRANSB(1)
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
!> SLA_TRSYLV2_L2_LOCAL_COPY solves a generalized Sylvester equation of the following forms
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
!>
!>
!> \par Optimizations:
!       ==============
!>
!> \li Replaced level-3 by level-2 and level-1 calls,
!> \li Reorder the solution order to column first style,
!> \li Replaced SAXPY operation by Fortran intrinsic
!> \li Replaced all BLAS calls except of STRMV and SGER by Fortran intrinsic,
!> \li Use local copies of A, B, C, D, and X (M, N <=128) .
!
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Dezember 2022
!> \ingroup sgltrsylv
!
SUBROUTINE SLA_TRSYLV2_L2_LOCAL_COPY ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    REAL SGN, SCALE
    INTEGER M, N, LDA, LDB,LDX, INFO
    REAL A(LDA, *), B(LDB, *) , X(LDX, *)
    REAL WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE, IT
    INTEGER DIMMAT, INFO1, IPIV(4), JPIV(4)
    REAL MAT(4,4), RHS(4)
    REAL A11, A12, A21, A22
    REAL B11, B12, B21, B22
    REAL SCAL

    LOGICAL BTRANSA, BTRANSB
    REAL ZERO , ONE
    PARAMETER (ZERO = 0.0, ONE = 1.0)

    INTEGER IONE,MAXNB
    PARAMETER ( IONE = 1,MAXNB=128)

    REAL AL(MAXNB,MAXNB), BL(MAXNB,MAXNB), XL(MAXNB,MAXNB)
    INTEGER LDAL, LDBL, LDXL

    ! EXTERNAL FUNCTIONS
    EXTERNAL SGER
    EXTERNAL SGESC2
    EXTERNAL SGETC2X
    EXTERNAL SLA_SMALL_SOLVE4
    EXTERNAL STRMV
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME


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
        CALL XERROR_HANDLER('SLA_TRSYLV2_L2_LOCAL_COPY', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1) THEN
        INFO = 2*MAX(M,N)
        RETURN
    END IF

    ! Quick Return
    SCALE = ONE
    IF ( M .EQ. 0 .OR. N .EQ. 0 ) THEN
        RETURN
    END IF

    ! Prepare
    INFO = 0

    IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB ) THEN
        ! ( N, N )
        AL(1:M,1:M) = A(1:M,1:M)
        BL(1:N,1:N) = B(1:N,1:N)
        XL(1:M,1:N) = X(1:M,1:N)


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

                        MAT(1,1) = BL(L,L) * AL(K,K) + SGN
                        RHS(1)   = XL(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)

                        B11 = BL(L,L)
                        B21 = BL(LH,L)
                        B12 = BL(L,LH)
                        B22 = BL(LH,LH)


                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B21*A11
                        MAT(2,1) = B12*A11
                        MAT(2,2) = B22*A11 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(K,LH)

                    ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)
                        A12 = AL(K,KH)
                        A21 = AL(KH,K)
                        A22 = AL(KH,KH)
                        B11 = BL(L,L)

                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B11*A12
                        MAT(2,1) = B11*A21
                        MAT(2,2) = B11*A22 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)

                    ELSE
                        DIMMAT = 4

                        A11 = AL(K,K)
                        A12 = AL(K,KH)
                        A21 = AL(KH,K)
                        A22 = AL(KH,KH)

                        B11 = BL(L,L)
                        B12 = BL(L,LH)
                        B21 = BL(LH, L)
                        B22 = BL(LH,LH)

                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B11*A12
                        MAT(1,3) = B21*A11
                        MAT(1,4) = B21*A12

                        MAT(2,1) = B11*A21
                        MAT(2,2) = B11*A22 + SGN
                        MAT(2,3) = B21*A21
                        MAT(2,4) = B21*A22

                        MAT(3,1) = B12*A11
                        MAT(3,2) = B12*A12
                        MAT(3,3) = B22*A11 + SGN
                        MAT(3,4) = B22*A12

                        MAT(4,1) = B12*A21
                        MAT(4,2) = B12*A22
                        MAT(4,3) = B22*A21
                        MAT(4,4) = B22*A22 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)
                        RHS(3) = XL(K, LH)
                        RHS(4) = XL(KH,LH)

                    END IF

                    IF ( IT .EQ. 0 ) THEN
                        CALL SLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                        IF ( INFO1 .EQ. 0 ) EXIT
                        IT = 1
                    ELSE
                        CALL SGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                        IF ( INFO1 .NE. 0 ) INFO = INFO1
                        CALL SGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
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
                        XL(1:K-1,L) = XL(1:K-1,L) - XL(K,L) * BL(L,L) *AL(1:K-1,K)
                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                        XL(1:K-1,L) = XL(1:K-1,L) - (XL(K,L) * BL(L,L)  + XL(K,LH) * BL(LH,L))*AL(1:K-1,K)
                        XL(1:K-1,LH) = XL(1:K-1,LH) - (XL(K,L) * BL(L,LH)  + XL(K,LH) * BL(LH,LH))*AL(1:K-1,K)
                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                        XL(1:K-1, L) = XL(1:K-1,L) - AL(1:K-1,K) *(XL(K,L) * BL(L,L)) - AL(1:K-1,KH) *(XL(KH,L) * BL(L,L))
                    ELSE
                        MAT(1,1) = XL(K,L) *BL(L,L)  + XL(K,LH)  * BL(LH,L)
                        MAT(2,1) = XL(KH,L)*BL(L,L)  + XL(KH,LH) * BL(LH,L)
                        MAT(3,1) = XL(K,L) *BL(L,LH) + XL(K,LH)  * BL(LH,LH)
                        MAT(4,1) = XL(KH,L)*BL(L,LH) + XL(KH,LH) * BL(LH,LH)

                        XL(1:K-1,L) = XL(1:K-1,L) - MAT(1,1) * AL(1:K-1,K) - MAT(2,1) * AL(1:K-1,KH)
                        XL(1:K-1,LH) = XL(1:K-1,LH) - MAT(3,1) * AL(1:K-1,K) - MAT(4,1) * AL(1:K-1,KH)

                    END IF

                END IF

                K = K - 1

            END DO

            ! Update January 2021
            IF ( LH .LT. N ) THEN
                IF ( LB .EQ. IONE ) THEN
                    WORK(1:M) = XL(1:M,L)
                    ! CALL SLACPY('All', M, LB, XL(1,L), LDXL,  WORK, M)
                    CALL STRMV('Upper', 'NoTrans','NoUnit', M, AL, LDAL, WORK, IONE)
                    DO IT = 1,M-1
                        IF (AL(IT+1,IT).NE.ZERO) THEN
                            WORK(IT+1) = WORK(IT+1) + AL(IT+1,IT) * XL(IT,L)
                        END IF
                    END DO

                    CALL SGER(M, N-LH, -ONE, WORK, IONE, BL(L,LH+1), LDBL, XL(1,LH+1), LDXL)

                ELSE
                    WORK(1:M) = XL(1:M,L)
                    WORK(M+1:2*M) = XL(1:M,LH)
                    CALL STRMV('Upper', 'NoTrans','NoUnit', M, AL, LDAL, WORK, IONE)
                    CALL STRMV('Upper', 'NoTrans','NoUnit', M, AL, LDAL, WORK(M+1), IONE)

                    DO IT = 1,M-1
                        IF (AL(IT+1,IT).NE.ZERO) THEN
                            WORK(IT+1) = WORK(IT+1) + AL(IT+1,IT) * XL(IT,L)
                            WORK(M+IT+1) = WORK(M+IT+1) + AL(IT+1,IT) * XL(IT,LH)

                        END IF
                    END DO

                    CALL SGER(M, N-LH, -ONE, WORK, IONE, BL(L,LH+1), LDBL, XL(1,LH+1), LDXL)
                    CALL SGER(M, N-LH, -ONE, WORK(M+1), IONE, BL(LH,LH+1), LDBL, XL(1,LH+1), LDXL)

                ENDIF


            END IF

            L = LH + 1
        END DO

        X(1:M,1:N) = XL(1:M,1:N)

    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N,T)
        AL(1:M,1:M) = A(1:M,1:M)
        BL(1:N,1:N) = B(1:N,1:N)
        XL(1:M,1:N) = X(1:M,1:N)

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

                        MAT(1,1) = BL(L,L) * AL(K,K) + SGN
                        RHS(1)   = XL(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)
                        B11 = BL(L,L)
                        B21 = BL(LH,L)
                        B12 = BL(L,LH)
                        B22 = BL(LH,LH)


                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B12*A11
                        MAT(2,1) = B21*A11
                        MAT(2,2) = B22*A11 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(K,LH)

                    ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)
                        A12 = AL(K,KH)
                        A21 = AL(KH,K)
                        A22 = AL(KH,KH)
                        B11 = BL(L,L)

                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B11*A12
                        MAT(2,1) = B11*A21
                        MAT(2,2) = B11*A22 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)

                    ELSE
                        DIMMAT = 4

                        A11 = AL(K,K)
                        A12 = AL(K,KH)
                        A21 = AL(KH,K)
                        A22 = AL(KH,KH)

                        B11 = BL(L,L)
                        B12 = BL(L,LH)
                        B21 = BL(LH, L)
                        B22 = BL(LH,LH)

                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B11*A12
                        MAT(1,3) = B12*A11
                        MAT(1,4) = B12*A12

                        MAT(2,1) = B11*A21
                        MAT(2,2) = B11*A22 + SGN
                        MAT(2,3) = B12*A21
                        MAT(2,4) = B12*A22

                        MAT(3,1) = B21*A11
                        MAT(3,2) = B21*A12
                        MAT(3,3) = B22*A11 + SGN
                        MAT(3,4) = B22*A12

                        MAT(4,1) = B21*A21
                        MAT(4,2) = B21*A22
                        MAT(4,3) = B22*A21
                        MAT(4,4) = B22*A22 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)
                        RHS(3) = XL(K, LH)
                        RHS(4) = XL(KH,LH)

                    END IF

                    IF ( IT .EQ. 0 ) THEN
                        CALL SLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                        IF ( INFO1 .EQ. 0 ) EXIT
                        IT = 1
                    ELSE
                        CALL SGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                        IF ( INFO1 .NE. 0 ) INFO = INFO1
                        CALL SGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
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
                        XL(1:KE, L ) = XL(1:KE, L) - (BL(L,L) * XL(K,L))* AL(1:KE, K)
                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                        XL(1:KE,L) = XL(1:KE,L) - (XL(K,L) * BL(L,L)  + XL(K,LH) * BL(L,LH))*AL(1:KE,K)
                        XL(1:KE,LH) = XL(1:KE,LH) - (XL(K,L) * BL(LH,L)  + XL(K,LH) * BL(LH,LH))*AL(1:KE,K)
                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                        XL(1:KE, L) = XL(1:KE,L) - AL(1:KE,K) *(XL(K,L) * BL(L,L)) - AL(1:KE,KH) *(XL(KH,L) * BL(L,L))

                    ELSE
                        MAT(1,1) = XL(K,L) *BL(L,L)  + XL(K,LH)  * BL(L,LH)
                        MAT(2,1) = XL(KH,L)*BL(L,L)  + XL(KH,LH) * BL(L,LH)
                        MAT(3,1) = XL(K,L) *BL(LH,L) + XL(K,LH)  * BL(LH,LH)
                        MAT(4,1) = XL(KH,L)*BL(LH,L) + XL(KH,LH) * BL(LH,LH)


                        XL(1:KE,L) = XL(1:KE,L) - MAT(1,1) * AL(1:KE,K) - MAT(2,1) * AL(1:KE,KH)
                        XL(1:KE,LH) = XL(1:KE,LH) - MAT(3,1) * AL(1:KE,K) - MAT(4,1) * AL(1:KE,KH)
                    END IF

                END IF

                K = K - 1
            END DO

            ! Update January 2021
            IF ( L .GT. 1 ) THEN
                IF ( LB .EQ. IONE ) THEN
                    WORK(1:M) = XL(1:M,L)
                    CALL STRMV('Upper', 'NoTrans','NoUnit', M, AL, LDAL, WORK, IONE)
                    DO IT = 1,M-1
                        IF (AL(IT+1,IT).NE.ZERO) THEN
                            WORK(IT+1) = WORK(IT+1) + AL(IT+1,IT) * XL(IT,L)
                        END IF
                    END DO
                    CALL SGER(M, LE, -ONE, WORK, IONE, BL(1,L), IONE, XL(1,1), LDXL)

                ELSE
                    WORK(1:M) = XL(1:M,L)
                    WORK(M+1:2*M) = XL(1:M,LH)
                    CALL STRMV('Upper', 'NoTrans','NoUnit', M, AL, LDAL, WORK, IONE)
                    CALL STRMV('Upper', 'NoTrans','NoUnit', M, AL, LDAL, WORK(M+1), IONE)

                    DO IT = 1,M-1
                        IF (AL(IT+1,IT).NE.ZERO) THEN
                            WORK(IT+1) = WORK(IT+1) + AL(IT+1,IT) * XL(IT,L)
                            WORK(M+IT+1) = WORK(M+IT+1) + AL(IT+1,IT) * XL(IT,LH)

                        END IF
                    END DO
                    CALL SGER(M, LE, -ONE, WORK, IONE, BL(1,L), IONE, XL(1,1), LDXL)
                    CALL SGER(M, LE, -ONE, WORK(M+1), IONE, BL(1,LH), IONE, XL(1,1), LDXL)
                ENDIF
            END IF
            L = L - 1
        END DO

        X(1:M,1:N) = XL(1:M,1:N)


    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
        ! (T, N)
        AL(1:M,1:M) = TRANSPOSE(A(1:M,1:M))
        BL(1:N,1:N) = B(1:N,1:N)
        XL(1:M,1:N) = X(1:M,1:N)

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

                        MAT(1,1) = BL(L,L) * AL(K,K) + SGN
                        RHS(1)   = XL(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)

                        B11 = BL(L,L)
                        B21 = BL(LH,L)
                        B12 = BL(L,LH)
                        B22 = BL(LH,LH)


                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B21*A11
                        MAT(2,1) = B12*A11
                        MAT(2,2) = B22*A11 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(K,LH)

                    ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)
                        A12 = AL(KH,K)
                        A21 = AL(K,KH)
                        A22 = AL(KH,KH)

                        B11 = BL(L,L)

                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B11*A21
                        MAT(2,1) = B11*A12
                        MAT(2,2) = B11*A22 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)

                    ELSE
                        DIMMAT = 4

                        A11 = AL(K,K)
                        A12 = AL(KH,K)
                        A21 = AL(K,KH)
                        A22 = AL(KH,KH)
                        B11 = BL(L,L)
                        B12 = BL(L,LH)
                        B21 = BL(LH, L)
                        B22 = BL(LH,LH)

                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B11*A21
                        MAT(1,3) = B21*A11
                        MAT(1,4) = B21*A21

                        MAT(2,1) = B11*A12
                        MAT(2,2) = B11*A22 + SGN
                        MAT(2,3) = B21*A12
                        MAT(2,4) = B21*A22

                        MAT(3,1) = B12*A11
                        MAT(3,2) = B12*A21
                        MAT(3,3) = B22*A11 + SGN
                        MAT(3,4) = B22*A21

                        MAT(4,1) = B12*A12
                        MAT(4,2) = B12*A22
                        MAT(4,3) = B22*A12
                        MAT(4,4) = B22*A22 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)
                        RHS(3) = XL(K, LH)
                        RHS(4) = XL(KH,LH)

                    END IF
                    IF ( IT .EQ. 0 ) THEN
                        CALL SLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                        IF ( INFO1 .EQ. 0 ) EXIT
                        IT = 1
                    ELSE
                        CALL SGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                        IF ( INFO1 .NE. 0 ) INFO = INFO1
                        CALL SGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
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
                        XL(KH+1:M, L) = XL(KH+1:M, L) - AL(KH+1:M,K) * (XL(K,L) * BL(L,L))

                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                        XL(KH+1:M, L) = XL(KH+1:M, L) - AL(KH+1:M,K) * (XL(K,L) * BL(L,L) + XL(K,LH) * BL(LH,L) )
                        XL(KH+1:M, LH) = XL(KH+1:M, LH) - AL(KH+1:M,K) * (XL(K,L) * BL(L,LH) + XL(K,LH) * BL(LH,LH) )

                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                        XL(KH+1:M,L) = XL(KH+1:M,L) - AL(KH+1:M,K)*(XL(K,L)*BL(L,L)) - AL(KH+1:M,KH)*(XL(KH,L)*BL(L,L))

                    ELSE
                        MAT(1,1) = XL(K,L)*BL(L,L) + XL(K,LH) * BL(LH,L)
                        MAT(2,1) = XL(KH,L)*BL(L,L) + XL(KH,LH) * BL(LH,L)
                        MAT(3,1) = XL(K,L)*BL(L,LH) + XL(K,LH) * BL(LH,LH)
                        MAT(4,1) = XL(KH,L)*BL(L,LH) + XL(KH,LH) * BL(LH,LH)

                        XL(KH+1:M,L) = XL(KH+1:M,L) - MAT(1,1)*AL(KH+1:M,K) - MAT(2,1) * AL(KH+1:M,KH)
                        XL(KH+1:M,LH) = XL(KH+1:M,LH) - MAT(3,1)*AL(KH+1:M,K) - MAT(4,1) * AL(KH+1:M,KH)
                    END IF
                END IF

                K = KH + 1
            END DO

            IF ( LH .LT. N ) THEN
                IF ( LB .EQ. IONE ) THEN
                    WORK(1:M) = XL(1:M, L)
                    CALL STRMV("Lower", "NoTrans", "NoUnit", M, AL, LDAL, WORK, IONE)
                    DO IT = 1, M-1
                        IF (AL(IT,IT+1).NE.ZERO) THEN
                            WORK(IT) = WORK(IT) + AL(IT,IT+1) * XL(IT+1,L)
                        END IF
                    END DO
                    CALL SGER(M, N-LH, -ONE, WORK, IONE, BL(L,LH+1), LDBL, XL(1,LH+1), LDXL)

                ELSE
                    WORK(1:M) = XL(1:M, L)
                    WORK(M+1:2*M) = XL(1:M, LH)
                    CALL STRMV("Lower", "NoTrans", "NoUnit", M, AL, LDAL, WORK, IONE)
                    CALL STRMV("Lower", "NoTrans", "NoUnit", M, AL, LDAL, WORK(M+1), IONE)

                    DO IT = 1, M-1
                        IF (AL(IT,IT+1).NE.ZERO) THEN
                            WORK(IT) = WORK(IT) +     AL(IT, IT+1) * XL(IT+1,L)
                            WORK(IT+M) = WORK(IT+M) + AL(IT, IT+1) * XL(IT+1,LH)
                        END IF
                    END DO
                    CALL SGER(M, N-LH, -ONE, WORK, IONE, BL(L,LH+1), LDBL, XL(1,LH+1), LDXL)
                    CALL SGER(M, N-LH, -ONE, WORK(M+1), IONE, BL(LH,LH+1), LDBL, XL(1,LH+1), LDXL)


                ENDIF
            END IF

            L = LH + 1
        END DO

        X(1:M,1:N) = XL(1:M,1:N)


    ELSE IF ( BTRANSA .AND. BTRANSB) THEN
        ! (T, T)
        AL(1:M,1:M) = TRANSPOSE(A(1:M,1:M))
        BL(1:N,1:N) = B(1:N,1:N)
        XL(1:M,1:N) = X(1:M,1:N)


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

                        MAT(1,1) = BL(L,L) * AL(K,K) + SGN
                        RHS(1)   = XL(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)

                        B11 = BL(L,L)
                        B21 = BL(LH,L)
                        B12 = BL(L,LH)
                        B22 = BL(LH,LH)


                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B12*A11
                        MAT(2,1) = B21*A11
                        MAT(2,2) = B22*A11 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(K,LH)

                    ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 2

                        A11 = AL(K,K)
                        A12 = AL(KH,K)
                        A21 = AL(K,KH)
                        A22 = AL(KH,KH)
                        B11 = BL(L,L)

                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B11*A21
                        MAT(2,1) = B11*A12
                        MAT(2,2) = B11*A22 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)

                    ELSE
                        DIMMAT = 4

                        A11 = AL(K,K)
                        A12 = AL(KH,K)
                        A21 = AL(K,KH)
                        A22 = AL(KH,KH)

                        B11 = BL(L,L)
                        B12 = BL(L,LH)
                        B21 = BL(LH, L)
                        B22 = BL(LH,LH)

                        MAT(1,1) = B11*A11 + SGN
                        MAT(1,2) = B11*A21
                        MAT(1,3) = B12*A11
                        MAT(1,4) = B12*A21

                        MAT(2,1) = B11*A12
                        MAT(2,2) = B11*A22 + SGN
                        MAT(2,3) = B12*A12
                        MAT(2,4) = B12*A22

                        MAT(3,1) = B21*A11
                        MAT(3,2) = B21*A21
                        MAT(3,3) = B22*A11 + SGN
                        MAT(3,4) = B22*A21

                        MAT(4,1) = B21*A12
                        MAT(4,2) = B21*A22
                        MAT(4,3) = B22*A12
                        MAT(4,4) = B22*A22 + SGN

                        RHS(1) = XL(K,L)
                        RHS(2) = XL(KH,L)
                        RHS(3) = XL(K, LH)
                        RHS(4) = XL(KH,LH)

                    END IF
                    IF ( IT .EQ. 0 ) THEN
                        CALL SLA_SMALL_SOLVE4( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                        IF ( INFO1 .EQ. 0 ) EXIT
                        IT = 1
                    ELSE
                        CALL SGETC2X(DIMMAT, MAT, 4, IPIV, JPIV, INFO1)
                        IF ( INFO1 .NE. 0 ) INFO = INFO1
                        CALL SGESC2(DIMMAT, MAT, 4, RHS, IPIV, JPIV, SCAL)
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
                        XL(KH+1:M, L) = XL(KH+1:M, L) - AL(KH+1:M,K) * (XL(K,L) * BL(L,L))
                    ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                        XL(KH+1:M, L) = XL(KH+1:M, L) - AL(KH+1:M,K) * (XL(K,L) * BL(L,L) + XL(K,LH) * BL(L,LH) )
                        XL(KH+1:M, LH) = XL(KH+1:M, LH) - AL(KH+1:M,K) * (XL(K,L) * BL(LH,L) + XL(K,LH) * BL(LH,LH) )
                    ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                        XL(KH+1:M,L) = XL(KH+1:M,L) - AL(KH+1:M,K)*(XL(K,L)*BL(L,L)) - AL(KH+1:M,KH)*(XL(KH,L)*BL(L,L))
                    ELSE
                        MAT(1,1) = XL(K,L)*BL(L,L) + XL(K,LH) * BL(L,LH)
                        MAT(2,1) = XL(KH,L)*BL(L,L) + XL(KH,LH) * BL(L,LH)
                        MAT(3,1) = XL(K,L)*BL(LH,L) + XL(K,LH) * BL(LH,LH)
                        MAT(4,1) = XL(KH,L)*BL(LH,L) + XL(KH,LH) * BL(LH,LH)


                        XL(KH+1:M,L) = XL(KH+1:M,L) - MAT(1,1)*AL(KH+1:M,K) - MAT(2,1) * AL(KH+1:M,KH)
                        XL(KH+1:M,LH) = XL(KH+1:M,LH) - MAT(3,1)*AL(KH+1:M,K) - MAT(4,1) * AL(KH+1:M,KH)
                    END IF
                END IF

                K = KH + 1
            END DO

            IF ( L .GT. IONE ) THEN
                IF ( LB .EQ. IONE ) THEN
                    WORK(1:M) = XL(1:M, L)
                    CALL STRMV("Lower", "NoTrans", "NoUnit", M, AL, LDAL, WORK, IONE)
                    DO IT = 1, M-1
                        IF (AL(IT,IT+1).NE.ZERO) THEN
                            WORK(IT) = WORK(IT) + AL(IT,IT+1) * XL(IT+1,L)
                        END IF
                    END DO

                    CALL SGER(M, LE, -ONE, WORK, IONE, BL(1,L), IONE, XL(1,1), LDXL)


                ELSE
                    WORK(1:M) = XL(1:M,L)
                    WORK(M+1:2*M) = XL(1:M,LH)
                    CALL STRMV('Lower', 'NoTrans','NoUnit', M, AL, LDAL, WORK, IONE)
                    CALL STRMV('Lower', 'NoTrans','NoUnit', M, AL, LDAL, WORK(M+1), IONE)

                    DO IT = 1,M-1
                        IF (AL(IT, IT+1).NE.ZERO) THEN
                            WORK(IT) = WORK(IT)     + AL(IT, IT+1) * XL(IT+1,L)
                            WORK(M+IT) = WORK(M+IT) + AL(IT, IT+1) * XL(IT+1,LH)
                        END IF
                    END DO
                    CALL SGER(M, LE, -ONE, WORK, IONE, BL(1,L), IONE, XL(1,1), LDXL)
                    CALL SGER(M, LE, -ONE, WORK(M+1), IONE, BL(1,LH), IONE, XL(1,1), LDXL)
                ENDIF
            END IF
            L = L - 1
        END DO

        X(1:M,1:N) = XL(1:M, 1:N)
    END IF

END SUBROUTINE



