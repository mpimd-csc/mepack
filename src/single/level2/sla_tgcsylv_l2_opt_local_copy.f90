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
!> \brief Level-2 Kagstroem-Westin/Bartels-Stewart Algorithm for the generalized coupled Sylvester equation (Optimized)
!
!  Definition:
!  ===========
!       SUBROUTINE SLA_TGCSYLV_L2_LOCAL_COPY ( TRANSA, TRANSB, SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
!                                       & E, LDE, F, LDF,  SCALE, WORK, INFO)
!           CHARACTER TRANSA(1), TRANSB(1)
!           REAL SGN1, SGN2, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
!           REAL A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF, *)
!           REAL WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TGCSYLV_L2_LOCAL_COPY solves a generalized coupled  Sylvester equation of the following form
!>
!>    op1(A) * R  + SGN1 * L  * op2(B) = SCALE * E                              (1)
!>    op1(C) * R  + SGN2 * L  * op2(D) = SCALE * F
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by SGGES from LAPACK.
!> \endverbatim
!> \remark The algorithm is implemented using BLAS level 2 operations and hand written optimizations.
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
!> \param[in] SGN1
!> \verbatim
!>          SGN1 is REAL, allowed values: +/-1
!>          Specifies the sign between in the first equation.
!> \endverbatim
!>
!> \param[in] SGN2
!> \verbatim
!>          SGN2 is REAL, allowed values: +/-1
!>          Specifies the sign between in the second equation.
!> \endverbatim
!>
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
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (LDC,M)
!>          The matrix C must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,M).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (LDD,N)
!>          The matrix D must be (quasi-) upper triangular. If the matrix B is already
!>          quasi-upper triangular the matrix D must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array D.  LDD >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (LDE,N)
!>          On input, the matrix E contains the right hand side E. !>          On output, the matrix E contains the solution R of Equation (1)
!>          On output, the matrix E contains the solution R of Equation (1)
!>          as selected by TRANSA, TRANSB, and SGN1/SGN2.
!>          Right hand side E and the solution R are M-by-N matrices.
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of the array E.  LDE >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] F
!> \verbatim
!>          F is REAL array, dimension (LDF,N)
!>          On input, the matrix F contains the right hand side F.
!>          On output, the matrix F contains the solution L of Equation (1)
!>          as selected by TRANSA, TRANSB, and SGN1/SGN2.
!>          Right hand side F and the solution L are M-by-N matrices.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of the array F.  LDE >= max(1,M).
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
!
!> \par Optimizations:
!       ==============
!>
!> \li Replaced level-3 by level-2 and level-1 calls,
!> \li Reorder the solution order to column first style,
!> \li Replaced SAXPY operation by Fortran intrinsic
!> \li Replaced all BLAS calls except of STRMV and SGER by Fortran intrinsic,
!> \li Use local copies of A, B, C, D, E, and F (M, N <=128) .

!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup sgltgsylv
!
SUBROUTINE SLA_TGCSYLV_L2_LOCAL_COPY ( TRANSA, TRANSB, SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
        & E, LDE, F, LDF,  SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
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
    INTEGER IONE, MAXNB
    PARAMETER(IONE = 1, MAXNB=128)

    REAL AL(MAXNB,MAXNB), BL(MAXNB,MAXNB), CL(MAXNB,MAXNB), DL(MAXNB,MAXNB), EL(MAXNB,MAXNB), FL(MAXNB, MAXNB)
    INTEGER LDAL, LDBL, LDCL, LDDL, LDEL, LDFL

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
        CALL XERROR_HANDLER('SLA_TGCSYLV_L2_LOCAL_COPY', -INFO)
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
        DO L = 1, M
            DO K = 1, M
                AL(K,L) = A(K,L)
                CL(K,L) = C(K,L)
            END DO
        END DO
        DO L = 1, N
            DO K = 1, N
                BL(K,L) = B(K,L)
                DL(K,L) = D(K,L)
            END DO
        END DO
        DO L = 1, N
            DO K = 1, M
                EL(K,L) = E(K,L)
                FL(K,L) = F(K,L)
            END DO
        END DO

        L = 1
        DO WHILE ( L .LE. N )
            LB = 1
            IF ( L .LT. N ) THEN
                IF ( B(L+1,L) .NE. ZERO ) THEN
                    LB = 2
                ELSE IF (D(L+1,L) .NE. ZERO )  THEN
                    LB = 2
                END IF
            END IF
            LH = L + LB - 1

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
                    ! Solve inner system
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
                    CALL SGER(M, N-LH, -SGN1*ONE, FL(1,L), IONE, BL(L,LH+1), LDBL, EL(1,LH+1), LDEL)
                    CALL SGER(M, N-LH, -SGN2*ONE, FL(1,L), IONE, DL(L,LH+1), LDDL, FL(1,LH+1), LDFL)

                ELSE
                    CALL SGER(M, N-LH, -SGN1*ONE, FL(1,L), IONE, BL(L,LH+1), LDBL, EL(1,LH+1), LDEL)
                    CALL SGER(M, N-LH, -SGN1*ONE, FL(1,LH), IONE, BL(LH,LH+1), LDBL, EL(1,LH+1), LDEL)

                    CALL SGER(M, N-LH, -SGN2*ONE, FL(1,L), IONE, DL(L,LH+1), LDDL, FL(1,LH+1), LDFL)
                    CALL SGER(M, N-LH, -SGN2*ONE, FL(1,LH), IONE, DL(LH,LH+1), LDDL, FL(1,LH+1), LDFL)
                ENDIF
            END IF

            L = LH + 1
        END DO

        DO L = 1, N
            DO K = 1, M
                E(K,L) = EL(K,L)
                F(K,L) = FL(K,L)
            END DO
        END DO

    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N,T)
        DO L = 1, M
            DO K = 1, M
                AL(K,L) = A(K,L)
                CL(K,L) = C(K,L)
            END DO
        END DO
        DO L = 1, N
            DO K = 1, N
                BL(K,L) = B(K,L)
                DL(K,L) = D(K,L)
            END DO
        END DO
        DO L = 1, N
            DO K = 1, M
                EL(K,L) = E(K,L)
                FL(K,L) = F(K,L)

            END DO
        END DO

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
                    IF (A(K,K-1) .NE. ZERO ) THEN
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
                    CALL SGER(M, LE, -SGN1*ONE, FL(1,L), IONE, BL(1,L), IONE, EL(1,1), LDEL)
                    CALL SGER(M, LE, -SGN2*ONE, FL(1,L), IONE, DL(1,L), IONE, FL(1,1), LDFL)
                ELSE
                    CALL SGER(M, LE, -SGN1*ONE, FL(1,L), IONE, BL(1,L), IONE, EL(1,1), LDEL)
                    CALL SGER(M, LE, -SGN1*ONE, FL(1,LH), IONE, BL(1,LH), IONE, EL(1,1), LDEL)
                    CALL SGER(M, LE, -SGN2*ONE, FL(1,L), IONE, DL(1,L), IONE, FL(1,1), LDFL)
                    CALL SGER(M, LE, -SGN2*ONE, FL(1,LH), IONE, DL(1,LH), IONE, FL(1,1), LDFL)
                ENDIF
            END IF
            L = L - 1
        END DO

        DO L = 1, N
            DO K = 1, M
                E(K,L) = EL(K,L)
                F(K,L) = FL(K,L)
            END DO
        END DO

    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
        ! (T, N)
        DO L = 1, M
            DO K = 1, M
                AL(K,L) = A(L,K)
                CL(K,L) = C(L,K)
            END DO
        END DO
        DO L = 1, N
            DO K = 1, N
                BL(K,L) = B(K,L)
                DL(K,L) = D(K,L)
            END DO
        END DO
        DO L = 1, N
            DO K = 1, M
                EL(K,L) = E(K,L)
                FL(K,L) = F(K,L)
            END DO
        END DO


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

            IF ( LH .LT. N ) THEN
                IF ( LB .EQ. IONE ) THEN
                    CALL SGER(M, N-LH, -SGN1*ONE, FL(1,L), IONE, BL(L,LH+1), LDBL, EL(1,LH+1), LDEL)
                    CALL SGER(M, N-LH, -SGN2*ONE, FL(1,L), IONE, DL(L,LH+1), LDDL, FL(1,LH+1), LDFL)

                ELSE
                    CALL SGER(M, N-LH, -SGN1*ONE, FL(1,L), IONE, BL(L,LH+1), LDBL, EL(1,LH+1), LDEL)
                    CALL SGER(M, N-LH, -SGN1*ONE, FL(1,LH), IONE, BL(LH,LH+1), LDBL, EL(1,LH+1), LDEL)

                    CALL SGER(M, N-LH, -SGN2*ONE, FL(1,L), IONE, DL(L,LH+1), LDDL, FL(1,LH+1), LDFL)
                    CALL SGER(M, N-LH, -SGN2*ONE, FL(1,LH), IONE, DL(LH,LH+1), LDDL, FL(1,LH+1), LDFL)
                ENDIF
            END IF

            L = LH + 1
        END DO

        DO L = 1, N
            DO K = 1, M
                E(K,L) = EL(K,L)
                F(K,L) = FL(K,L)
            END DO
        END DO

    ELSE IF ( BTRANSA .AND. BTRANSB) THEN
        ! (T, T)

        DO L = 1, M
            DO K = 1, M
                AL(K,L) = A(L,K)
                CL(K,L) = C(L,K)
            END DO
        END DO
        DO L = 1, N
            DO K = 1, N
                BL(K,L) = B(K,L)
                DL(K,L) = D(K,L)
            END DO
        END DO
        DO L = 1, N
            DO K = 1, M
                EL(K,L) = E(K,L)
                FL(K,L) = F(K,L)
            END DO
        END DO


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
                    IF (( A(K+1,K) .NE. ZERO )) THEN
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
                    CALL SGER(M, LE, -SGN1*ONE, FL(1,L), IONE, BL(1,L), IONE, EL(1,1), LDEL)
                    CALL SGER(M, LE, -SGN2*ONE, FL(1,L), IONE, DL(1,L), IONE, FL(1,1), LDFL)
                ELSE
                    CALL SGER(M, LE, -SGN1*ONE, FL(1,L), IONE, BL(1,L), IONE, EL(1,1), LDEL)
                    CALL SGER(M, LE, -SGN1*ONE, FL(1,LH), IONE, BL(1,LH), IONE, EL(1,1), LDEL)
                    CALL SGER(M, LE, -SGN2*ONE, FL(1,L), IONE, DL(1,L), IONE, FL(1,1), LDFL)
                    CALL SGER(M, LE, -SGN2*ONE, FL(1,LH), IONE, DL(1,LH), IONE, FL(1,1), LDFL)
                ENDIF
            END IF

            L = L - 1
        END DO

        DO L = 1, N
            DO K = 1, M
                E(K,L) = EL(K,L)
                F(K,L) = FL(K,L)
            END DO
        END DO

    END IF
END SUBROUTINE



