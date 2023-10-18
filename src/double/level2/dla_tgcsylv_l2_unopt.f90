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


!> \brief Level-2 Kagstroem-Westin/Bartels-Stewart Algorithm for the generalized coupled Sylvester equation (unoptimized)
!
!  Definition:
!  ===========
!       SUBROUTINE DLA_TGCSYLV_L2_UNOPT ( TRANSA, TRANSB, SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
!                                       & E, LDE, F, LDF,  SCALE, WORK, INFO)
!           CHARACTER TRANSA(1), TRANSB(1)
!           DOUBLE PRECISION SGN1, SGN2, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
!           DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF, *)
!           DOUBLE PRECISION WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLA_TGCSYLV_L2_UNOPT solves a generalized coupled  Sylvester equation of the following form
!>
!>    op1(A) * R  + SGN1 * L  * op2(B) = SCALE * E                              (1)
!>    op1(C) * R  + SGN2 * L  * op2(D) = SCALE * F
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by DGGES from LAPACK.
!> \endverbatim
!> \attention The algorithm is implemented using BLAS level 2 operations without further optimizations. For a faster implementation
!> see DLA_TGSYLV_L2 .
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
!>          SGN1 is DOUBLE PRECISION, allowed values: +/-1
!>          Specifies the sign between in the first equation.
!> \endverbatim
!>
!> \param[in] SGN2
!> \verbatim
!>          SGN2 is DOUBLE PRECISION, allowed values: +/-1
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
!>          A is DOUBLE PRECISION array, dimension (LDA,M)
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
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
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
!>          C is DOUBLE PRECISION array, dimension (LDC,M)
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
!>          D is DOUBLE PRECISION array, dimension (LDD,N)
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
!>          E is DOUBLE PRECISION array, dimension (LDE,N)
!>          On input, the matrix E contains the right hand side E.
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
!>          F is DOUBLE PRECISION array, dimension (LDF,N)
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
!>          SCALE is DOUBLE PRECISION
!>          SCALE is a scaling factor to prevent the overflow in the result.
!>          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
!>          could not be solved correctly, 0 < SCALE <= 1 holds true.
!> \endverbatim
!>
!> \param[in] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension LWORK
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
!> \li Nothing.
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date October 2023
!> \ingroup dbltgsylv
!
SUBROUTINE DLA_TGCSYLV_L2_UNOPT ( TRANSA, TRANSB, SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
        & E, LDE, F, LDF,  SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    DOUBLE PRECISION SGN1, SGN2, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF,*)
    DOUBLE PRECISION WORK(*)


    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE, IT
    INTEGER DIMMAT, INFO1, IPIV(8), JPIV(8)
    DOUBLE PRECISION MAT(8,8), RHS(8)
    DOUBLE PRECISION SCAL

    LOGICAL BTRANSA, BTRANSB
    DOUBLE PRECISION ZERO , ONE
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

    ! EXTERNAL FUNCTIONS
    EXTERNAL DGEMM
    EXTERNAL DGESC2
    EXTERNAL DGETC2X
    EXTERNAL DLA_SMALL_SOLVE8
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL  LSAME

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
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -5
    ELSE IF ( N .LT. 0 ) THEN
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
        CALL XERROR_HANDLER('DLA_TGCSYLV_L2_UNOPT', -INFO)
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
        ! ( N, N )
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

            L = 1
            DO WHILE ( L .LE. N )
                LB = 1
                IF ( L .LT. N ) THEN
                    IF ( ( B(L+1,L) .NE. ZERO .OR. D(L+1,L) .NE. ZERO ) ) THEN
                        LB = 2
                    END IF
                END IF
                LH = L + LB - 1
                IT = 0
                DO
                    ! Solve inner system
                    IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 2
                        MAT(1,1) = A(K,K)
                        MAT(1,2) = SGN1 * B(L,L)
                        MAT(2,1) = C(K,K)
                        MAT(2,2) = SGN2 * D(L,L)
                        RHS(1) = E(K,L)
                        RHS(2) = F(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 4
                        MAT(1:4,1:4) = ZERO
                        MAT(1,1) = A(K,K)
                        MAT(1,3) = SGN1 * B(L,L)
                        MAT(1,4) = SGN1 * B(LH, L)

                        MAT(2,2) = A(K,K)
                        MAT(2,3) = SGN1 * B(L,LH)
                        MAT(2,4) = SGN1 * B(LH, LH)

                        MAT(3,1) = C(K,K)
                        MAT(3,3) = SGN2 * D(L,L)
                        MAT(3,4) = SGN2 * D(LH, L)

                        MAT(4,2) = C(K,K)
                        MAT(4,3) = SGN2 * D(L,LH)
                        MAT(4,4) = SGN2 * D(LH, LH)

                        RHS(1) = E(K, L)
                        RHS(2) = E(K, LH)
                        RHS(3) = F(K, L)
                        RHS(4) = F(K, LH)

                    ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 4
                        MAT(1:4,1:4) = ZERO

                        MAT(1,1) = A(K,K)
                        MAT(1,2) = A(K,KH)
                        MAT(1,3) = SGN1 * B(L,L)

                        MAT(2,1) = A(KH,K)
                        MAT(2,2) = A(KH, KH)
                        MAT(2,4) = SGN1 * B(L,L)

                        MAT(3,1) = C(K,K)
                        MAT(3,2) = C(K,KH)
                        MAT(3,3) = SGN2 * D(L,L)

                        MAT(4,2) = C(KH, KH)
                        MAT(4,4) = SGN2 * D(L,L)

                        RHS(1) = E(K, L)
                        RHS(2) = E(KH, L)
                        RHS(3) = F(K, L)
                        RHS(4) = F(KH, L)
                    ELSE
                        DIMMAT = 8

                        MAT(1:8,1:8) = ZERO
                        MAT(1,1) = A(K,K)
                        MAT(1,2) = A(K,KH)
                        MAT(1,5) = SGN1 * B(L,L)
                        MAT(1,7) = SGN1 * B(LH, L)

                        MAT(2,1) = A(KH,K)
                        MAT(2,2) = A(KH,KH)
                        MAT(2,6) = SGN1 * B(L,L)
                        MAT(2,8) = SGN1 * B(LH, L)

                        MAT(3,3) = A(K,K)
                        MAT(3,4) = A(K,KH)
                        MAT(3,5) = SGN1 * B(L,LH)
                        MAT(3,7) = SGN1 * B(LH,LH)

                        MAT(4,3) = A(KH,K)
                        MAT(4,4) = A(KH,KH)
                        MAT(4,6) = SGN1 * B(L,LH)
                        MAT(4,8) = SGN1 * B(LH,LH)

                        MAT(5,1) = C(K,K)
                        MAT(5,2) = C(K,KH)
                        MAT(5,5) = SGN2 * D(L,L)
                        MAT(5,7) = SGN2 * D(LH, L)

                        MAT(6,2) = C(KH,KH)
                        MAT(6,6) = SGN2 * D(L,L)
                        MAT(6,8) = SGN2 * D(LH, L)

                        MAT(7,3) = C(K,K)
                        MAT(7,4) = C(K,KH)
                        MAT(7,5) = SGN2 * D(L,LH)
                        MAT(7,7) = SGN2 * D(LH,LH)

                        MAT(8,4) = C(KH,KH)
                        MAT(8,6) = SGN2 * D(L,LH)
                        MAT(8,8) = SGN2 * D(LH,LH)


                        RHS(1) = E(K,L)
                        RHS(2) = E(KH,L)
                        RHS(3) = E(K,LH)
                        RHS(4) = E(KH,LH)

                        RHS(5) = F(K,L)
                        RHS(6) = F(KH,L)
                        RHS(7) = F(K,LH)
                        RHS(8) = F(KH,LH)
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
                            E(1:M, 1:N) = SCAL*E(1:M, 1:N)
                            F(1:M, 1:N) = SCAL*F(1:M, 1:N)

                            SCALE = SCALE * SCAL
                            INFO = 1
                        END IF
                        EXIT
                    END IF

                END DO

                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    E(K,L) = RHS(1)
                    F(K,L) = RHS(2)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN

                    E(K, L)  = RHS(1)
                    E(K, LH)  = RHS(2)
                    F(K, L)  = RHS(3)
                    F(K, LH)  = RHS(4)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    E(K, L)  = RHS(1)
                    E(KH, L)  = RHS(2)
                    F(K, L)  = RHS(3)
                    F(KH, L)  = RHS(4)
                ELSE
                    E(K,L) = RHS(1)
                    E(KH,L)  = RHS(2)
                    E(K,LH)  = RHS(3)
                    E(KH,LH)  = RHS(4)

                    F(K,L) = RHS(5)
                    F(KH,L)  = RHS(6)
                    F(K,LH)  = RHS(7)
                    F(KH,LH)  = RHS(8)
                END IF

                ! Update January 2021
                ! in current row
                IF ( LH .LT. N ) THEN
                    CALL DGEMM('N','N', KB, N-LH, LB, -SGN1*ONE, F(K,L), LDF, B(L,LH+1), LDB, ONE, E(K,LH+1), LDE)
                    CALL DGEMM('N','N', KB, N-LH, LB, -SGN2*ONE, F(K,L), LDF, D(L,LH+1), LDD, ONE, F(K,LH+1), LDF)
                END IF

                L = LH + 1
            END DO

            ! Update January 2021
            IF ( K .GT. 1 ) THEN
                CALL DGEMM('N', 'N', KE, N, KB, -ONE, A(1,K), LDA, E(K,1), LDE, ONE, E(1,1), LDE)
                CALL DGEMM('N', 'N', KE, N, KB, -ONE, C(1,K), LDA, E(K,1), LDE, ONE, F(1,1), LDF)
            END IF
            K = K - 1
        END DO

    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N,T)
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

            L = N
            DO WHILE ( L .GT. 0 )
                LB = 1
                LH = L
                IF ( L .GT. 1) THEN
                    IF ( (B(L,L-1) .NE. ZERO .OR. D(L,L-1) .NE. ZERO) ) THEN
                        LB = 2
                        L = L - 1
                    END IF
                END IF
                LE = L - 1

                IT = 0
                DO

                    ! Solve Inner System
                    IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 2
                        MAT(1,1) = A(K,K)
                        MAT(1,2) = SGN1 * B(L,L)
                        MAT(2,1) = C(K,K)
                        MAT(2,2) = SGN2 * D(L,L)
                        RHS(1) = E(K,L)
                        RHS(2) = F(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 4
                        MAT(1:4,1:4) = ZERO
                        MAT(1,1) = A(K,K)
                        MAT(1,3) = SGN1 * B(L,L)
                        MAT(1,4) = SGN1 * B(L, LH)

                        MAT(2,2) = A(K,K)
                        MAT(2,3) = SGN1 * B(LH,L)
                        MAT(2,4) = SGN1 * B(LH, LH)

                        MAT(3,1) = C(K,K)
                        MAT(3,3) = SGN2 * D(L,L)
                        MAT(3,4) = SGN2 * D(L, LH)

                        MAT(4,2) = C(K,K)
                        MAT(4,3) = SGN2 * D(LH, L)
                        MAT(4,4) = SGN2 * D(LH, LH)

                        RHS(1) = E(K, L)
                        RHS(2) = E(K, LH)
                        RHS(3) = F(K, L)
                        RHS(4) = F(K, LH)

                    ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 4
                        MAT(1:4,1:4) = ZERO

                        MAT(1,1) = A(K,K)
                        MAT(1,2) = A(K,KH)
                        MAT(1,3) = SGN1 * B(L,L)

                        MAT(2,1) = A(KH,K)
                        MAT(2,2) = A(KH, KH)
                        MAT(2,4) = SGN1 * B(L,L)

                        MAT(3,1) = C(K,K)
                        MAT(3,2) = C(K,KH)
                        MAT(3,3) = SGN2 * D(L,L)

                        MAT(4,2) = C(KH, KH)
                        MAT(4,4) = SGN2 * D(L,L)

                        RHS(1) = E(K, L)
                        RHS(2) = E(KH, L)
                        RHS(3) = F(K, L)
                        RHS(4) = F(KH, L)
                    ELSE
                        DIMMAT = 8
                        MAT(1:8,1:8) = ZERO

                        MAT(1,1) = A(K,K)
                        MAT(1,2) = A(K,KH)
                        MAT(1,5) = SGN1 * B(L,L)
                        MAT(1,7) = SGN1 * B(L, LH)

                        MAT(2,1) = A(KH,K)
                        MAT(2,2) = A(KH,KH)
                        MAT(2,6) = SGN1 * B(L,L)
                        MAT(2,8) = SGN1 * B(L, LH)

                        MAT(3,3) = A(K,K)
                        MAT(3,4) = A(K,KH)
                        MAT(3,5) = SGN1 * B(LH,L)
                        MAT(3,7) = SGN1 * B(LH,LH)

                        MAT(4,3) = A(KH,K)
                        MAT(4,4) = A(KH,KH)
                        MAT(4,6) = SGN1 * B(LH,L)
                        MAT(4,8) = SGN1 * B(LH,LH)

                        MAT(5,1) = C(K,K)
                        MAT(5,2) = C(K,KH)
                        MAT(5,5) = SGN2 * D(L,L)
                        MAT(5,7) = SGN2 * D(L, LH)

                        MAT(6,2) = C(KH,KH)
                        MAT(6,6) = SGN2 * D(L,L)
                        MAT(6,8) = SGN2 * D(L, LH)

                        MAT(7,3) = C(K,K)
                        MAT(7,4) = C(K,KH)
                        MAT(7,5) = SGN2 * D(LH,L)
                        MAT(7,7) = SGN2 * D(LH,LH)

                        MAT(8,4) = C(KH,KH)
                        MAT(8,6) = SGN2 * D(LH,L)
                        MAT(8,8) = SGN2 * D(LH,LH)


                        RHS(1) = E(K,L)
                        RHS(2) = E(KH,L)
                        RHS(3) = E(K,LH)
                        RHS(4) = E(KH,LH)

                        RHS(5) = F(K,L)
                        RHS(6) = F(KH,L)
                        RHS(7) = F(K,LH)
                        RHS(8) = F(KH,LH)

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
                            E(1:M, 1:N) = SCAL*E(1:M, 1:N)
                            F(1:M, 1:N) = SCAL*F(1:M, 1:N)
                            SCALE = SCALE * SCAL
                            INFO = 1
                        END IF
                        EXIT
                    END IF

                END DO


                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    E(K,L) = RHS(1)
                    F(K,L) = RHS(2)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN

                    E(K, L)  = RHS(1)
                    E(K, LH)  = RHS(2)
                    F(K, L)  = RHS(3)
                    F(K, LH)  = RHS(4)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    E(K, L)  = RHS(1)
                    E(KH, L)  = RHS(2)
                    F(K, L)  = RHS(3)
                    F(KH, L)  = RHS(4)
                ELSE
                    E(K,L) = RHS(1)
                    E(KH,L)  = RHS(2)
                    E(K,LH)  = RHS(3)
                    E(KH,LH)  = RHS(4)

                    F(K,L) = RHS(5)
                    F(KH,L)  = RHS(6)
                    F(K,LH)  = RHS(7)
                    F(KH,LH)  = RHS(8)
                END IF


                ! Update January 2021
                IF ( L .GT. 0 ) THEN
                    CALL DGEMM('N','T', KB, LE, LB, -SGN1*ONE, F(K,L), LDF, B(1,L), LDB, ONE, E(K,1), LDE)
                    CALL DGEMM('N','T', KB, LE, LB, -SGN2*ONE, F(K,L), LDF, D(1,L), LDD, ONE, F(K,1), LDF)
                END IF

                L = L - 1
            END DO
            ! Update January 2021
            IF ( K .GT. 1 ) THEN
                CALL DGEMM('N', 'N', KE, N, KB, -ONE, A(1,K), LDA, E(K,1), LDE, ONE, E(1,1), LDE)
                CALL DGEMM('N', 'N', KE, N, KB, -ONE, C(1,K), LDC, E(K,1), LDE, ONE, F(1,1), LDF)
            END IF

            K = K - 1
        END DO

    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
        ! (T, N)
        K = 1
        DO WHILE ( K .LE. M )
            KB = 1
            IF ( K .LT. M ) THEN
                IF ( ( A(K+1,K) .NE. ZERO )) THEN
                    KB = 2
                END IF
            END IF
            KH = K + KB - 1

            L = 1
            DO WHILE ( L .LE. N )
                LB = 1
                IF ( L .LT. N ) THEN
                    IF ( ( B(L+1,L) .NE. ZERO .OR. D(L+1,L) .NE. ZERO ) ) THEN
                        LB = 2
                    END IF
                END IF
                LH = L + LB - 1

                IT = 0
                DO

                    ! Solve Inner System
                    IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 2
                        MAT(1,1) = A(K,K)
                        MAT(1,2) = SGN1 * B(L,L)
                        MAT(2,1) = C(K,K)
                        MAT(2,2) = SGN2 * D(L,L)
                        RHS(1) = E(K,L)
                        RHS(2) = F(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 4
                        MAT(1:4,1:4) = ZERO

                        MAT(1,1) = A(K,K)
                        MAT(1,3) = SGN1 * B(L,L)
                        MAT(1,4) = SGN1 * B(LH, L)

                        MAT(2,2) = A(K,K)
                        MAT(2,3) = SGN1 * B(L,LH)
                        MAT(2,4) = SGN1 * B(LH, LH)

                        MAT(3,1) = C(K,K)
                        MAT(3,3) = SGN2 * D(L,L)
                        MAT(3,4) = SGN2 * D(LH, L)

                        MAT(4,2) = C(K,K)
                        MAT(4,3) = SGN2 * D(L,LH)
                        MAT(4,4) = SGN2 * D(LH, LH)

                        RHS(1) = E(K, L)
                        RHS(2) = E(K, LH)
                        RHS(3) = F(K, L)
                        RHS(4) = F(K, LH)

                    ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 4
                        MAT(1:4,1:4) = ZERO

                        MAT(1,1) = A(K,K)
                        MAT(1,2) = A(KH,K)
                        MAT(1,3) = SGN1 * B(L,L)

                        MAT(2,1) = A(K,KH)
                        MAT(2,2) = A(KH, KH)
                        MAT(2,4) = SGN1 * B(L,L)

                        MAT(3,1) = C(K,K)
                        MAT(3,3) = SGN2 * D(L,L)

                        MAT(4,1) = C(K, KH)
                        MAT(4,2) = C(KH, KH)
                        MAT(4,4) = SGN2 * D(L,L)

                        RHS(1) = E(K, L)
                        RHS(2) = E(KH, L)
                        RHS(3) = F(K, L)
                        RHS(4) = F(KH, L)

                    ELSE
                        DIMMAT = 8
                        MAT(1:8,1:8) = ZERO

                        MAT(1,1) = A(K,K)
                        MAT(1,2) = A(KH,K)
                        MAT(1,5) = SGN1 * B(L,L)
                        MAT(1,7) = SGN1 * B(LH, L)

                        MAT(2,1) = A(K,KH)
                        MAT(2,2) = A(KH,KH)
                        MAT(2,6) = SGN1 * B(L,L)
                        MAT(2,8) = SGN1 * B(LH, L)

                        MAT(3,3) = A(K,K)
                        MAT(3,4) = A(KH,K)
                        MAT(3,5) = SGN1 * B(L,LH)
                        MAT(3,7) = SGN1 * B(LH,LH)

                        MAT(4,3) = A(K,KH)
                        MAT(4,4) = A(KH,KH)
                        MAT(4,6) = SGN1 * B(L,LH)
                        MAT(4,8) = SGN1 * B(LH,LH)

                        MAT(5,1) = C(K,K)
                        MAT(5,5) = SGN2 * D(L,L)
                        MAT(5,7) = SGN2 * D(LH, L)

                        MAT(6,1) = C(K,KH)
                        MAT(6,2) = C(KH,KH)
                        MAT(6,6) = SGN2 * D(L,L)
                        MAT(6,8) = SGN2 * D(LH, L)

                        MAT(7,3) = C(K,K)
                        MAT(7,5) = SGN2 * D(L,LH)
                        MAT(7,7) = SGN2 * D(LH,LH)

                        MAT(8,3) = C(K,KH)
                        MAT(8,4) = C(KH,KH)
                        MAT(8,6) = SGN2 * D(L,LH)
                        MAT(8,8) = SGN2 * D(LH,LH)


                        RHS(1) = E(K,L)
                        RHS(2) = E(KH,L)
                        RHS(3) = E(K,LH)
                        RHS(4) = E(KH,LH)

                        RHS(5) = F(K,L)
                        RHS(6) = F(KH,L)
                        RHS(7) = F(K,LH)
                        RHS(8) = F(KH,LH)

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
                            E(1:M, 1:N) = SCAL*E(1:M, 1:N)
                            F(1:M, 1:N) = SCAL*F(1:M, 1:N)
                            SCALE = SCALE * SCAL
                            INFO = 1
                        END IF
                        EXIT
                    END IF

                END DO

                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    E(K,L) = RHS(1)
                    F(K,L) = RHS(2)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN

                    E(K, L)  = RHS(1)
                    E(K, LH)  = RHS(2)
                    F(K, L)  = RHS(3)
                    F(K, LH)  = RHS(4)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    E(K, L)  = RHS(1)
                    E(KH, L)  = RHS(2)
                    F(K, L)  = RHS(3)
                    F(KH, L)  = RHS(4)
                ELSE
                    E(K,L) = RHS(1)
                    E(KH,L)  = RHS(2)
                    E(K,LH)  = RHS(3)
                    E(KH,LH)  = RHS(4)

                    F(K,L) = RHS(5)
                    F(KH,L)  = RHS(6)
                    F(K,LH)  = RHS(7)
                    F(KH,LH)  = RHS(8)
                END IF


                ! in current row
                IF ( LH .LT. N ) THEN
                    CALL DGEMM('N','N', KB, N-LH, LB, -SGN1*ONE, F(K,L), LDF, B(L,LH+1), LDB, ONE, E(K,LH+1), LDE)
                    CALL DGEMM('N','N', KB, N-LH, LB, -SGN2*ONE, F(K,L), LDF, D(L,LH+1), LDD, ONE, F(K,LH+1), LDF)
                END IF

                L = LH + 1
            END DO

            IF ( KH .LT. M ) THEN
                CALL DGEMM('T', 'N', M-KH, N, KB, -ONE, A(K,KH+1), LDA, E(K,1), LDE, ONE, E(KH+1,1), LDE)
                CALL DGEMM('T', 'N', M-KH, N, KB, -ONE, C(K,KH+1), LDC, E(K,1), LDE, ONE, F(KH+1,1), LDF)
            END IF

            K = KH + 1
        END DO
        !

    ELSE IF ( BTRANSA .AND. BTRANSB) THEN
        ! (T, T)
        K = 1
        DO WHILE ( K .LE. M )
            KB = 1
            IF ( K .LT. M ) THEN
                IF ( ( A(K+1,K) .NE. ZERO )) THEN
                    KB = 2
                END IF
            END IF
            KH = K + KB - 1

            L = N
            DO WHILE ( L .GT. 0 )
                LB = 1
                LH = L
                IF ( L .GT. 1 ) THEN
                    IF (  (B(L,L-1) .NE. ZERO .OR. D(L,L-1) .NE. ZERO) ) THEN
                        LB = 2
                        L = L - 1
                    END IF
                END IF
                LE = L - 1

                IT = 0
                DO
                    ! Solve Inner System
                    IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 2
                        MAT(1,1) = A(K,K)
                        MAT(1,2) = SGN1 * B(L,L)
                        MAT(2,1) = C(K,K)
                        MAT(2,2) = SGN2 * D(L,L)
                        RHS(1) = E(K,L)
                        RHS(2) = F(K,L)

                    ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                        DIMMAT = 4
                        MAT(1:4,1:4) = ZERO

                        MAT(1,1) = A(K,K)
                        MAT(1,3) = SGN1 * B(L,L)
                        MAT(1,4) = SGN1 * B(L, LH)

                        MAT(2,2) = A(K,K)
                        MAT(2,3) = SGN1 * B(LH,L)
                        MAT(2,4) = SGN1 * B(LH, LH)

                        MAT(3,1) = C(K,K)
                        MAT(3,3) = SGN2 * D(L,L)
                        MAT(3,4) = SGN2 * D(L, LH)

                        MAT(4,2) = C(K,K)
                        MAT(4,3) = SGN2 * D(LH,L)
                        MAT(4,4) = SGN2 * D(LH, LH)

                        RHS(1) = E(K, L)
                        RHS(2) = E(K, LH)
                        RHS(3) = F(K, L)
                        RHS(4) = F(K, LH)

                    ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                        DIMMAT = 4
                        MAT(1:4,1:4) = ZERO

                        MAT(1,1) = A(K,K)
                        MAT(1,2) = A(KH,K)
                        MAT(1,3) = SGN1 * B(L,L)

                        MAT(2,1) = A(K,KH)
                        MAT(2,2) = A(KH, KH)
                        MAT(2,4) = SGN1 * B(L,L)

                        MAT(3,1) = C(K,K)
                        MAT(3,3) = SGN2 * D(L,L)

                        MAT(4,1) = C(K, KH)
                        MAT(4,2) = C(KH, KH)
                        MAT(4,4) = SGN2 * D(L,L)

                        RHS(1) = E(K, L)
                        RHS(2) = E(KH, L)
                        RHS(3) = F(K, L)
                        RHS(4) = F(KH, L)

                    ELSE
                        DIMMAT = 8
                        MAT(1:8,1:8) = ZERO

                        MAT(1,1) = A(K,K)
                        MAT(1,2) = A(KH,K)
                        MAT(1,5) = SGN1 * B(L,L)
                        MAT(1,7) = SGN1 * B(L, LH)

                        MAT(2,1) = A(K,KH)
                        MAT(2,2) = A(KH,KH)
                        MAT(2,6) = SGN1 * B(L,L)
                        MAT(2,8) = SGN1 * B(L, LH)

                        MAT(3,3) = A(K,K)
                        MAT(3,4) = A(KH,K)
                        MAT(3,5) = SGN1 * B(LH,L)
                        MAT(3,7) = SGN1 * B(LH,LH)

                        MAT(4,3) = A(K,KH)
                        MAT(4,4) = A(KH,KH)
                        MAT(4,6) = SGN1 * B(LH,L)
                        MAT(4,8) = SGN1 * B(LH,LH)

                        MAT(5,1) = C(K,K)
                        MAT(5,5) = SGN2 * D(L,L)
                        MAT(5,7) = SGN2 * D(L, LH)

                        MAT(6,1) = C(K,KH)
                        MAT(6,2) = C(KH,KH)
                        MAT(6,6) = SGN2 * D(L,L)
                        MAT(6,8) = SGN2 * D(L, LH)

                        MAT(7,3) = C(K,K)
                        MAT(7,5) = SGN2 * D(LH,L)
                        MAT(7,7) = SGN2 * D(LH,LH)

                        MAT(8,3) = C(K,KH)
                        MAT(8,4) = C(KH,KH)
                        MAT(8,6) = SGN2 * D(LH,L)
                        MAT(8,8) = SGN2 * D(LH,LH)


                        RHS(1) = E(K,L)
                        RHS(2) = E(KH,L)
                        RHS(3) = E(K,LH)
                        RHS(4) = E(KH,LH)

                        RHS(5) = F(K,L)
                        RHS(6) = F(KH,L)
                        RHS(7) = F(K,LH)
                        RHS(8) = F(KH,LH)

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
                            E(1:M, 1:N) = SCAL*E(1:M, 1:N)
                            F(1:M, 1:N) = SCAL*F(1:M, 1:N)

                            SCALE = SCALE * SCAL
                            INFO = 1
                        END IF
                        EXIT
                    END IF

                END DO

                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    E(K,L) = RHS(1)
                    F(K,L) = RHS(2)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN

                    E(K, L)  = RHS(1)
                    E(K, LH)  = RHS(2)
                    F(K, L)  = RHS(3)
                    F(K, LH)  = RHS(4)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    E(K, L)  = RHS(1)
                    E(KH, L)  = RHS(2)
                    F(K, L)  = RHS(3)
                    F(KH, L)  = RHS(4)
                ELSE
                    E(K,L) = RHS(1)
                    E(KH,L)  = RHS(2)
                    E(K,LH)  = RHS(3)
                    E(KH,LH)  = RHS(4)

                    F(K,L) = RHS(5)
                    F(KH,L)  = RHS(6)
                    F(K,LH)  = RHS(7)
                    F(KH,LH)  = RHS(8)
                END IF

                ! in current row
                IF ( L .GT. 1  ) THEN
                    CALL DGEMM('N','T', KB, LE, LB, -SGN1*ONE, F(K,L), LDF, B(1,L), LDB, ONE, E(K,1), LDE)
                    CALL DGEMM('N','T', KB, LE, LB, -SGN2*ONE, F(K,L), LDF, D(1,L), LDD, ONE, F(K,1), LDF)
                END IF

                L = L - 1
            END DO

            ! Update January 2021
            IF ( KH .LT. M ) THEN
                CALL DGEMM('T', 'N', M-KH, N, KB, -ONE, A(K,KH+1), LDA, E(K,1), LDE, ONE, E(KH+1,1), LDE)
                CALL DGEMM('T', 'N', M-KH, N, KB, -ONE, C(K,KH+1), LDC, E(K,1), LDE, ONE, F(KH+1,1), LDF)
            END IF

            K = KH + 1
        END DO

    END IF
END SUBROUTINE



