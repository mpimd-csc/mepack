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


!> \brief Level-3 Bartels-Stewart Algorithm for the generalized Sylvester equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TGSYLV_L3_ROWWISE( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANSA(1), TRANSB(1)
!           DOUBLE PRECISION SGN, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
!           DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
!           DOUBLE PRECISION WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLA_TGSYLV_L3_ROWWISE solves a generalized Sylvester equation of the following forms
!>
!>    op1(A) * X * op2(B) + op1(C) * X * op2(D) = SCALE * Y                              (1)
!>
!> or
!>
!>    op1(A) * X * op2(B) - op1(C) * X * op2(D) = SCALE * Y                              (2)
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by DGGES from LAPACK.
!> \endverbatim
!> \attention The algorithm is implemented using BLAS level 3 operations.
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
!>          SGN is DOUBLE PRECISION, allowed values: +/-1
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
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,N)
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
!
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Dezember 2022
!> \ingroup dbltgsylv
!
SUBROUTINE DLA_TGSYLV_L3_ROWWISE ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    USE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
    USE MEPACK_OPTIONS_ISOLVER
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    DOUBLE PRECISION SGN, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
    DOUBLE PRECISION WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE
    INTEGER INFO1
    DOUBLE PRECISION SCAL
    INTEGER NB,MB, ISOLVER, MAXNB

    LOGICAL BTRANSA, BTRANSB
    DOUBLE PRECISION ZERO , ONE
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
    INTEGER IZERO
    PARAMETER(IZERO = 0)

    ! EXTERNAL FUNCTIONS
    EXTERNAL DGEMM
    EXTERNAL DLACPY
    EXTERNAL DLA_QTRMM2
    EXTERNAL DLA_TGSYLV_L2
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_128
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_32
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_64
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_96
    EXTERNAL DLA_TGSYLV_L2_REORDER
    EXTERNAL DLA_TGSYLV_L2_ROWWISE
    EXTERNAL DLA_TGSYLV_RECURSIVE
    EXTERNAL DTRMM
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME


    MB = TGSYLV_BLOCKSIZE_MB(M,N)
    NB = TGSYLV_BLOCKSIZE_NB(M,N)
    ISOLVER = TGSYLV_ISOLVER()
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
    ELSE IF ( LDC .LT. MAX(1, M)) THEN
        INFO = -11
    ELSE IF ( LDD .LT. MAX(1, N)) THEN
        INFO = -13
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -15
    ELSE IF ( MB .LT. 2 ) THEN
        INFO = -30
    ELSE IF ( NB .LT. 2 ) THEN
        INFO = -31
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('DLA_TGSYLV_L3_ROWWISE', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1 ) THEN
        INFO =  MAX(2*MAX(NB,MB),MAX(MB,NB)**2, M*MAX(NB,MB))
        RETURN
    END IF



    ! Quick Return
    SCALE = ONE
    SCAL = ONE
    INFO = 0
    INFO1 = 0

    IF ( M .EQ. 0 .OR. N .EQ. 0 ) THEN
        RETURN
    END IF

    IF ( M .LE. MB .AND. N .LE. NB ) THEN
        KB = M
        LB = N
        SELECT CASE(ISOLVER)
        CASE(1)
            IF ( MAXNB .LE. 32 ) THEN
                CALL DLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 64) THEN
                CALL DLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 96) THEN
                CALL DLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 128) THEN
                CALL DLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            ELSE
                CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            END IF
        CASE(2)
            IF ( M .GT. 128 .OR. N .GT. 128) THEN
                CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            ELSE
                CALL DLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            END IF
        CASE(3)
            CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
        CASE(4)
            CALL DLA_TGSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
        CASE(5)
            CALL DLA_TGSYLV_L2_ROWWISE(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
        CASE(6)
            CALL DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
        END SELECT
        IF ( INFO1.NE. 0 ) THEN
            INFO = INFO1
        END IF
        IF ( SCAL .NE. ONE) THEN
            SCALE = SCAL
        END IF

        RETURN
    END IF

    ! Prepare
    INFO = 0

    IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB ) THEN
        ! ( N, N )
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
                IF (A(K,K-1) .NE. ZERO ) THEN
                    IF ( KB.EQ.1) THEN
                        KB = KB + 1
                        K = K - 1
                    ELSE
                        KB = KB - 1
                        K = K + 1
                    END IF
                END IF
            END IF
            KE = K - 1


            L = 1
            LH = 1
            DO WHILE ( LH .LE. N )
                L = LH
                LH = MIN(N, L + NB - 1)
                LB = LH - L + 1
                IF ( LH .LT. N ) THEN
                    IF ( B(LH+1,LH) .NE. ZERO .OR. D(LH+1,LH) .NE. ZERO )  THEN
                        LB = LB - 1
                        LH = LH - 1
                    END IF
                END IF

                SELECT CASE(ISOLVER)
                CASE(1)
                    IF ( MAXNB .LE. 32 ) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 64) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 96) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 128) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    END IF
                CASE(2)
                    IF ( LB.GT.128 .OR. KB.GT.128) THEN
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE
                        CALL DLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    END IF
                CASE(3)
                    CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(4)
                    CALL DLA_TGSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(5)
                    CALL DLA_TGSYLV_L2_ROWWISE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(6)
                    CALL DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                END SELECT

                IF ( INFO1 .LT. 0 ) THEN
                    INFO = 100
                END IF

                IF ( SCAL .NE. ONE ) THEN
                    CALL DLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:N) = X(1:M, 1:N) * SCAL
                    SCALE = SCALE * SCAL
                    CALL DLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF

                ! Update January 2021
                ! in current row
                IF ( LH .LT. N ) THEN
                    CALL DGEMM('N','N', KB, LB, KB, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, WORK, KB )
                    CALL DGEMM('N','N', KB, N-LH, LB, -ONE, WORK, KB, B(L,LH+1), LDB, ONE, X(K,LH+1), LDX)

                    CALL DGEMM('N','N', KB, LB, KB, ONE, C(K,K), LDC, X(K,L), LDX, ZERO, WORK, KB )
                    CALL DGEMM('N','N', KB, N-LH, LB, -ONE*SGN, WORK, KB, D(L,LH+1), LDD, ONE, X(K,LH+1), LDX)

                END IF

                LH = LH + 1
            END DO

            ! Update January 2021
            IF ( K .GT. 1 ) THEN
                CALL DLACPY('All', KB, N, X(K,1), LDX,  WORK, KB)
                CALL DTRMM('Right', 'Upper', 'NoTrans', 'NoUnit', KB, N, ONE, B, LDB, WORK, KB)
                CALL DLA_QTRMM2('Right','NoTrans',KB, N, ONE, B, LDB, X(K,1), LDX, WORK, KB)
                CALL DGEMM('N', 'N', KE, N, KB, -ONE, A(1,K), LDA, WORK, KB, ONE, X(1,1), LDX)

                CALL DLACPY('All', KB, N, X(K,1), LDX, WORK, KB)
                CALL DTRMM('Right', 'Upper', 'NoTrans', 'NoUnit', KB, N, ONE, D, LDD, WORK, KB)
                CALL DLA_QTRMM2('Right','NoTrans',KB, N, ONE, D, LDD, X(K,1), LDX, WORK, KB)
                CALL DGEMM('N', 'N', KE, N, KB, -ONE*SGN, C(1,K), LDC, WORK, KB, ONE, X(1,1), LDX)
            END IF

            K = K - 1
        END DO

    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N,T)
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
            KE = K - 1


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
                    IF ( (B(L,L-1) .NE. ZERO .OR. D(L,L-1) .NE. ZERO) ) THEN
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

                SELECT CASE(ISOLVER)
                CASE(1)
                    IF ( MAXNB .LE. 32 ) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 64) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 96) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 128) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    END IF
                CASE(2)
                    IF ( LB.GT.128 .OR. KB.GT.128) THEN
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE
                        CALL DLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    END IF
                CASE(3)
                    CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(4)
                    CALL DLA_TGSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(5)
                    CALL DLA_TGSYLV_L2_ROWWISE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(6)
                    CALL DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                END SELECT

                IF ( INFO1 .LT. 0 ) THEN
                    INFO = 100
                END IF
                IF ( SCAL .NE. ONE ) THEN
                    CALL DLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:N) = X(1:M, 1:N) * SCAL
                    SCALE = SCALE * SCAL
                    CALL DLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF


                ! Update January 2021
                IF ( L .GT. 0 ) THEN
                    CALL DGEMM('N','N', KB, LB, KB, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, WORK, KB )
                    CALL DGEMM('N','T', KB, LE, LB, -ONE, WORK, KB, B(1,L), LDB, ONE, X(K,1), LDX)

                    CALL DGEMM('N','N', KB, LB, KB, ONE, C(K,K), LDC, X(K,L), LDX, ZERO, WORK, KB )
                    CALL DGEMM('N','T', KB, LE, LB, -ONE*SGN, WORK, KB, D(1,L), LDD, ONE, X(K,1), LDX)
                END IF

                L = L - 1
            END DO

            ! Update January 2021
            IF ( K .GT. 1 ) THEN
                CALL DLACPY('All', KB, N, X(K,1), LDX,  WORK, KB)
                CALL DTRMM('Right', 'Upper', 'Trans', 'NoUnit', KB, N, ONE, B, LDB, WORK, KB)
                CALL DLA_QTRMM2('Right','Trans',KB, N, ONE, B, LDB, X(K,1), LDX, WORK, KB)
                CALL DGEMM('N', 'N', KE, N, KB, -ONE, A(1,K), LDA, WORK, KB, ONE, X(1,1), LDX)

                CALL DLACPY('All', KB, N, X(K,1), LDX, WORK, KB)
                CALL DTRMM('Right', 'Upper', 'Trans', 'NoUnit', KB, N, ONE, D, LDD, WORK, KB)
                CALL DLA_QTRMM2('Right','Trans',KB, N, ONE, D, LDD, X(K,1), LDX, WORK, KB)
                CALL DGEMM('N', 'N', KE, N, KB, -ONE*SGN, C(1,K), LDC, WORK, KB, ONE, X(1,1), LDX)
            END IF

            K = K - 1
        END DO

    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
        ! (T, N)
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


            L = 1
            LH = 1
            DO WHILE ( LH .LE. N )
                L = LH
                LH = MIN(N, L + NB - 1)
                LB = LH - L + 1
                IF ( LH .LT. N ) THEN
                    IF ( ( B(LH+1,LH) .NE. ZERO .OR. D(LH+1,LH) .NE. ZERO ) ) THEN
                        LB = LB - 1
                        LH = LH - 1
                    END IF
                END IF


                SELECT CASE(ISOLVER)
                CASE(1)
                    IF ( MAXNB .LE. 32 ) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 64) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 96) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 128) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    END IF
                CASE(2)
                    IF ( LB.GT.128 .OR. KB.GT.128) THEN
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE
                        CALL DLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    END IF

                CASE(3)
                    CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(4)
                    CALL DLA_TGSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(5)
                    CALL DLA_TGSYLV_L2_ROWWISE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(6)
                    CALL DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                END SELECT

                IF ( SCAL .NE. ONE ) THEN
                    CALL DLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:N) = X(1:M, 1:N) * SCAL
                    SCALE = SCALE * SCAL
                    CALL DLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF


                ! in current row
                IF ( LH .LT. N ) THEN
                    CALL DGEMM('T','N', KB, LB, KB, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, WORK, KB )
                    CALL DGEMM('N','N', KB, N-LH, LB, -ONE, WORK, KB, B(L,LH+1), LDB, ONE, X(K,LH+1), LDX)

                    CALL DGEMM('T','N', KB, LB, KB, ONE, C(K,K), LDC, X(K,L), LDX, ZERO, WORK, KB )
                    CALL DGEMM('N','N', KB, N-LH, LB, -ONE*SGN, WORK, KB, D(L,LH+1), LDD, ONE, X(K,LH+1), LDX)

                END IF


                LH = LH + 1
            END DO

            IF ( KH .LT. M ) THEN
                CALL DLACPY('All', KB, N, X(K,1), LDX,  WORK, KB)
                CALL DTRMM('Right', 'Upper', 'NoTrans', 'NoUnit', KB, N, ONE, B, LDB, WORK, KB)
                CALL DLA_QTRMM2('Right','NoTrans',KB, N, ONE, B, LDB, X(K,1), LDX, WORK, KB)
                CALL DGEMM('T', 'N', M-KH, N, KB, -ONE, A(K,KH+1), LDA, WORK, KB, ONE, X(KH+1,1), LDX)

                CALL DLACPY('All', KB, N, X(K,1), LDX,  WORK, KB)
                CALL DTRMM('Right', 'Upper', 'NoTrans', 'NoUnit', KB, N, ONE, D, LDD, WORK, KB)
                CALL DLA_QTRMM2('Right','NoTrans',KB, N, ONE, D, LDD, X(K,1), LDX, WORK, KB)
                CALL DGEMM('T', 'N', M-KH, N, KB, -ONE*SGN, C(K,KH+1), LDC, WORK, KB, ONE, X(KH+1,1), LDX)

            END IF

            KH = KH + 1
        END DO


    ELSE IF ( BTRANSA .AND. BTRANSB) THEN
        ! (T, T)
        K = 1
        KH = 1
        DO WHILE ( KH .LE. M )
            K = KH
            KH = MIN(M, K + MB - 1)
            KB = KH - K + 1
            IF ( KH .LT. M) THEN
                IF ( A(KH+1, KH) .NE. ZERO ) THEN
                    KB = KB - 1
                    KH = KH - 1
                END IF
            END IF


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
                    IF ( (B(L,L-1) .NE. ZERO .OR. D(L,L-1) .NE. ZERO) ) THEN
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

                SELECT CASE(ISOLVER)
                CASE(1)
                    IF ( MAXNB .LE. 32 ) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 64) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 96) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE IF ( MAXNB .LE. 128) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    END IF
                CASE(2)
                    IF ( LB.GT.128 .OR. KB.GT.128) THEN
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    ELSE
                        CALL DLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                            & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                    END IF
                CASE(3)
                    CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(4)
                    CALL DLA_TGSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(5)
                    CALL DLA_TGSYLV_L2_ROWWISE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                CASE(6)
                    CALL DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                        & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)
                END SELECT

                IF ( INFO1 .LT. 0 ) THEN
                    INFO = 100
                END IF
                IF ( SCAL .NE. ONE ) THEN
                    CALL DLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:N) = X(1:M, 1:N) * SCAL
                    SCALE = SCALE * SCAL
                    CALL DLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF



                ! in current row
                IF ( L .GT. 1  ) THEN
                    CALL DGEMM('T','N', KB, LB, KB, ONE, A(K,K), LDA, X(K,L), LDX, ZERO, WORK, KB )
                    CALL DGEMM('N','T', KB, LE, LB, -ONE, WORK, KB, B(1,L), LDB, ONE, X(K,1), LDX)

                    CALL DGEMM('T','N', KB, LB, KB, ONE, C(K,K), LDC, X(K,L), LDX, ZERO, WORK, KB )
                    CALL DGEMM('N','T', KB, LE, LB, -ONE*SGN, WORK, KB, D(1,L), LDD, ONE, X(K,1), LDX)

                END IF

                L = L - 1
            END DO

            ! Update January 2021
            IF ( KH .LT. M ) THEN
                CALL DLACPY('All', KB, N, X(K,1), LDX,  WORK, KB)
                CALL DTRMM('Right', 'Upper', 'Trans', 'NoUnit', KB, N, ONE, B, LDB, WORK, KB)
                CALL DLA_QTRMM2('Right','Trans',KB, N, ONE, B, LDB, X(K,1), LDX, WORK, KB)
                CALL DGEMM('T', 'N', M-KH, N, KB, -ONE, A(K,KH+1), LDA, WORK, KB, ONE, X(KH+1,1), LDX)

                CALL DLACPY('All', KB, N, X(K,1), LDX,  WORK, KB)
                CALL DTRMM('Right', 'Upper', 'Trans', 'NoUnit', KB, N, ONE, D, LDD, WORK, KB)
                CALL DLA_QTRMM2('Right','Trans',KB, N, ONE, D, LDD, X(K,1), LDX, WORK, KB)
                CALL DGEMM('T', 'N', M-KH, N, KB, -ONE*SGN, C(K,KH+1), LDC, WORK, KB, ONE, X(KH+1,1), LDX)

            END IF

            KH = KH + 1
        END DO

    END IF
END SUBROUTINE



