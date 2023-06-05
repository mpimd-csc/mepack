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


!> \brief Level-3 Bartels-Stewart Algorithm for the generalized Sylvester equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TGSYLV_L3_2S ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANSA(1), TRANSB(1)
!           REAL SGN, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
!           REAL A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
!           REAL WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TGSYLV_L3_2S solves a generalized Sylvester equation of the following forms
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
!> created by SGGES from LAPACK.
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
!>          The leading dimension of the array X.  LDX >= max(1,M).
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
!
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
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date June 2023
!> \ingroup sgltgsylv
!
SUBROUTINE SLA_TGSYLV_L3_2S ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE MEPACK_OPTIONS_ISOLVER
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    REAL SGN, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
    REAL A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
    REAL WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE
    INTEGER INFO1
    REAL SCAL
    INTEGER MB, NB, ISOLVER, MAXNB
    INTEGER BIGNB
    INTEGER DAG_INFO

    LOGICAL BTRANSA, BTRANSB
    REAL ZERO , ONE
    PARAMETER (ZERO = 0.0, ONE = 1.0)
    INTEGER IZERO, IONE
    PARAMETER(IZERO = 0, IONE = 1 )

    ! EXTERNAL FUNCTIONS
    EXTERNAL SGEMM
    EXTERNAL SLACPY
    EXTERNAL SLA_STRMM_FIX
    EXTERNAL SLA_TGSYLV_DAG
    EXTERNAL SLA_TGSYLV_L2
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_128
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_32
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_64
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_96
    EXTERNAL SLA_TGSYLV_L2_REORDER
    EXTERNAL SLA_TGSYLV_L2_UNOPT
    EXTERNAL SLA_TGSYLV_RECURSIVE
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME

    MB = TGSYLV_BLOCKSIZE_MB(M,N)
    NB = TGSYLV_BLOCKSIZE_NB(M,N)
    ISOLVER = TGSYLV_ISOLVER()
    MAXNB = MAX(MB, NB)
    BIGNB = TGSYLV_BLOCKSIZE_2STAGE(M,N)


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
    ELSE IF ( BIGNB .LT. MAXNB) THEN
        INFO = -33
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_TGSYLV_L3_2S', -INFO)
        RETURN
    END IF


    IF ( ININFO .EQ. -1 ) THEN
        ! MEMORY DAG
        DAG_INFO = -1
        CALL SLA_TGSYLV_DAG(TRANSA, TRANSB, ONE, BIGNB, BIGNB,  A, BIGNB, B, BIGNB, &
                        & C, BIGNB, D, BIGNB,  X, BIGNB, SCAL, WORK, DAG_INFO)
        INFO = MAX( DAG_INFO, BIGNB**2, M*BIGNB, 4*NB)
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
                CALL SLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 64) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 96) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 128) THEN
                CALL SLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            ELSE
                CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            END IF
        CASE(2)
            IF ( M .GT. 128 .OR. N .GT. 128) THEN
                CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            ELSE
                CALL SLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
            END IF
        CASE(3)
            CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
        CASE(4)
            CALL SLA_TGSYLV_L2(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
        CASE(5)
            CALL SLA_TGSYLV_L2_UNOPT(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  X, LDX, SCAL, WORK, INFO1)
        CASE(6)
            CALL SLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, KB, LB,  A, LDA, B, LDB, &
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
        L = 1
        LH = 1
        DO WHILE ( LH .LE. N )
            L = LH
            LH = MIN(N, L + BIGNB - 1)
            LB = LH - L + 1
            IF ( LH .LT. N ) THEN
                IF ( ( B(LH+1,LH) .NE. ZERO .OR. D(LH+1,LH) .NE. ZERO ) ) THEN
                    LB = LB - 1
                    LH = LH - 1
                END IF
            END IF

            K = M
            DO WHILE (K .GT. 0)
                IF (K .EQ. M .AND. MOD(M,BIGNB) .NE. IZERO) THEN
                    KH = K
                    K = (M/BIGNB)*BIGNB + 1
                ELSE
                    KH = K
                    K = MAX(1, KH - BIGNB + 1)
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

                CALL SLA_TGSYLV_DAG(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)


                IF ( INFO1 .LT. 0 ) THEN
                    INFO = 100
                END IF

                IF ( SCAL .NE. ONE ) THEN
                    CALL SLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:N) = X(1:M, 1:N) * SCAL
                    SCALE = SCALE * SCAL
                    CALL SLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF

                ! Update January 2021
                ! in current row
                IF ( K .GT. 1 ) THEN
                    CALL SGEMM('N', 'N', KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK, KB)
                    CALL SGEMM('N', 'N', K-1, LB, KB, -ONE, A(1,K), LDA, WORK, KB, ONE, X(1,L), LDX)

                    CALL SGEMM('N', 'N', KB, LB, LB, ONE, X(K,L), LDX, D(L,L), LDD, ZERO, WORK, KB)
                    CALL SGEMM('N', 'N', K-1, LB, KB, -SGN*ONE, C(1,K), LDC, WORK, KB, ONE, X(1,L), LDX)
                END IF

                K = K - 1
            END DO

            ! Update January 2021
            IF ( LH .LT. N  ) THEN
                ! CALL SLACPY('All', M, LB, X(1,L), LDX,  WORK, M)
                ! CALL STRMM('Left', 'Upper', 'NoTrans', 'NoUnit', M, LB, ONE, A, LDA, WORK, M)
                ! CALL SLA_QTRMM2('Left', 'NoTrans', M, LB, ONE, A, LDA, X(1,L), LDX, WORK, M)
                CALL SLA_STRMM_FIX("N", M, LB, ONE, A, LDA, X(1,L), LDX, ZERO, WORK, M)

                CALL SGEMM("N", "N", M, N-LH, LB, -ONE, WORK, M, B(L,LH+1), LDB, ONE, X(1,LH+1), LDX)

                ! CALL SLACPY('All', M, LB, X(1,L), LDX,  WORK, M)
                ! CALL STRMM('Left', 'Upper', 'NoTrans', 'NoUnit', M, LB, ONE, C, LDC, WORK, M)
                CALL SLA_STRMM_FIX("N", M, LB, ONE, C, LDC, X(1,L), LDX, ZERO, WORK, M)

                CALL SGEMM("N", "N", M, N-LH, LB, -SGN*ONE, WORK, M, D(L,LH+1), LDD, ONE, X(1,LH+1), LDX)

            END IF

            LH = LH + 1
        END DO

    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N,T)
        L = N
        DO WHILE ( L .GT. 0 )
            IF (L .EQ. N .AND. MOD(N,BIGNB) .NE. IZERO) THEN
                LH = L
                L = (L/BIGNB)*BIGNB + 1
            ELSE
                LH = L
                L = MAX(1, LH - BIGNB + 1)
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


            K = M
            DO WHILE (K .GT. 0)
                IF (K .EQ. M .AND. MOD(M,BIGNB) .NE. IZERO) THEN
                    KH = K
                    K = (M/BIGNB)*BIGNB + 1
                ELSE
                    KH = K
                    K = MAX(1, KH - BIGNB + 1)
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

                CALL SLA_TGSYLV_DAG(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)

                IF ( INFO1 .LT. 0 ) THEN
                    INFO = 100
                END IF

                IF ( SCAL .NE. ONE ) THEN
                    CALL SLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:N) = X(1:M, 1:N) * SCAL
                    SCALE = SCALE * SCAL
                    CALL SLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF

                ! Update January 2021
                IF ( K .GT. 1 ) THEN
                    CALL SGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK, KB)
                    CALL SGEMM('N', 'N', KE, LB, KB, -ONE, A(1,K), LDA, WORK, KB, ONE, X(1,L), LDX)

                    CALL SGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, D(L,L), LDD, ZERO, WORK, KB)
                    CALL SGEMM('N', 'N', KE, LB, KB, -SGN*ONE, C(1,K), LDC, WORK, KB, ONE, X(1,L), LDX)

                END IF

                K = K - 1
            END DO

            ! Update January 2021
            IF ( L .GT. 1 ) THEN
                ! CALL SLACPY('All', M, LB, X(1,L), LDX, WORK, M)
                ! CALL STRMM('Left', 'Upper', 'NoTrans', 'NoUnit', M, LB, ONE, A, LDA, WORK, M)
                ! CALL SLA_QTRMM2('Left', 'NoTrans', M, LB, ONE, A, LDA, X(1,L), LDX, WORK, M)
                CALL SLA_STRMM_FIX("N", M, LB, ONE, A, LDA, X(1,L), LDX, ZERO, WORK, M)

                CALL SGEMM('N', 'T', M, LE, LB, -ONE, WORK, M, B(1,L), LDB, ONE, X(1,1), LDX)

                ! CALL SLACPY('All', M, LB, X(1,L), LDX, WORK, M)
                ! CALL STRMM('Left', 'Upper', 'NoTrans', 'NoUnit', M, LB, ONE, C, LDC, WORK, M)
                CALL SLA_STRMM_FIX("N", M, LB, ONE, C, LDC, X(1,L), LDX, ZERO, WORK, M)

                CALL SGEMM('N', 'T', M, LE, LB, -SGN*ONE, WORK, M, D(1,L), LDD, ONE, X(1,1), LDX)
            END IF

            L = L - 1
        END DO

    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
        ! (T, N)

        L = 1
        LH = 1
        DO WHILE ( LH .LE. N )
            L = LH
            LH = MIN(N, L + BIGNB - 1)
            LB = LH - L + 1
            IF ( LH .LT. N ) THEN
                IF ( ( B(LH+1,LH) .NE. ZERO .OR. D(LH+1,LH) .NE. ZERO ) ) THEN
                    IF ( LB.EQ.1) THEN
                        LB = LB + 1
                        LH = LH + 1
                    ELSE
                        LB = LB - 1
                        LH = LH - 1
                    END IF
                END IF
            END IF

            K = 1
            KH = 1
            DO WHILE ( KH .LE. M )
                K = KH
                KH = MIN(M, K + BIGNB - 1)
                KB = KH - K + 1
                IF ( KH .LT. M ) THEN
                    IF (  A(KH+1, KH) .NE. ZERO ) THEN
                        IF ( KB .EQ. 1) THEN
                            KB = KB + 1
                            KH = KH + 1
                        ELSE
                            KB = KB - 1
                            KH = KH - 1
                        END IF
                    END IF
                END IF


                CALL SLA_TGSYLV_DAG(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)

                IF ( INFO1 .LT. 0 ) THEN
                    INFO = 100
                END IF

                IF ( SCAL .NE. ONE ) THEN
                    CALL SLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:N) = X(1:M, 1:N) * SCAL
                    SCALE = SCALE * SCAL
                    CALL SLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF


                ! in current row
                IF ( KH  .LT. M ) THEN
                    CALL SGEMM("N", "N", KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK, KB)
                    CALL SGEMM("T", "N", M-KH, LB, KB, -ONE, A(K,KH+1), LDA, WORK, KB, ONE, X(KH+1, L), LDX)

                    CALL SGEMM("N", "N", KB, LB, LB, ONE, X(K,L), LDX, D(L,L), LDD, ZERO, WORK, KB)
                    CALL SGEMM("T", "N", M-KH, LB, KB, -SGN*ONE, C(K,KH+1), LDC, WORK, KB, ONE, X(KH+1, L), LDX)


                END IF

                KH = KH + 1

            END DO

            IF ( LH .LT. N ) THEN
                ! CALL SGEMM("T", "N", M, LB, M, ONE, A, LDA, X(1,L), LDX, ZERO, WORK, M)
                ! CALL SLACPY('All', M, LB, X(1,L), LDX,  WORK, M)
                ! CALL STRMM('Left', 'Upper', 'Trans', 'NoUnit', M, LB, ONE, A, LDA, WORK, M)
                ! CALL SLA_QTRMM2('Left', 'Trans', M, LB, ONE, A, LDA, X(1,L), LDX, WORK, M)

                CALL SLA_STRMM_FIX("T", M, LB, ONE, A, LDA, X(1,L), LDX, ZERO, WORK, M)

                CALL SGEMM("N", "N", M, N-LH, LB, -ONE, WORK, M, B(L,LH+1), LDB, ONE, X(1,LH+1), LDX)

                ! CALL SLACPY('All', M, LB, X(1,L), LDX,  WORK, M)
                ! CALL STRMM('Left', 'Upper', 'Trans', 'NoUnit', M, LB, ONE, C, LDC, WORK, M)
                CALL SLA_STRMM_FIX("T", M, LB, ONE, C, LDC, X(1,L), LDX, ZERO, WORK, M)

                ! CALL SGEMM("T", "N", M, LB, M, ONE, C, LDC, X(1,L), LDX, ZERO, WORK, M)
                CALL SGEMM("N", "N", M, N-LH, LB, -SGN*ONE, WORK, M, D(L,LH+1), LDD, ONE, X(1,LH+1), LDX)

            END IF

            LH = LH + 1
        END DO


    ELSE IF ( BTRANSA .AND. BTRANSB) THEN
        ! (T, T)
        L = N
        DO WHILE ( L .GT. 0 )
            IF (L .EQ. N .AND. MOD(N,BIGNB) .NE. IZERO) THEN
                LH = L
                L = (L/BIGNB)*BIGNB + 1
            ELSE
                LH = L
                L = MAX(1, LH - BIGNB + 1)
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

            K = 1
            KH = 1
            DO WHILE ( KH .LE. M )
                K = KH
                KH = MIN(M, K + BIGNB - 1)
                KB = KH - K + 1
                IF ( KH .LT. M ) THEN
                    IF ( A(KH+1, KH) .NE. ZERO ) THEN
                        IF ( KB .EQ. 1 ) THEN
                            KB = KB + 1
                            KH = KH + 1
                        ELSE
                            KB = KB - 1
                            KH = KH - 1
                        END IF
                    END IF
                END IF

                CALL SLA_TGSYLV_DAG(TRANSA, TRANSB, SGN, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & C(K,K), LDC, D(L,L), LDD,  X(K,L), LDX, SCAL, WORK, INFO1)

                IF ( INFO1 .LT. 0 ) THEN
                    INFO = 100
                END IF

                IF ( SCAL .NE. ONE ) THEN
                    CALL SLACPY("All", KB, LB, X(K,L), LDX, WORK, KB)
                    X(1:M, 1:N) = X(1:M, 1:N) * SCAL
                    SCALE = SCALE * SCAL
                    CALL SLACPY("All", KB, LB, WORK, KB, X(K,L), LDX)
                END IF


                IF ( KH .LT. M  ) THEN
                    CALL SGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, B(L,L), LDB, ZERO, WORK, KB)
                    CALL SGEMM('T', 'N', M-KH, LB, KB, -ONE, A(K, KH+1), LDA, WORK, KB, ONE, X(KH+1,L),LDX)

                    CALL SGEMM('N', 'T', KB, LB, LB, ONE, X(K,L), LDX, D(L,L), LDB, ZERO, WORK, KB)
                    CALL SGEMM('T', 'N', M-KH, LB, KB, -SGN*ONE, C(K, KH+1), LDC, WORK, KB, ONE, X(KH+1,L),LDX)

                END IF
                KH = KH + 1

            END DO

            ! Update January 2021
            IF ( L .GT. IONE ) THEN
                ! CALL SLACPY('All', M, LB, X(1,L), LDX,  WORK, M)
                ! CALL STRMM('Left', 'Upper', 'Trans', 'NoUnit', M, LB, ONE, A, LDA, WORK, M)
                ! CALL SLA_QTRMM2('Left','Trans',M, LB, ONE, A, LDA, X(1,L), LDX, WORK, M)
                CALL SLA_STRMM_FIX("T", M, LB, ONE, A, LDA, X(1,L), LDX, ZERO, WORK, M)

                CALL SGEMM('N', 'T', M, LE, LB, -ONE, WORK, M, B(1, L), LDB, ONE, X(1,1), LDX)

                ! CALL SLACPY('All', M, LB, X(1,L), LDX,  WORK, M)
                ! CALL STRMM('Left', 'Upper', 'Trans', 'NoUnit', M, LB, ONE, C, LDC, WORK, M)
                CALL SLA_STRMM_FIX("T", M, LB, ONE, C, LDC, X(1,L), LDX, ZERO, WORK, M)

                CALL SGEMM('N', 'T', M, LE, LB, -SGN*ONE, WORK, M, D(1, L), LDD, ONE, X(1,1), LDX)


            END IF

            L = L - 1
        END DO

    END IF
END SUBROUTINE



