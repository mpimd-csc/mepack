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
!

!> \brief Level-3 Kagstroem-Westin/Bartels-Stewart Algorithm with sub blocking for the dual generalized coupled Sylvester equation
!
!  Definition:
!  ===========
!       SUBROUTINE DLA_TGCSYLV_DUAL_L3_2S ( TRANSA, TRANSB, SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
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
!> DLA_TGCSYLV_DUAL_L3_2S solves a generalized coupled  Sylvester equation of the following form
!>
!>    op1(A)**T * R  + op1(C)**T * L               = SCALE * E                             (1)
!>    SGN1 * R * op2(B)**T + SGN2 * L * op2(D)** T = SCALE * F
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand sides E, F and the solutions R, L are
!> M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by DGGES from LAPACK.
!> The equation (1) is the dual to the generalized coupled Sylvester equation
!>
!>    op1(A) * R + SGN1 * L * op2(B)  = SCALE * E                                          (2)
!>    op1(C) * R + SGN2 * L * op2(D)  = SCALE * F
!>
!> The equation (1) is the dual one to equation (2) with respect to the underlying linear system.
!> Let Z be the matrix formed by rewriting (2) into its Kronecker form. This yields
!>
!>         | kron(I, op1(A))   SGN1*kron(op2(B)**T, I) | | Vec R |   | Vec E |
!>   Z X = |                                           |*|       | = |       |
!>         | kron(I, op1(C))   SGN2*kron(op2(D)**T, I) | | Vec L |   | Vec F |
!>
!> Regarding Z**T one obtains
!>
!>            | kron(I, op1(A)**T )    kron(I, op1(C)**T)   | | Vec R |   | Vec E |
!>   Z**T X = |                                             |*|       | = |       |
!>            | SGN1*kron(op2(B), I)   SGN2*kron(op2(D), I) | | Vec L |   | Vec F |
!>
!> which belongs to the Sylvester equation (1). For this reason the parameters TRANSA and TRANSB
!> are expressed in terms of the Sylvester equation (2).
!>
!> \endverbatim
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
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Januar 2023
!> \ingroup dbltgsylv
!
SUBROUTINE DLA_TGCSYLV_DUAL_L3_2S ( TRANSA, TRANSB, SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
        & E, LDE, F, LDF,  SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    USE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
    USE MEPACK_OPTIONS_ISOLVER
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    DOUBLE PRECISION SGN1, SGN2, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF,*)
    DOUBLE PRECISION WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KE
    INTEGER INFO1
    DOUBLE PRECISION SCAL
    INTEGER NB,MB, ISOLVER, MAXNB, BIGNB

    LOGICAL BTRANSA, BTRANSB
    DOUBLE PRECISION ZERO , ONE
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
    INTEGER IZERO
    PARAMETER(IZERO = 0)

    ! EXTERNAL FUNCTIONS
    EXTERNAL DGEMM
    EXTERNAL DLA_TGCSYLV_DUAL_DAG
    EXTERNAL DLA_TGCSYLV_DUAL_L2
    EXTERNAL DLA_TGCSYLV_DUAL_L2_LOCAL_COPY
    EXTERNAL DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128
    EXTERNAL DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32
    EXTERNAL DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64
    EXTERNAL DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96
    EXTERNAL DLA_TGCSYLV_DUAL_RECURSIVE
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME


    MB = TGCSYLV_DUAL_BLOCKSIZE_MB(M,N)
    NB = TGCSYLV_DUAL_BLOCKSIZE_NB(M,N)
    ISOLVER = TGCSYLV_DUAL_ISOLVER()
    MAXNB = MAX(MB, NB)
    BIGNB = TGCSYLV_DUAL_BLOCKSIZE_2STAGE(M,N)

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
        CALL XERROR_HANDLER('DLA_TGCSYLV_DUAL_L3_2S', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1 ) THEN
        INFO = 1
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
                CALL DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 64) THEN
                CALL DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 96) THEN
                CALL DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 128) THEN
                CALL DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            ELSE
                CALL DLA_TGCSYLV_DUAL_L2(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            END IF
        CASE(2)
            IF ( M .GT. 128 .OR. N .GT. 128) THEN
                CALL DLA_TGCSYLV_DUAL_L2(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            ELSE
                CALL DLA_TGCSYLV_DUAL_L2_LOCAL_COPY(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            END IF
        CASE(3)
            CALL DLA_TGCSYLV_DUAL_L2(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
        CASE(4)
            CALL DLA_TGCSYLV_DUAL_L2(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
        CASE(5)
            CALL DLA_TGCSYLV_DUAL_L2(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
        CASE(6)
            CALL DLA_TGCSYLV_DUAL_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
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
        ! ( N, N ) - in the standard equation. This means
        ! ( T, T ) for the dual one.
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
                        KB = KB - 1
                        KH = KH - 1
                    END IF
                END IF

                CALL DLA_TGCSYLV_DUAL_DAG(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK, INFO1)

                IF ( KH .LT. M  ) THEN
                    CALL DGEMM('T', 'N', M-KH, LB, KB, -ONE, A(K, KH+1), LDA, E(K,L), LDE, ONE, E(KH+1,L),LDE)
                    CALL DGEMM('T', 'N', M-KH, LB, KB, -ONE, C(K, KH+1), LDC, F(K,L), LDF, ONE, E(KH+1,L),LDE)
                END IF
                KH = KH + 1

            END DO

            ! Update January 2021
            IF ( L .GT. 1) THEN
                CALL DGEMM('N', 'T', M, LE, LB, -ONE*SGN1, E(1,L), LDE, B(1, L), LDB, ONE, F(1,1), LDF)
                CALL DGEMM('N', 'T', M, LE, LB, -ONE*SGN2, F(1,L), LDF, D(1, L), LDD, ONE, F(1,1), LDF)
            END IF

            L = L - 1
        END DO


    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! ( N, T ) - in the standard equation. This means
        ! ( T, N ) for the dual one.
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

            K = 1
            KH = 1
            DO WHILE ( KH .LE. M )
                K = KH
                KH = MIN(M, K + BIGNB - 1)
                KB = KH - K + 1
                IF ( KH .LT. M ) THEN
                    IF ( A(KH+1, KH) .NE. ZERO ) THEN
                        KB = KB - 1
                        KH = KH - 1
                    END IF
                END IF

                CALL DLA_TGCSYLV_DUAL_DAG(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK, INFO1)

                ! in current row
                IF ( KH  .LT. M ) THEN
                    CALL DGEMM("T", "N", M-KH, LB, KB, -ONE, A(K,KH+1), LDA, E(K,L), LDE, ONE, E(KH+1, L), LDE)
                    CALL DGEMM("T", "N", M-KH, LB, KB, -ONE, C(K,KH+1), LDC, F(K,L), LDF, ONE, E(KH+1, L), LDE)
                END IF

                KH = KH + 1

            END DO

            IF ( LH .LT. N ) THEN
                CALL DGEMM("N", "N", M, N-LH, LB, -ONE*SGN1, E(1,L), LDE, B(L,LH+1), LDB, ONE, F(1,LH+1), LDF)
                CALL DGEMM("N", "N", M, N-LH, LB, -ONE*SGN2, F(1,L), LDF, D(L,LH+1), LDD, ONE, F(1,LH+1), LDF)
            END IF

            LH = LH + 1
        END DO


    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
        ! ( T, N ) - in the standard equation. This means
        ! ( N, T ) for the dual one.
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
                IF ( K .GT. 1) THEN
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

                CALL DLA_TGCSYLV_DUAL_DAG(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK, INFO1)

                ! Update January 2021
                IF ( K .GT. 1 ) THEN
                    CALL DGEMM('N', 'N', KE, LB, KB, -ONE, A(1,K), LDA, E(K,L), LDE, ONE, E(1,L), LDE)
                    CALL DGEMM('N', 'N', KE, LB, KB, -ONE, C(1,K), LDC, F(K,L), LDF, ONE, E(1,L), LDE)
                END IF

                K = K - 1
            END DO

            ! Update January 2021
            IF ( L .GT. 1 ) THEN
                CALL DGEMM('N', 'T', M, LE, LB, -ONE*SGN1, E(1,L), LDE, B(1,L), LDB, ONE, F(1,1), LDF)
                CALL DGEMM('N', 'T', M, LE, LB, -ONE*SGN2, F(1,L), LDF, D(1,L), LDD, ONE, F(1,1), LDF)
            END IF

            L = L - 1
        END DO


    ELSE IF ( BTRANSA .AND. BTRANSB) THEN
        ! ( T, T ) - in the standard equation. This means
        ! ( N, N ) for the dual one.
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

                CALL DLA_TGCSYLV_DUAL_DAG(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                    & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK, INFO1)



                ! Update January 2021
                ! in current row
                IF ( K .GT. 1 ) THEN
                    CALL DGEMM('N', 'N', K-1, LB, KB, -ONE, A(1,K), LDA, E(K,L), LDE, ONE, E(1,L), LDE)
                    CALL DGEMM('N', 'N', K-1, LB, KB, -ONE, C(1,K), LDC, F(K,L), LDF, ONE, E(1,L), LDE)
                END IF

                K = K - 1
            END DO

            ! Update January 2021
            IF ( LH .LT. N  ) THEN
                CALL DGEMM("N", "N", M, N-LH, LB, -ONE*SGN1, E(1,L), LDE, B(L,LH+1), LDB, ONE, F(1,LH+1), LDF)
                CALL DGEMM("N", "N", M, N-LH, LB, -ONE*SGN2, F(1,L), LDF, D(L,LH+1), LDD, ONE, F(1,LH+1), LDF)
            END IF

            LH = LH + 1
        END DO


    END IF
    RETURN
END SUBROUTINE




