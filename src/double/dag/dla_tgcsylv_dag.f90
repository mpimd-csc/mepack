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

!> \brief Level-3 Kagstroem-Westin/Bartels-Stewart Algorithm for the generalized coupled Sylvester equation with DAG scheduling
!
!  Definition:
!  ===========
!       SUBROUTINE DLA_TGCSYLV_DAG ( TRANSA, TRANSB, SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
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
!> DLA_TGCSYLV_DAG solves a generalized coupled  Sylvester equation of the following form
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
!> \remark The Algorithm uses OpenMP DAG scheduling
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
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Januar 2023
!> \ingroup dbltgsylv
!
SUBROUTINE DLA_TGCSYLV_DAG ( TRANSA, TRANSB, SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
        & E, LDE, F, LDF,  SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    DOUBLE PRECISION SGN1, SGN2, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF,*)
    DOUBLE PRECISION WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER K, L, KB, LB, KH, LE, LH, KOLD, LOLD
    INTEGER INFO1
    INTEGER J, JB, JH
    DOUBLE PRECISION SCAL
    INTEGER NB,MB, ISOLVER, MAXNB

    LOGICAL BTRANSA, BTRANSB
    DOUBLE PRECISION ZERO , ONE
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
    INTEGER IZERO,IONE
    PARAMETER(IZERO = 0, IONE = 1)

    ! EXTERNAL FUNCTIONS
    EXTERNAL DGEMM
    EXTERNAL DLA_TGCSYLV_L2
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_128
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_32
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_64
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_96
    EXTERNAL DLA_TGCSYLV_L2_REORDER
    EXTERNAL DLA_TGCSYLV_L2_UNOPT
    EXTERNAL DLA_TGCSYLV_RECURSIVE
    EXTERNAL IDLA_TGCSYLV_BLOCK
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME

    MB = TGCSYLV_BLOCKSIZE_MB(M,N)
    NB = TGCSYLV_BLOCKSIZE_NB(M,N)
    ISOLVER = TGCSYLV_ISOLVER()
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
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32

    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('DLA_TGCSYLV_DAG', -INFO)
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
                CALL DLA_TGCSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 64) THEN
                CALL DLA_TGCSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 96) THEN
                CALL DLA_TGCSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            ELSE IF ( MAXNB .LE. 128) THEN
                CALL DLA_TGCSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            ELSE
                CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            END IF
        CASE(2)
            IF ( M .GT. 128 .OR. N .GT. 128) THEN
                CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            ELSE
                CALL DLA_TGCSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
            END IF
        CASE(3)
            CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
        CASE(4)
            CALL DLA_TGCSYLV_L2(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
        CASE(5)
            CALL DLA_TGCSYLV_L2_UNOPT(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, INFO1)
        CASE(6)
            CALL DLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A, LDA, B, LDB, &
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
        ! ( N, N )

        !$omp parallel default(shared)
        !$omp master

        LOLD = 1
        KOLD = 2
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

                IF ( KH .EQ. M .AND. L .EQ. IONE ) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: F(K,L),E(K,L)) default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1 )

                    !$omp end task
                ELSE IF ( KH .EQ. M .AND. L .GT. IONE) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: F(K,L),E(K,L)) depend(in: F(K,LOLD),E(K,LOLD)) default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task
                ELSE IF ( KH .LT. M .AND. L .EQ. IONE) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: E(K,L),F(K,L)) depend(in: E(KOLD,L),F(KOLD,L))  default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task
                ELSE
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(in: E(KOLD,L), E(K,LOLD),F(KOLD,L),F(K,LOLD)) depend(out: E(K,L),F(K,L))  default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task
                END IF

                ! Update January 2021
                ! in current row
                IF ( K .GT. 1 ) THEN
                    J = K - 1
                    DO WHILE (J .GT. 0)
                        JH = J
                        J  = MAX(1, JH - MB + 1)
                        JB = JH - J + 1
                        IF ( J .GT. 1 ) THEN
                            IF ( A(J,J-1) .NE. ZERO ) THEN
                                JB = JB - 1
                                J = J + 1
                            END IF
                        END IF
                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:E(K,L)) depend(inout:E(J,L))
                        CALL DGEMM('N', 'N', JB , LB, KB, -ONE, A(J,K), LDA, E(K,L), LDE, ONE, E(J,L), LDE)
                        !$omp end task

                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:E(K,L)) depend(inout:F(J,L))
                        CALL DGEMM('N', 'N', JB , LB, KB, -ONE, C(J,K), LDC, E(K,L), LDE, ONE, F(J,L), LDF)
                        !$omp end task

                        J = J - 1
                    END DO
                END IF

                IF ( LH .LT. N  ) THEN
                    JH = LH+1
                    DO WHILE ( JH .LE. N )
                        J = JH
                        JH = MIN(N, J + NB - 1)
                        JB = JH - J + 1
                        IF ( JH .LT. N ) THEN
                            IF ( B(JH+1,JH) .NE. ZERO .OR. D(JH+1,JH) .NE. ZERO )  THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF

                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:F(K,L)) depend(inout:E(K,J))
                        CALL DGEMM("N", "N", KB, JB, LB, -ONE*SGN1, F(K,L), LDF, B(L,J), LDB, ONE, E(K,J), LDE)
                        !$omp end task

                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:F(K,L)) depend(inout:F(K,J))
                        CALL DGEMM("N", "N", KB, JB, LB, -ONE*SGN2, F(K,L), LDF, D(L,J), LDD, ONE, F(K,J), LDF)
                        !$omp end task

                        JH = JH + 1
                    END DO
                END IF
                KOLD = K
                K = K - 1
            END DO

            ! Update January 2021

            LOLD = L
            LH = LH + 1
        END DO
        !$omp end master
        !$omp taskwait
        !$omp end parallel

    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N,T)
        !$omp parallel default(shared)
        !$omp master
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
                IF (B(L,L-1) .NE. ZERO .OR. D(L,L-1) .NE. ZERO)  THEN
                    IF ( LB .EQ. 1 ) THEN
                        LB = LB + 1
                        L = L - 1
                    ELSE
                        LB = LB - 1
                        L = L + 1
                    END IF
                END IF
            END IF


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


                IF ( KH.EQ.M .AND. LH.EQ.N) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: F(K,L),E(K,L)) default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task
                ELSE IF ( KH.EQ.M .AND. LH .LT. N) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: F(K,L),E(K,L)) depend(in: F(K,LOLD),E(K,LOLD)) default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task

                ELSE IF ( KH.LT.M .AND. LH .EQ. N) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: E(K,L),F(K,L)) depend(in: E(KOLD,L),F(KOLD,L))  default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task

                ELSE
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(in: E(KOLD,L), E(K,LOLD),F(KOLD,L),F(K,LOLD)) depend(out: E(K,L),F(K,L))  default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task
                END IF




                ! Update January 2021
                IF ( K .GT. 1 ) THEN
                    J = K - 1
                    DO WHILE (J .GT. 0)
                        JH = J
                        J  = MAX(1, JH - MB + 1)
                        JB = JH - J + 1
                        IF ( J .GT. 1 ) THEN
                            IF (A(J,J-1) .NE. ZERO ) THEN
                                JB = JB - 1
                                J = J + 1
                            END IF
                        END IF

                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:E(K,L)) depend(inout:E(J,L))
                        CALL DGEMM('N', 'N', JB, LB, KB, -ONE, A(J,K), LDA, E(K,L), LDE, ONE, E(J,L), LDE)
                        !$omp end task

                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:E(K,L)) depend(inout:F(J,L))
                        CALL DGEMM('N', 'N', JB, LB, KB, -ONE, C(J,K), LDC, E(K,L), LDE, ONE, F(J,L), LDF)
                        !$omp end task

                        J = J - 1
                    END DO
                END IF

                IF ( L .GT. 1 ) THEN
                    J = L - 1
                    DO WHILE (J .GT. 0)
                        JH = J
                        J  = MAX(1, JH - NB + 1)
                        JB = JH - J + 1
                        IF ( J .GT. 1 ) THEN
                            IF ( (B(J,J-1) .NE. ZERO .OR. D(J,J-1).NE.ZERO) ) THEN
                                JB = JB - 1
                                J = J + 1
                            END IF
                        END IF
                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:F(K,L)) depend(inout:E(K,J))

                        CALL DGEMM('N', 'T', KB, JB, LB, -ONE*SGN1, F(K,L), LDF, B(J,L), LDB, ONE, E(K,J), LDE)
                        !$omp end task

                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:F(K,L)) depend(inout:F(K,J))

                        CALL DGEMM('N', 'T', KB, JB, LB, -ONE*SGN2, F(K,L), LDF, D(J,L), LDD, ONE, F(K,J), LDF)
                        !$omp end task

                        J = J - 1
                    END DO
                END IF

                KOLD = K
                K = K - 1
            END DO

            LOLD = L
            L = L - 1
        END DO
        !$omp end master
        !$omp taskwait
        !$omp end parallel


    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
        ! (T, N)

        !$omp parallel default(shared)
        !$omp master
        KOLD  =  1
        LOLD =  1
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

                IF ( K.EQ.IONE .AND. L.EQ.IONE) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: F(K,L),E(K,L)) default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)


                    !$omp end task
                ELSE IF ( K.EQ.IONE .AND. L.GT.IONE) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: F(K,L),E(K,L)) depend(in: F(K,LOLD),E(K,LOLD)) default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)


                    !$omp end task

                ELSE IF ( K.GT. IONE .AND. L.EQ.IONE) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: E(K,L),F(K,L)) depend(in: E(KOLD,L),F(KOLD,L))  default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task

                ELSE
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(in: E(KOLD,L), E(K,LOLD),F(KOLD,L),F(K,LOLD)) depend(out: E(K,L),F(K,L))  default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task
                END IF

                ! in current row
                IF ( KH  .LT. M ) THEN
                    JH = KH+1
                    DO WHILE ( JH .LE. M )
                        J = JH
                        JH = MIN(M, J + MB - 1)
                        JB = JH - J + 1
                        IF ( JH .LT. M ) THEN
                            IF ( A(JH+1,JH) .NE. ZERO )  THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF
                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:E(K,L)) depend(inout:E(J,L))

                        CALL DGEMM("T", "N", JB, LB, KB, -ONE, A(K,J), LDA, E(K,L), LDE, ONE, E(J, L), LDE)
                        !$omp end task
                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:E(K,L)) depend(inout:F(J,L))

                        CALL DGEMM("T", "N", JB, LB, KB, -ONE, C(K,J), LDC, E(K,L), LDE, ONE, F(J, L), LDF)
                        !$omp end task

                        JH = JH + 1
                    END DO
                END IF

                IF ( LH .LT. N ) THEN
                    JH = LH+1
                    DO WHILE ( JH .LE. N )
                        J = JH
                        JH = MIN(N, J + NB - 1)
                        JB = JH - J + 1
                        IF ( JH .LT. N ) THEN
                            IF ( B(JH+1,JH) .NE. ZERO .OR. D(JH+1,JH) .NE. ZERO )  THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF
                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:F(K,L)) depend(inout:E(K,J))

                        CALL DGEMM("N", "N", KB, JB, LB, -ONE*SGN1,  F(K,L), LDF, B(L,J), LDB, ONE, E(K,J), LDE)
                        !$omp end task

                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:F(K,L)) depend(inout:F(K,J))

                        CALL DGEMM("N", "N", KB, JB, LB, -ONE*SGN2, F(K,L), LDF, D(L,J), LDD, ONE, F(K,J), LDF)
                        !$omp end task

                        JH = JH + 1
                    END DO
                END IF

                KOLD = K
                KH = KH + 1

            END DO


            LOLD = L
            LH = LH + 1
        END DO
        !$omp end master
        !$omp taskwait
        !$omp end parallel


    ELSE IF ( BTRANSA .AND. BTRANSB) THEN
        ! (T, T)
        !$omp parallel default(shared)
        !$omp master

        KOLD = 1
        LOLD = 1
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
                IF  (B(L,L-1) .NE. ZERO .OR. D(L,L-1) .NE. ZERO)  THEN
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
                KH = MIN(M, K + MB - 1)
                KB = KH - K + 1
                IF ( KH .LT. M ) THEN
                    IF ( A(KH+1, KH) .NE. ZERO ) THEN
                        KB = KB - 1
                        KH = KH - 1
                    END IF
                END IF

                IF ( K .EQ. IONE .AND. LH .EQ. N ) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: F(K,L),E(K,L)) default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task
                ELSE IF (K .EQ. IONE .AND. LH .LT. N ) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: F(K,L),E(K,L)) depend(in: F(K,LOLD),E(K,LOLD)) default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)
                    !$omp end task

                ELSE IF (K .GT. IONE .AND. LH .EQ. N ) THEN
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(out: E(K,L),F(K,L)) depend(in: E(KOLD,L),F(KOLD,L))  default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)
                    !$omp end task

                ELSE
                    !$omp task firstprivate(K,KH,KB,L,LH,LB,SCAL,INFO1) &
                    !$omp& depend(in: E(KOLD,L), E(K,LOLD),F(KOLD,L),F(K,LOLD)) depend(out: E(K,L),F(K,L))  default(shared)
                    CALL IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N, K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
                        & E, LDE, F, LDF, SCAL,  INFO1)

                    !$omp end task

                END IF



                IF ( KH .LT. M  ) THEN
                    JH = KH+1
                    DO WHILE ( JH .LE. M )
                        J = JH
                        JH = MIN(M, J + MB - 1)
                        JB = JH - J + 1
                        IF ( JH .LT. M ) THEN
                            IF (( A(JH+1,JH) .NE. ZERO ) ) THEN
                                JB = JB - 1
                                JH = JH - 1
                            END IF
                        END IF
                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:E(K,L)) depend(inout:E(J,L))
                        CALL DGEMM('T', 'N', JB, LB, KB, -ONE, A(K, J), LDA, E(K,L), LDE, ONE, E(J,L),LDE)
                        !$omp end task

                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:E(K,L)) depend(inout:F(J,L))
                        CALL DGEMM('T', 'N', JB, LB, KB, -ONE, C(K, J), LDC, E(K,L), LDE, ONE, F(J,L),LDF)
                        !$omp end task

                        JH = JH +1
                    END DO
                END IF

                ! Update January 2021
                IF ( L .GT. 1) THEN
                    J = L - 1
                    DO WHILE (J .GT. 0)
                        JH = J
                        J  = MAX(1, JH - NB + 1)
                        JB = JH - J + 1
                        IF ( J .GT. 1 ) THEN
                            IF ((B(J,J-1) .NE. ZERO .OR. D(J,J-1).NE.ZERO) ) THEN
                                JB = JB - 1
                                J = J + 1
                            END IF
                        END IF
                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:F(K,L)) depend(inout:E(K,J))
                        CALL DGEMM('N', 'T', KB, JB, LB, -ONE*SGN1, F(K,L), LDF, B(J, L), LDB, ONE, E(K,J), LDE)
                        !$omp end task
                        !$omp task firstprivate(K,L,J,KB,LB,JB) default(shared) &
                        !$omp& depend(in:F(K,L)) depend(inout:F(K,J))
                        CALL DGEMM('N', 'T', KB, JB, LB, -ONE*SGN2, F(K,L), LDF, D(J, L), LDD, ONE, F(K,J), LDF)
                        !$omp end task

                        J = J - 1
                    END DO
                END IF

                KOLD = 1
                KH = KH + 1

            END DO


            LOLD = 1
            L = L - 1
        END DO
        !$omp end master
        !$omp taskwait
        !$omp end parallel


    END IF
    RETURN

END SUBROUTINE


SUBROUTINE IDLA_TGCSYLV_BLOCK(TRANSA, TRANSB, SGN1, SGN2, M, N ,K, L, KB, LB, A, LDA, B, LDB, C, LDC, D, LDD, &
        & E, LDE, F, LDF, SCAL, INFO1)
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
    IMPLICIT NONE
    CHARACTER TRANSA(1), TRANSB(1)
    DOUBLE PRECISION SGN1, SGN2, SCAL
    INTEGER LDA, LDB, LDC, LDD, LDE, LDF, M, N, K,L,KB, LB, INFO1
    DOUBLE PRECISION A(LDA,*), B(LDB, *), C(LDC,*), D(LDD,*), E(LDE, *), F(LDF,*)
    INTEGER MAXNB,ISOLVER
    DOUBLE PRECISION WORK(1)

    EXTERNAL DLA_TGCSYLV_L2
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_128
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_32
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_64
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_96
    EXTERNAL DLA_TGCSYLV_L2_REORDER
    EXTERNAL DLA_TGCSYLV_L2_UNOPT
    EXTERNAL DLA_TGCSYLV_RECURSIVE


    MAXNB = MAX(KB, LB)
    ISOLVER = TGCSYLV_ISOLVER()


    SELECT CASE(ISOLVER)
    CASE(1)
        IF ( MAXNB .LE. 32 ) THEN
            CALL DLA_TGCSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
        ELSE IF ( MAXNB .LE. 64) THEN
            CALL DLA_TGCSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
        ELSE IF ( MAXNB .LE. 96) THEN
            CALL DLA_TGCSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
        ELSE IF ( MAXNB .LE. 128) THEN
            CALL DLA_TGCSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
        ELSE
            CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
        END IF
    CASE(2)
        IF ( M .GT. 128 .OR. N .GT. 128) THEN
            CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
        ELSE
            CALL DLA_TGCSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
                & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
        END IF
    CASE(3)
        CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
            & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
    CASE(4)
        CALL DLA_TGCSYLV_L2(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
            & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
    CASE(5)
        CALL DLA_TGCSYLV_L2_UNOPT(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
            & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
    CASE(6)
        CALL DLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, KB, LB,  A(K,K), LDA, B(L,L), LDB, &
            & C(K,K), LDC, D(L,L), LDD,  E(K,L), LDE, F(K,L), LDF, SCAL, WORK(1), INFO1)
    END SELECT
    IF ( INFO1 .LT. 0 ) THEN
        INFO1 = 100
    END IF

END SUBROUTINE






