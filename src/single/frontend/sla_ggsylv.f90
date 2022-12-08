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

!> \brief Frontend for the solution of Generalized Sylvester Equations.
!
!  Definition:
!  ===========
!
!   SUBROUTINE SLA_GGSYLV(FACTA, FACTB, TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD,
!                            QA, LDQA, ZA, LDZA, QB, LDQB, ZB, LDZB,
!                            X, LDX, SCALE, WORK, LDWORK, INFO)
!           CHARACTER FACTA(1), FACTB(1), TRANSA(1), TRANSB(1)
!           INTEGER M , LDA, LDB, LDC, LDD, LDQA, LDQB, LDZA, LDZB, LDX, INFO, LDWORK
!           REAL SCALE, SGN
!           REAL A(LDA,*), B(LDB, *), C(LDC, *), D(LDD,*)
!           REAL QA(LDQA, *), ZA(LDZA, *), QB(LDQB, *), ZB(LDZB, *), X(LDX, *), WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> SLA_GGSYLV solves a generalized Sylvester equation of the following forms
!>
!>    op1(A) * X * op2(B) + op1(C) * X * op2(D) = SCALE * Y                              (1)
!>
!> or
!>
!>    op1(A) * X * op2(B) - op1(C) * X * op2(D) = SCALE * Y                              (2)
!>
!> where (A,C) is a M-by-M matrix pencil and (B,D) is a N-by-N matrix pencil.
!> The right hand side Y and the solution X M-by-N matrices. The pencils (A,C)
!> and (B,D) can be either given as general unreduced matrices or in terms of
!> their generalized Schur decomposition. If they are given as general matrices
!> their generalized Schur decomposition will be computed.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] FACTA
!> \verbatim
!>          FACTA is CHARACTER
!>          Specifies how the matrix pencil (A,C) is given.
!>          == 'N':  The matrix pencil (A,C) is given as a general matrices and its Schur decomposition
!>                  A = QA*S*ZA**T, C = QA*R*ZA**T will be computed.
!>          == 'F':  The matrix pencil (A,C) is already in generalized Schur form and S, R, QA, and ZA
!>                  are given.
!> \endverbatim
!
!> \param[in] FACTB
!> \verbatim
!>          FACTB is CHARACTER
!>          Specifies how the matrix pencil (B,D) is given.
!>          == 'N':  The matrix pencil (B,D) is given as a general matrices and its Schur decomposition
!>                  B = QB*U*ZB**T, D = QB*V*ZB**T will be computed.
!>          == 'F':  The matrix pencil (B,D) is already in generalized Schur form and U, V, QB, and ZB
!>                  are given.
!> \endverbatim
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER
!>          Specifies the form of the system of equations with respect to A and C :
!>          == 'N':  op1(A) = A
!>          == 'T':  op1(A) = A**T
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER
!>          Specifies the form of the system of equations with respect to B and D:
!>          == 'N':  op2(B) = B,
!>          == 'T':  op2(B) = B**T
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
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,M)
!>          If FACT == "N", the matrix A is a general matrix and it is overwritten with the
!>          (quasi-) upper triangular factor S of the Schur decomposition of (A,C).
!>          If FACT == "F", the matrix A contains its (quasi-) upper triangular matrix S of
!>          the Schur decomposition of (A,C).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
!>          If FACT == "N",  the matrix B is a general matrix and it is overwritten with the
!>          (quasi-) upper triangular factor U of the Schur decomposition of (B,D).
!>          If FACT == "F", the matrix B contains its (quasi-) upper triangular matrix U of
!>          the Schur decomposition of (B,D).
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is REAL array, dimension (LDC,M)
!>          If FACT == "N", the matrix C is a general matrix and it is overwritten with the
!>          upper triangular factor R of the Schur decomposition of (A,C).
!>          If FACT == "F", the matrix C contains its upper triangular matrix R of
!>          the Schur decomposition of (A,C).
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,M).
!> \endverbatim
!>
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (LDD,N)
!>          If FACT == "N",  the matrix D is a general matrix and it is overwritten with the
!>          upper triangular factor V of the Schur decomposition of (B,D).
!>          If FACT == "F", the matrix D contains its upper triangular matrix V of
!>          the Schur decomposition of (B,D).
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array D.  LDD >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] QA
!> \verbatim
!>          QA is REAL array, dimension (LDQA,M)
!>          If FACT == "N", the matrix QA is an empty M-by-M matrix on input and contains the
!>          left Schur vectors of (A,C) on output.
!>          If FACT == "F", the matrix QA contains the left Schur vectors of (A,C).
!> \endverbatim
!>
!> \param[in] LDQA
!> \verbatim
!>          LDQA is INTEGER
!>          The leading dimension of the array QA.  LDQA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] ZA
!> \verbatim
!>          ZA is REAL array, dimension (LDZA,M)
!>          If FACT == "N", the matrix ZA is an empty M-by-M matrix on input and contains the
!>          right Schur vectors of (A,C) on output.
!>          If FACT == "F", the matrix ZA contains the right Schur vectors of (A,C).
!> \endverbatim
!>
!> \param[in] LDZA
!> \verbatim
!>          LDZA is INTEGER
!>          The leading dimension of the array ZA.  LDZA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] QB
!> \verbatim
!>          QB is REAL array, dimension (LDQB,N)
!>          If FACT == "N", the matrix QB is an empty N-by-N matrix on input and contains the
!>          left Schur vectors of (B,D) on output.
!>          If FACT == "F", the matrix QB contains the left Schur vectors of (B,D).
!> \endverbatim
!>
!> \param[in] LDQB
!> \verbatim
!>          LDQB is INTEGER
!>          The leading dimension of the array QB.  LDQB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] ZB
!> \verbatim
!>          ZB is REAL array, dimension (LDZB,N)
!>          If FACT == "N", the matrix ZB is an empty N-by-N matrix on input and contains the
!>          right Schur vectors of (B,D) on output.
!>          If FACT == "F", the matrix ZB contains the right Schur vectors of (B,D).
!> \endverbatim
!>
!> \param[in] LDZB
!> \verbatim
!>          LDZB is INTEGER
!>          The leading dimension of the array ZB.  LDZB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is REAL array, dimension (LDX,N)
!>          On input, the matrix X contains the right hand side Y.
!>          On output, the matrix X contains the solution of Equation (1) or (2)
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
!> \param[in,out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LDWORK))
!>          Workspace for the algorithm. The optmimal workspace is given either by \ref mepack_memory_frontend
!>          or a previous call to the this routine with LDWORK === -1.
!> \endverbatim
!>
!> \param[in,out] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          Size of the workspace for the algorithm. This can be determined by a call \ref mepack_memory_frontend .
!>          Alternatively, if LDWORK == -1 on input, the subroutine will return the required size of the workspace in LDWORK
!>          without performing any computations.
!> \endverbatim
!>
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          == 0:  successful exit
!>          = 1:  SGGES failed
!>          = 2:  SLA_SORT_GEV failed
!>          = 3:  Inner solver failed
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!>
!
!> \see SLA_TGSYLV_DAG
!> \see SLA_TGSYLV_L3_COLWISE
!> \see SLA_TGSYLV_L3_2S
!> \see SLA_TGSYLV_L2_REORDER
!> \see SLA_TGSYLV_L2
!> \see SLA_TGSYLV_L2_COLWISE
!> \see SLA_TGSYLV_L2_LOCAL_COPY_32
!> \see SLA_TGSYLV_L2_LOCAL_COPY_64
!> \see SLA_TGSYLV_L2_LOCAL_COPY_96
!> \see SLA_TGSYLV_L2_LOCAL_COPY_128
!> \see SLA_TGSYLV_L2_LOCAL_COPY
!> \see SLA_TGSYLV_GARDINER_LAUB
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Dezember 2022
!> \ingroup sglggsylv
!
SUBROUTINE SLA_GGSYLV(FACTA, FACTB, TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
        & QA, LDQA, ZA, LDZA, QB, LDQB, ZB, LDZB, &
        & X, LDX, SCALE, WORK, LDWORK, INFO)
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_SORTEV
    USE MEPACK_OPTIONS_FRONTEND_SOLVER
    IMPLICIT NONE

    ! Arguments
    CHARACTER FACTA(1), FACTB(1), TRANSA(1), TRANSB(1)
    INTEGER M, N , LDA, LDB, LDC, LDD, LDQA, LDQB, LDZA, LDZB, LDX, INFO, LDWORK
    REAL SCALE, SGN
    REAL A(LDA,*), B(LDB, *), C(LDC, *), D(LDD,*)
    REAL QA(LDQA, *), QB(LDQB, *), ZA(LDZA, *), ZB(LDZB, *), X(LDX, *), WORK(*)



    ! Local Variables
    INTEGER ININFO
    INTEGER LDWORKI, IINFO
    LOGICAL BTRANSA, BTRANSB
    LOGICAL BFACTA, BFACTB, DUMMY
    INTEGER SDIM, MB, ISOLVER, NB
    INTEGER SORTEV
    INTEGER SOLVER
    INTEGER MAXNB
    INTEGER BIGNB

    REAL ONE, ZERO
    PARAMETER(ONE = 1.0, ZERO = 0.0)


    ! External Functions

    EXTERNAL SLA_TRANSFORM_GENERAL
    EXTERNAL SLA_TGSYLV_DAG
    EXTERNAL SLA_TGSYLV_L3_COLWISE
    EXTERNAL SLA_TGSYLV_L3_2S
    EXTERNAL SLA_TGSYLV_L2
    EXTERNAL SLA_TGSYLV_L2_REORDER
    EXTERNAL SLA_TGSYLV_L2_COLWISE
    EXTERNAL SLA_TGSYLV_RECURSIVE
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_32
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_64
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_96
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY_128
    EXTERNAL SLA_TGSYLV_L2_LOCAL_COPY
    EXTERNAL SLA_TGSYLV_GARDINER_LAUB

    EXTERNAL SGGES
    EXTERNAL SLA_SORT_GEV
    EXTERNAL LSAME
    EXTERNAL XERROR_HANDLER
    LOGICAL LSAME


    ! Check Input
    BTRANSA = LSAME(TRANSA, 'N')
    BTRANSB = LSAME(TRANSB, 'N')
    BFACTA  = LSAME(FACTA, 'F')
    BFACTB  = LSAME(FACTB, 'F')

    MB = TGSYLV_BLOCKSIZE_MB(M,N)
    NB = TGSYLV_BLOCKSIZE_NB(M,N)
    ISOLVER = TGSYLV_ISOLVER()
    SORTEV = SORTEV_ENABLED()
    SOLVER = TGSYLV_FRONTEND_SOLVER()
    BIGNB = TGSYLV_BLOCKSIZE_2STAGE(M,N)

    MAXNB = MAX(M,N)

    ININFO = INFO
    INFO = 0
    IF ( .NOT. BFACTA .AND. .NOT. LSAME(FACTA, 'N')) THEN
        INFO = -1
    ELSE IF ( .NOT. BFACTB .AND. .NOT. LSAME(FACTB, 'N')) THEN
        INFO = -2
    ELSE IF ( .NOT. BTRANSA .AND. .NOT. LSAME(TRANSA, 'T')) THEN
        INFO = -3
    ELSE IF ( .NOT. BTRANSB .AND. .NOT. LSAME(TRANSB, 'T')) THEN
        INFO = -4
    ELSE IF ( SGN .NE. -ONE .AND. SGN .NE. ONE ) THEN
        INFO = -5
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -6
    ELSE IF ( N .LT. 0 ) THEN
        INFO = -7
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -9
    ELSE IF ( LDB .LT. MAX(1, N)) THEN
        INFO = -11
    ELSE IF ( LDC .LT. MAX(1, M)) THEN
        INFO = -13
    ELSE IF ( LDD .LT. MAX(1, N)) THEN
        INFO = -15
    ELSE IF ( LDQA .LT. MAX(1, M)) THEN
        INFO = -17
    ELSE IF ( LDZA .LT. MAX(1, M)) THEN
        INFO = -19
    ELSE IF ( LDQB .LT. MAX(1, N)) THEN
        INFO = -21
    ELSE IF ( LDZB .LT. MAX(1, N)) THEN
        INFO = -23
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -25
    ELSE IF ( MB .LT. 2 ) THEN
        INFO = -30
    ELSE IF ( NB .LT. 2 ) THEN
        INFO = -31
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32
    ELSE IF ( SOLVER .LT. 1 .OR. SOLVER .GT. 6 .OR. ( SOLVER.EQ.FRONTEND_SOLVER_2STAGE .AND. BIGNB .LT. MIN(MB, NB)) ) THEN
        INFO = -33
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_GGSYLV', -INFO)
        RETURN
    END IF

    !
    ! Compute the required workspace
    !
    LDWORKI = 1
    IF ( .NOT. BFACTA) THEN
        LDWORKI = MAX(LDWORKI, 3*M + 8*M+16)
    END IF
    IF ( .NOT. BFACTB) THEN
        LDWORKI = MAX(LDWORKI, 3*N + 8*N+16)
    END IF
    IINFO = -1
    SELECT CASE(SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL SLA_TGSYLV_L3_COLWISE(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_LEVEL2)
            SELECT CASE(ISOLVER)
            CASE(1)
                IF ( MAXNB .LE. 32 ) THEN
                    CALL SLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                ELSE IF ( MAXNB .LE. 64) THEN
                    CALL SLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                ELSE IF ( MAXNB .LE. 96) THEN
                    CALL SLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                ELSE IF ( MAXNB .LE. 128) THEN
                    CALL SLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                ELSE
                    CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                END IF
            CASE(2)
                IF ( M .GT. 128 .OR. N .GT. 128) THEN
                    CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                ELSE
                    CALL SLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                END IF
            CASE(3)
                CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
            CASE(4)
                CALL SLA_TGSYLV_L2(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
            CASE(5)
                CALL SLA_TGSYLV_L2_COLWISE(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
            CASE(6)
                CALL SLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)


            END SELECT
        CASE (FRONTEND_SOLVER_DAG)
            CALL SLA_TGSYLV_DAG(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_2STAGE)
            CALL SLA_TGSYLV_L3_2S(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_RECURSIVE)
            CALL SLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_GARDINER_LAUB)
            CALL SLA_TGSYLV_GARDINER_LAUB(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, IINFO)
    END SELECT

    IF ( IINFO .LT. 0 ) THEN
        INFO = -28
        CALL XERROR_HANDLER("SLA_GGSYLV", -INFO)
        RETURN
    END IF
    LDWORKI = MAX(LDWORKI, IINFO, M*N)

    IF(LDWORK .EQ. -1) THEN
        LDWORK = LDWORKI
        RETURN
    ELSE IF (LDWORKI .GT. LDWORK) THEN
        INFO = -28
        CALL XERROR_HANDLER("SLA_GGSYLV", -INFO)
        RETURN
    END IF


    !
    ! Quick Return
    !
    INFO = 0
    IF ( M .EQ. 0 ) RETURN
    IF ( N .EQ. 0 ) RETURN

    !
    ! Factorize A
    !
    IF (.NOT. BFACTA) THEN
        CALL SGGES("Vectors", "Vectors", "NoSort", DUMMY, M, A, LDA, C, LDC, SDIM, WORK(1), WORK(M+1), WORK(2*M+1), &
            & QA, LDQA,  ZA, LDZA,  WORK(3*M+1), LDWORKI-3*M, DUMMY, IINFO)
        IF ( IINFO .NE. 0 ) THEN
            INFO = 1
            RETURN
        END IF

        IF ( SORTEV .NE. 0 ) THEN
            CALL SLA_SORT_GEV(M, A, LDA, C, LDC, QA, LDQA, ZA, LDZA, MB, WORK, LDWORKI, IINFO)
            IF ( IINFO .NE. 0 ) THEN
                INFO = 2
                RETURN
            END IF
        END IF
    END IF

    !
    ! Factorize B
    !
    IF (.NOT. BFACTB) THEN
        CALL SGGES("Vectors", "Vectors", "NoSort", DUMMY, N, B, LDB, D, LDD, SDIM, WORK(1), WORK(N+1), WORK(2*N+1), &
            & QB, LDQB,  ZB, LDZB,  WORK(3*N+1), LDWORKI-3*N, DUMMY, IINFO)
        IF ( IINFO .NE. 0 ) THEN
            INFO = 1
            RETURN
        END IF

        IF ( SORTEV .NE. 0 ) THEN
            CALL SLA_SORT_GEV(N, B, LDB, D, LDD, QB, LDQB, ZB, LDZB, NB, WORK, LDWORKI, IINFO)
            IF ( IINFO .NE. 0 ) THEN
                INFO = 2
                RETURN
            END IF
        END IF

    END IF

    ! Transform Y
    !
    ! TRANS A      TRANS B         Y <-
    !   N            N             QA**T*Y*ZB
    !   N            T             QA**T*Y*QB
    !   T            N             ZA**T*Y*ZB
    !   T            T             ZA**T*Y*QB
    IF ( BTRANSA .AND. BTRANSB ) THEN   !(N,N)
        CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, X, LDX, QA, LDQA, ZB, LDZB, WORK)
    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB ) THEN   !(N,T)
        CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, X, LDX, QA, LDQA, QB, LDQB, WORK)
    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN   !(T,N)
        CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, X, LDX, ZA, LDZA, ZB, LDZB, WORK)
    ELSE !(T,T)
        CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, X, LDX, ZA, LDZA, QB, LDQB, WORK)
    END IF

    ! Solve
    IINFO = 0
    SELECT CASE(SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL SLA_TGSYLV_L3_COLWISE(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_LEVEL2)
            SELECT CASE(ISOLVER)
            CASE(1)
                IF ( MAXNB .LE. 32 ) THEN
                    CALL SLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                ELSE IF ( MAXNB .LE. 64) THEN
                    CALL SLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                ELSE IF ( MAXNB .LE. 96) THEN
                    CALL SLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                ELSE IF ( MAXNB .LE. 128) THEN
                    CALL SLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                ELSE
                    CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                END IF
            CASE(2)
                IF ( M .GT. 128 .OR. N .GT. 128) THEN
                    CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                ELSE
                    CALL SLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                END IF
            CASE(3)
                CALL SLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
            CASE(4)
                CALL SLA_TGSYLV_L2(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
            CASE(5)
                CALL SLA_TGSYLV_L2_COLWISE(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
            CASE(6)
                CALL SLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                    & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)


            END SELECT
        CASE (FRONTEND_SOLVER_DAG)
            CALL SLA_TGSYLV_DAG(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_2STAGE)
            CALL SLA_TGSYLV_L3_2S(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_RECURSIVE)
            CALL SLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_GARDINER_LAUB)
            CALL SLA_TGSYLV_GARDINER_LAUB(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, IINFO)
    END SELECT



    IF ( IINFO .NE. 0) THEN
        INFO = 3
        RETURN
    END IF

    ! Transform X
    !
    ! TRANS A      TRANS B         X <-
    !   N            N             ZA*Y*QB**T
    !   N            T             ZA*Y*ZB**T
    !   T            N             QA*Y*QB**T
    !   T            T             QA*Y*ZB**T
    IF ( BTRANSA .AND. BTRANSB ) THEN   !(N,N)
        CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, X, LDX, ZA, LDZA, QB, LDQB, WORK)
    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB ) THEN   !(N,T)
        CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, X, LDX, ZA, LDZA, ZB, LDZB, WORK)
    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN   !(T,N)
        CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, X, LDX, QA, LDQA, QB, LDQB, WORK)
    ELSE !(T,T)
        CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, X, LDX, QA, LDQA, ZB, LDZB, WORK)
    END IF

    RETURN
END SUBROUTINE

